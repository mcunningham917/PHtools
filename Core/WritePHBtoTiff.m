%% Generate list of PHBs and output PHB layers
% 
% Description
% 
% Find steps in PH output, and writes a geotiff of elevation bands around
% each PH mode. Each band interatively "fills in" a blank copy of the supercatment.
% Final output is a geotiff with each PHB band.
% 
%% Set variables


Defaults;
addpath(topoToolboxFilePath); 

minBenchLength = 3; % In catchments, consecutive number of catchments with shared PH mode 
spillOverElevations = 25; % Elevation above and below PH modal elevation
areaThreshPixelNum = 2777; % In pixels, must be same as in PH
cuSumThresh = .02; % For PH bench separation
pixelLength = 30;
nanFlag = -32768;

%% 
for count = supercatchmentNum

    
    streamSupercatchment = count;

    supercatchmentFilePath = fullfile(phDataFilePath,groupArea,'Supercatchments');
    supercatchmentFileName = ['Supercatchment', num2str(streamSupercatchment)];
    supercatchmentDemName = [groupArea, 'Supercatchment', num2str(streamSupercatchment),'.tif']
    PHBOutputFileName = [groupArea, 'Supercatchment', num2str(count), '_allPHBs.tif'];

    %Input list of (outlet, mode) pairs
    progressivePourPointSubcatchmentFilePath = fullfile(phAnalysisFilePath,groupArea,'Subcatchments','25mStep', num2str(supercatchmentFileName));
    
    %Output file path for PHB layer
    allSupercatchmentPHBfilePath = fullfile(phAnalysisFilePath,groupArea,'PHBs', 'Cusum02_BenchLength3Steps','AllSupercatchmentsTiffs');
    
   %Ouput file path for polygon layer
    %allSupercatchmentPHBforPolygons = fullfile(phAnalysisFilePath,'Cusum02_BenchLength3Steps', 'AllSupercatchments_ForPolygons');

    mkdir(allSupercatchmentPHBfilePath);
    %mkdir(allSupercatchmentPHBforPolygons);


    
    %% Delineate channel network for supercatchment
    
    [supercatchmentDemArray, supercatchmentDemGrid, supercatchmentDemInfo, supercatchmentGeospatialReferenceArray] =...
        loadDEM(supercatchmentFilePath, supercatchmentDemName, nanFlag);
    
    supercatchmentDEMArrayForPHB = supercatchmentDemArray;
    supercatchmentDEMArrayForPHBIX = find(~isnan(supercatchmentDemArray));
    supercatchmentDEMArrayForPHB(supercatchmentDEMArrayForPHBIX)=1;
    

    demSinksFilled = fillsinks(supercatchmentDemGrid);
    flowDirectionObj = FLOWobj(demSinksFilled);
    flowAccumulationStruct  = flowacc(flowDirectionObj);
    streamNodeThreshNum=areaThreshPixelNum; 
    streamNodeThreshold = flowAccumulationStruct>streamNodeThreshNum;
    streamNetworkStreamObj = STREAMobj(flowDirectionObj,streamNodeThreshold);
    streamNetworkStruct=STREAMobj2mapstruct(streamNetworkStreamObj);
    

    %Get list of chains
    clear firstOrderStreamPath firstOrderStreamList
    
    firstOrderStreamPath = dir(progressivePourPointSubcatchmentFilePath);

        for i = 1:length(firstOrderStreamPath)
    
            realFolder = strfind(firstOrderStreamPath(i).name,'Super')
    
            if(~isempty(realFolder))
                firstOrderStreamList(i) = i;
            end
            
        end
    
    firstOrderStreamList=firstOrderStreamList(firstOrderStreamList>0);


    
%% Find PHBs for each first order tributary in the supercatchment

    % Bring in (outlet, mode) pair for each PH chain
    
    for highTributaryNum = 1:length(firstOrderStreamList)
        clear allSubcatchmentDataArray hypsoPeakList outletList
        
        chainNum = firstOrderStreamList(highTributaryNum)
        chainNumName = fullfile(progressivePourPointSubcatchmentFilePath, firstOrderStreamPath(chainNum).name)
        allSubcatchmentDataArray= dlmread(chainNumName);
            
        streamNumFile = firstOrderStreamPath(chainNum).name;
        word1Location = strfind(streamNumFile, 'Num');
        word2Location = strfind(streamNumFile, '_Mode');
        inBetweenText = streamNumFile(word1Location+3:word2Location-1);
        streamNum = str2double(inBetweenText);

        [outletList, orderedList] = sort(allSubcatchmentDataArray(:,1), 'ascend');
        hypsoPeakList = allSubcatchmentDataArray(orderedList,2);
       

        % Find PHBs
        clear benchStruct benchEndPoints
        
        % Get endpoints of PH benches
        benchEndPoints = ProminentHypsoID(hypsoPeakList, cuSumThresh, minBenchLength);
    
        if(benchEndPoints>0)
    
            % Create a structure array of all PHBs
            benchStruct = PHBID(benchEndPoints, outletList, hypsoPeakList, minBenchLength);
    
        if(~isempty(benchStruct))
    
        % Use PHB structure array to pull out all subcatchments

            for benchListNum = 1: length(benchStruct)
        
                hypsoPeakElevation = round(min(benchStruct(benchListNum).HypsoPeaks));
                benchOutletElevation = round(min(benchStruct(benchListNum).BenchOutlets));
                maxBenchOutlet = max(benchStruct(benchListNum).BenchOutlets);
       
                benchOutletElevation = round(benchOutletElevation);
                targetPourPointElevation = benchOutletElevation;
                
    
                [profileX, profileY, profileZ,~, ~, ~] =...
                    longProfileGenerator(streamNetworkStruct, streamNetworkStreamObj,supercatchmentDemGrid, streamNum, streamNodeThreshNum);
        
                % First, check that the target pour point lies below the
                % artificial channel head
        
                    if(targetPourPointElevation< max(profileZ))
        
                        [~, closestPourPointIndex] = min(abs(profileZ - targetPourPointElevation));
 
                        subcatchmentOutlet = [profileX(closestPourPointIndex), profileY(closestPourPointIndex)]; 
                        subcatchmentOutletElevation = round(profileZ(closestPourPointIndex));

                        % Extract catchment above h_change for given PHB
                        finalDrainageBasinArray = ExtractSubcatchmentDEM(flowDirectionObj, supercatchmentDemArray,...
                            subcatchmentOutlet(1,1), subcatchmentOutlet(1,2));
         
                        finalDrainageBasinArrayArea = finalDrainageBasinArray(~isnan(finalDrainageBasinArray));
                       
                        % Determine which pixels will be labeled with PH
                        % modal elevation (specified elevation thresh above
                        % and below)
                        
                        hypsoPeakElevationIndices = find(finalDrainageBasinArray>=(hypsoPeakElevation-spillOverElevations) &...
                            finalDrainageBasinArray<=(hypsoPeakElevation+spillOverElevations));
                        
                        supercatchmentDEMArrayForPHB(hypsoPeakElevationIndices)=hypsoPeakElevation;

                    else
            
                    % If the catchment is above the defined channel, an
                    % the path of maximum flow accumulation needs to be
                    % tracked to the main divide.
                    [streamHeadElevation, streamHeadIndex] = max(profileZ);
                     streamHeadX = profileX(streamHeadIndex);
                     streamHeadY = profileY(streamHeadIndex);
            
                    highSubcatchmentGrid = ExtractSubcatchmentGRID(supercatchmentDemArray,...
                        flowDirectionObj,streamHeadX, streamHeadY);
            
                    highSubcatchmentFilled = fillsinks(highSubcatchmentGrid);
                    highFlowDirectionObj = FLOWobj(highSubcatchmentFilled);
                    highFlowAccumulationStruct  = flowacc(highFlowDirectionObj);
            
                    [aboveHeadZ, upstreamNodes]= TrackToDivide(highSubcatchmentGrid, highFlowDirectionObj, highFlowAccumulationStruct)
            
                    [~, closestPourPointIndex] = min(abs(aboveHeadZ - targetPourPointElevation));
                    subcatchmentOutletElevation = round(aboveHeadZ(closestPourPointIndex));
                    deltaH = hypsoPeakElevation-subcatchmentOutletElevation;
            
                    upstreamNodeArray = [aboveHeadZ, upstreamNodes]
                    outletNode = upstreamNodeArray(closestPourPointIndex,2);
                    outletElevation = upstreamNodeArray(closestPourPointIndex,1);
                    [subcatchmentOutletX, subcatchmentOutletY] = ind2coord(highSubcatchmentGrid, outletNode);
            
            
                    finalDrainageBasinArray = ExtractSubcatchmentDEM(flowDirectionObj, supercatchmentDemArray,...
                    subcatchmentOutletX, subcatchmentOutletY);
        
    
    
                    hypsoPeakElevationIndices = find(finalDrainageBasinArray>=(hypsoPeakElevation-spillOverElevations) & finalDrainageBasinArray<=(hypsoPeakElevation+spillOverElevations));
 
                        
                    supercatchmentDEMArrayForPHB(hypsoPeakElevationIndices)=hypsoPeakElevation;
                    
                    
                   

                        
                        
                    end
                    
                    %supercatchmentDEMArrayForPHB(supercatchmentDEMArrayForPHB==1) = NaN;  
                    %supercatchmentDEMArrayForPHB = isfloat(supercatchmentDEMArrayForPHB);
                    
                    fullOutputFileForMaster = fullfile(allSupercatchmentPHBfilePath, PHBOutputFileName);
                    
                    SubcatchmentWrite(supercatchmentDEMArrayForPHB, supercatchmentGeospatialReferenceArray,supercatchmentDemInfo,...
                                nanFlag, fullOutputFileForMaster); 

            end
        end
        end
    end
end





