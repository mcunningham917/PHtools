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

minBenchLength; % In catchments, consecutive number of catchments with shared PH mode 
spillOverElevations = phbBandThickness; % Elevation above and below PH modal elevation
areaThreshPixelNum = Ac; % In pixels, must be same as in PH
cuSumThresh = .02; % For PH bench separation
pixelLength = 30;
nanFlag = -32768;

%% 
for count = supercatchmentNum

    
    streamSupercatchment = count;

    supercatchmentFilePath = fullfile(phDataFilePath,groupArea,'Supercatchments');
    supercatchmentFileName = ['Supercatchment', num2str(streamSupercatchment)];
    supercatchmentDemName = [groupArea, 'Supercatchment', num2str(streamSupercatchment),'.tif']
    deltaHOutputFileName = [groupArea, 'Supercatchment', num2str(count), '_allPHBs_deltaH.tif'];
    hBenchOutputFileName = [groupArea, 'Supercatchment', num2str(count), '_allPHBs_hBench.tif'];

    %Input list of (outlet, mode) pairs
    subcatchmentFolderName = [num2str(phStepLength),'mStep'];
    progressivePourPointSubcatchmentFilePath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'Subcatchments', subcatchmentFolderName, num2str(supercatchmentFileName));
    
    %Output file path for PHB layer
    allSupercatchmentPHBfilePath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'PHBs', 'Cusum02_BenchLength3Steps','Maps');
    allSupercatchmentPHBTablePath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'PHBs', 'Cusum02_BenchLength3Steps','Tables');
    supercatchmentTableName = ['Supercatchment', num2str(count), '_allOutletModePairs.txt'];
   %Ouput file path for polygon layer
    %allSupercatchmentPHBforPolygons = fullfile(phAnalysisFilePath,'Cusum02_BenchLength3Steps', 'AllSupercatchments_ForPolygons');

    mkdir(allSupercatchmentPHBfilePath);
    %mkdir(allSupercatchmentPHBforPolygons);


    
    %% Delineate channel network for supercatchment
    
    [supercatchmentDemArray, supercatchmentDemGrid, supercatchmentDemInfo, supercatchmentGeospatialReferenceArray] =...
        loadDEM(supercatchmentFilePath, supercatchmentDemName, nanFlag);
    
    supercatchmentDEMArrayForPHB = supercatchmentDemArray;
    supercatchmentDEMArrayForPHBIX = find(~isnan(supercatchmentDemArray));
    supercatchmentDEMArrayForPHB(supercatchmentDEMArrayForPHBIX)=NaN;
    supercatchmentDEMArrayForDeltaH=supercatchmentDEMArrayForPHB;
    

    demSinksFilled = fillsinks(supercatchmentDemGrid);
    flowDirectionObj = FLOWobj(demSinksFilled);
    flowAccumulationStruct  = flowacc(flowDirectionObj);
    streamNodeThreshNum=areaThreshPixelNum; 
    streamNodeThreshold = flowAccumulationStruct>streamNodeThreshNum;
    streamNetworkStreamObj = STREAMobj(flowDirectionObj,streamNodeThreshold);
    streamNetworkStruct=STREAMobj2mapstruct(streamNetworkStreamObj);
    

    phbOutletArray = dlmread(fullfile(allSupercatchmentPHBTablePath, supercatchmentTableName),'\t',2);

    modeOutletArray = phbOutletArray(:,[1,2]);
    modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
    zeroIndices = find(modeOutletArray(:,2)==0);
    modeOutletArray(zeroIndices,:)=[];
   
    [~, ia, ic]  = unique(modeOutletArray, 'rows');
    modeOutletArrayUnique = phbOutletArray(ia,[1,2]);
    phbOutletArrayUnique = phbOutletArray(ia,:)

    
%% Find PHBs for each first order tributary in the supercatchment

    % Bring in (outlet, mode) pair for each PH chain
    
   
        for benchListNum = 1:length(phbOutletArrayUnique);
       

                streamNum = phbOutletArrayUnique(benchListNum,4)
        
                hypsoPeakElevation = round(min(phbOutletArrayUnique(benchListNum,2)));
                benchOutletElevation = round(min(phbOutletArrayUnique(benchListNum,1)));
                
                deltaH = hypsoPeakElevation-benchOutletElevation;
       
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
                        supercatchmentDEMArrayForDeltaH(hypsoPeakElevationIndices)=deltaH;

                    else
            
                    % If the catchment is above the defined channel, an
                    % the path of maximum flow accumulation needs to be
                    % tracked to the main divide.
                    [streamHeadElevation, streamHeadIndex] = max(profileZ);
                     streamHeadX = profileX(streamHeadIndex);
                     streamHeadY = profileY(streamHeadIndex);
            
                    highSubcatchmentGrid = ExtractSubcatchmentGRID(supercatchmentDemArray,...
                        flowDirectionObj,streamHeadX, streamHeadY);
                    
                        subcatchmentDemArray = highSubcatchmentGrid.Z;
                        subcatchmentDEMVec = subcatchmentDemArray(~isnan(subcatchmentDemArray));
                        superDemVec = supercatchmentDemArray(~isnan(supercatchmentDemArray));
               
                        if(((length(subcatchmentDEMVec)/length(superDemVec))<0.25))
            
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
                    supercatchmentDEMArrayForDeltaH(hypsoPeakElevationIndices)=deltaH;
                    
                        end
                   

                        
                        
                    end
                    
                    %supercatchmentDEMArrayForPHB(supercatchmentDEMArrayForPHB==1) = NaN;  
                    %supercatchmentDEMArrayForPHB = isfloat(supercatchmentDEMArrayForPHB);
                    
                    hBenchfullOutputFileForMaster = fullfile(allSupercatchmentPHBfilePath, hBenchOutputFileName);
                    deltaHfullOutputFileForMaster = fullfile(allSupercatchmentPHBfilePath, deltaHOutputFileName);
                    
                    SubcatchmentWrite(supercatchmentDEMArrayForPHB, supercatchmentGeospatialReferenceArray,supercatchmentDemInfo,...
                                nanFlag, hBenchfullOutputFileForMaster);
                   SubcatchmentWrite(supercatchmentDEMArrayForDeltaH, supercatchmentGeospatialReferenceArray,supercatchmentDemInfo,...
                                nanFlag, deltaHfullOutputFileForMaster); 

            end
        end






