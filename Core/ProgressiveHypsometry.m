%% Progressive hypsometry
%
% Description
% 
% Within a supercatchment, identify chains (channels), progressively step
% upstream along each channel, extract the catchment draining to each step
% elevation, and compute the modal elevation.
%
% Output
%
% A .txt file is written for each chain that includes (outlet, mode) pairs

%% Set variables

pourPointElevationStep = 25; % Step elevation height--vertical spacing of catchment elevations
pourPointElevationSensitivity = 10; % Variability around each step height
streamNodeThreshNum=Ac; %In pixels, area threshold required to define channel head (at 30 m pixel size, 2777 is 2.5 sq. km)
nanFlag = -32768;
   
supercatchmentFilePath =...
    fullfile(phDataFilePath,groupArea,'Supercatchments');

supercatchmentList = dir(supercatchmentFilePath)
    
highSubcatchmentDumpFile = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'HighSubcatchmentDumpTemp');

mkdir(highSubcatchmentDumpFile);
        
%% Progressive hypsometry 

% Loop through each supercatchment in file list

for count =supercatchmentNum
    
    clear streamProfileStruct
    
    
    supercatchmentFileName = [groupArea,'Supercatchment',num2str(count),'.tif'];
    supercatchmentFolderName = ['Supercatchment', num2str(count)];
    outputFilePath=fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'Subcatchments','25mStep', supercatchmentFolderName);
    mkdir(outputFilePath);
 
%% Bring in supercatchment DEM, write to GridObj, and delineate channel network

    [demArray, demGrid, demInfo, geospatialReferenceArray] =...
        loadDEM(supercatchmentFilePath, supercatchmentFileName, nanFlag);
   
    %Use topotoolbox to delineate stream network
    demSinksFilled = fillsinks(demGrid);
    flowDirectionObj = FLOWobj(demSinksFilled);
    flowAccumulationStruct  = flowacc(flowDirectionObj);
    streamNodeThreshold = flowAccumulationStruct>streamNodeThreshNum;
    streamNetworkStreamObj = STREAMobj(flowDirectionObj,streamNodeThreshold);
    streamNetworkStruct=STREAMobj2mapstruct(streamNetworkStreamObj);

%% Build structure array of x,y,z of each chain (channel)


for streamCount= 1:length(streamNetworkStruct)
    
    if(streamNetworkStruct(streamCount).streamorder==1)
        
        streamIndex = streamNetworkStruct(streamCount).IX;
        
        %Extract their profiles using longProfileGenerator
        %Use only raw Z values (not smoothed)
        [streamX, streamY, streamZ, ~, ~, ~] =...
            longProfileGenerator(streamNetworkStruct, streamNetworkStreamObj, demGrid, streamIndex, streamNodeThreshNum);
     
        %Build struct with X,Y,Z for each streamline
        streamProfileStruct(streamCount) = ...
            struct('XCoords',streamX,'YCoords',streamY, 'Elevation', streamZ, 'StreamNum', streamIndex);
    end
end


   
%% Progressive hypsometry core

% Iteratively extract catchments along progressively rising channel chain
% Loop nested in supercatchment loop

    streamCount = 1; % Reset counter for each supercatchment

%%

    for streamCount = 1:length(streamProfileStruct)
        clear targetPourPointElevationList remainingStepsList allSubcatchmentDataArray
    
    
        if(~isempty(streamProfileStruct(streamCount).StreamNum))
 
            % Find min and max of each channel chain
            streamNum = streamProfileStruct(streamCount).StreamNum;
            streamLineZ = streamProfileStruct(streamCount).Elevation;
            maxStreamElevation = max(streamLineZ); % Elevation of defined channel head
            minStreamElevation = min(streamLineZ);

            % Make a list of target pourpoints based on total height drop of 
            numPourPointSteps = (maxStreamElevation - minStreamElevation)/pourPointElevationStep;
            targetPourPointElevationList = (linspace(minStreamElevation, maxStreamElevation,numPourPointSteps ))';
            targetPourPointElevationList = round(targetPourPointElevationList);
    
            if(~isempty(targetPourPointElevationList))

                for stepNum = 1: length(targetPourPointElevationList)
    
                    targetPourPointElevation = targetPourPointElevationList(stepNum,1)
                    [closest, closestPourPointIndex] = ...
                        min(abs(streamLineZ - targetPourPointElevationList(stepNum)));
    
                    subcatchmentOutlet =...
                        [streamProfileStruct(streamCount).XCoords(closestPourPointIndex),...
                        streamProfileStruct(streamCount).YCoords(closestPourPointIndex)]; 
                    
                    finalOutletElevation = round(streamLineZ(closestPourPointIndex));

                    %Extract DEM of drainage basins based list of evenly spaced pour points

                    finalDrainageBasinArray =...
                        ExtractSubcatchmentDEM(flowDirectionObj, demArray,subcatchmentOutlet(:,1), subcatchmentOutlet(:,2));
    
                    % Write drainage basin DEM to vector
                    nY=size(finalDrainageBasinArray,2);
                    nX = size(finalDrainageBasinArray,1);
                    maxElevation = max(max(finalDrainageBasinArray));
                    minElevation = min(min(finalDrainageBasinArray));
                    demCoordIndex = find(~isnan(finalDrainageBasinArray));
                    demReshapedArray=reshape(finalDrainageBasinArray,nX*nY,1);
                    demReshapedArray = demReshapedArray(~isnan(demReshapedArray));
                    
                   
                    
                    % Calculate elevation PDF using KS density
                    % Uses default bandwidth 
                   
                        
                        [elevationProbability,binCenters]=ksdensity(demReshapedArray);
               
                    % Interpolate peak of the elevation PDF
                        elevationRange = minElevation:1:maxElevation;
                        splineInterpolatedPDF = spline(binCenters,elevationProbability,elevationRange);
                        [splineMax, splineMaxIndex] = max(splineInterpolatedPDF);
                    
                        hypsoPeak= elevationRange(splineMaxIndex); % Modal elevation, PDF peak
                        outlet=min(demReshapedArray); % Approximate outlet as minimum elevation of basin
                
                    % Round each modal elevation and outlet to single 
                        hypsoPeak = round(hypsoPeak);
                        outlet = round(outlet);

                    %allSubcatchmentStruct(streamCount) = struct('Outlets', outlet.stepNum, 'Mode',  hypsoPeak.stepNum)
                        allSubcatchmentDataArray(stepNum,1) = outlet;
                        allSubcatchmentDataArray(stepNum,2) = hypsoPeak;
        
                
                        % Find catchment that drains to defined channel
                        % head, write it out as DEM in temp folder
                    
                        
                        if(stepNum == length(targetPourPointElevationList))
          
                        highDumpName = ['Supercatchment',num2str(count),'StreamNum',num2str(streamNum),...
                            'PourPointElevation',num2str(finalOutletElevation),'.tif'];
                        highDumpFilePath=fullfile(highSubcatchmentDumpFile, highDumpName);
                        SubcatchmentWrite(finalDrainageBasinArray, geospatialReferenceArray,demInfo,nanFlag,highDumpFilePath ) 
                        end
                   
                   
                end
                
                % Once the loop ends, go to the most recently delineated catchment
                % Use its flow accumulation structure to walk up to the divide
        
                previousCount = stepNum;
        
                % Bring in catchment above channel head
        
                [subcatchmentDemArray,subcatchmentDemGrid, subcatchmentdemInfo, subcatchmentgeospatialReferenceArray] =...
                    loadDEM(highSubcatchmentDumpFile, highDumpName, nanFlag);
        
                % Delineate drainage for this small catchment
                

                subcatchmentDemSinksFilled = fillsinks(subcatchmentDemGrid);
                subcatchmentFlowDirectionObj = FLOWobj(subcatchmentDemSinksFilled);
                subcatchmentFlowAccumulationStruct  = flowacc(subcatchmentFlowDirectionObj);

                % Track upstream flow accumulation
        
                [aboveHeadZ, upstreamNodes] =...
                    TrackToDivide(subcatchmentDemGrid, subcatchmentFlowDirectionObj,subcatchmentFlowAccumulationStruct); 
    
                % Create vec of steps for catchment above channel head
                upstreamNodeArray = [aboveHeadZ, upstreamNodes];
                maxZAboveHead = max(aboveHeadZ);
                minZAboveHead = targetPourPointElevation;
                remainingSteps = (maxZAboveHead -  minZAboveHead)/pourPointElevationStep;
                remainingStepsNum = round(remainingSteps);

                remainingStepsList = (linspace(minZAboveHead, maxZAboveHead,remainingStepsNum ))';

                aboveHeadElevationList = round(remainingStepsList);
    
                % Avoid making subcatchment of only 1 pixel
                aboveHeadElevationListLength = length(aboveHeadElevationList);
                aboveHeadElevationList = aboveHeadElevationList(1:(aboveHeadElevationListLength-1));
                countNum = previousCount+1;
                
                % Repeat progressive hypsometry routine for above channel head catchments
            
                for stepNum = 2: length(aboveHeadElevationList)

                    % Get closest elevation above head to next step
                    [~, closestPourPointIndex] = min(abs(aboveHeadZ - aboveHeadElevationList(stepNum)));
  
                    aboveHeadPourPointElevation = aboveHeadElevationList(stepNum);
                    outletNode = upstreamNodeArray(closestPourPointIndex,2);
        
                     % Save elevation of closest step to finalOutletElevation
                     finalOutletElevation = round(aboveHeadZ(closestPourPointIndex));
        
                    [subcatchmentOutletX, subcatchmentOutletY] = ind2coord(subcatchmentDemGrid, outletNode);

                    finalDrainageBasinArray =...
                        ExtractSubcatchmentDEM(subcatchmentFlowDirectionObj,...
                        subcatchmentDemArray,subcatchmentOutletX, subcatchmentOutletY);

                    % Write basin DEM to vec
                
                    nY=size(finalDrainageBasinArray,2);
                    nX = size(finalDrainageBasinArray,1);
                    maxElevation = max(max(finalDrainageBasinArray));
                    minElevation = min(min(finalDrainageBasinArray));
                    demCoordIndex = find(~isnan(finalDrainageBasinArray));
                    demReshapedArray=reshape(finalDrainageBasinArray,nX*nY,1);
                    demReshapedArray = demReshapedArray(~isnan(demReshapedArray));
                    
                    
                
                    % Calculate elevation PDF
                    % Use default bandwidth
                    
                        [elevationProbability,binCenters]=ksdensity(demReshapedArray);
                
    
                        elevationRange = minElevation:1:maxElevation;
                        splineInterpolatedPDF = spline(binCenters,elevationProbability,elevationRange);
                        [splineMax, splineMaxIndex] = max(splineInterpolatedPDF);
                        hypsoPeak= elevationRange(splineMaxIndex);
                        outlet=min(demReshapedArray); % Approximate outlet elevation as min of drainage basin  
                
                        outlet = round(outlet);
                        hypsoPeak = round(hypsoPeak);

    
                        allSubcatchmentDataArray(countNum,1) = outlet;
                        allSubcatchmentDataArray(countNum,2) = hypsoPeak;

                        countNum=countNum+1;
                   
                end

                [outletList, orderedList] = sort(allSubcatchmentDataArray(:,1), 'ascend');
                hypsoPeakList = allSubcatchmentDataArray(orderedList,2);
                modeOutletArray = [outletList,hypsoPeakList];
    
                % Name output .txt file by Supercatchment number and chain number
    
                outputFileName = ['Supercatchment',num2str(count),'StreamNum',num2str(streamNum), '_ModeOutletArray.txt'];
                fullOutputFile = fullfile(outputFilePath, outputFileName);
                dlmwrite(fullOutputFile,modeOutletArray)

                
            end 
        end 
    
   
    
    
    end % Total progressive hypsometry loop




end % End supercatchment loop
