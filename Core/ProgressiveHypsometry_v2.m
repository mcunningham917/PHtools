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

pourPointElevationStep = phStepLength; % Step elevation height--vertical spacing of catchment elevations
pourPointElevationSensitivity = 5; % Variability around each step height
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
    
    clear streamProfileStruct branchHypsometryStruct
    
    
    supercatchmentFileName = [groupArea,'Supercatchment',num2str(count),'.tif'];
    supercatchmentFolderName = ['Supercatchment', num2str(count)];
    subcatchmentFolderName = [num2str(pourPointElevationStep),'mStep'];
    outputFilePath=fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'Subcatchments',subcatchmentFolderName, supercatchmentFolderName);
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


% for streamCount= 1:length(streamNetworkStruct)
%     
%     if(streamNetworkStruct(streamCount).streamorder==1)
%         
%         streamIndex = streamNetworkStruct(streamCount).IX;
%         
%         %Extract their profiles using longProfileGenerator
%         %Use only raw Z values (not smoothed)
%         [streamX, streamY, streamZ, ~, ~, ~] =...
%             longProfileGenerator(streamNetworkStruct, streamNetworkStreamObj, demGrid, streamIndex, streamNodeThreshNum);
%      
%         %Build struct with X,Y,Z for each streamline
%         streamProfileStruct(streamCount) = ...
%             struct('XCoords',streamX,'YCoords',streamY, 'Elevation', streamZ, 'StreamNum', streamIndex);
%     end
% end


%% Calculate hypsometry along each Branch

for allStreamCount = 1:length(streamNetworkStruct)
    
        clear branchZ branchX branchY branchIX targetPourPointElevationList
    
        branchX=streamNetworkStruct(allStreamCount).X';
        branchY=streamNetworkStruct(allStreamCount).Y';
        branchX = branchX(~isnan(branchX));
        branchY = branchY(~isnan(branchY));
        branchIX = coord2ind(demGrid,branchX,branchY);
        branchIX = branchIX(1:(length(branchIX)));
        branchZ = demArray(branchIX);
        
        maxStreamElevation = max(branchZ); % Elevation of defined channel head
        minStreamElevation = min(branchZ);
        
        % Make a list of target pourpoints based on total height drop of 
        numPourPointSteps = (maxStreamElevation - minStreamElevation)/pourPointElevationStep;
        targetPourPointElevationList = (linspace(minStreamElevation, maxStreamElevation,numPourPointSteps ))';
        targetPourPointElevationList = round(targetPourPointElevationList);
        
        
        for stepNum = 1: length(targetPourPointElevationList)
    
                    targetPourPointElevation = targetPourPointElevationList(stepNum,1)
                    [closest, closestPourPointIndex] = ...
                        min(abs(branchZ - targetPourPointElevationList(stepNum)));
    
                    subcatchmentOutlet =...
                        [branchX(closestPourPointIndex),branchY(closestPourPointIndex)]; 
                    
                    finalOutletElevation = round(branchZ(closestPourPointIndex));

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
                        branchOutlets = round(hypsoPeak);
                        branchModes = round(outlet);
                        
                        allSubcatchmentDataArray(stepNum,1) = outlet;
                        allSubcatchmentDataArray(stepNum,2) = hypsoPeak;
   
    
                        
                       

                        

        end
       
        if(exist('allSubcatchmentDataArray','var'))
         branchHyspometryStruct(allStreamCount) = struct('BranchId', streamNetworkStruct(allStreamCount).IX, 'Outlets', allSubcatchmentDataArray(:,1), 'Modes', allSubcatchmentDataArray(:,2))
         clear allSubcatchmentDataArray;
        else
             branchHyspometryStruct(allStreamCount) = struct('BranchId', streamNetworkStruct(allStreamCount).IX, 'Outlets', NaN, 'Modes', NaN)
        end
end



%%

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

%%

for streamCount = 1:length(streamProfileStruct)
    
    clear steamLineZ
        
    streamNum = streamProfileStruct(streamCount).StreamNum;
    
    if(~isempty(streamNum) && ~isempty(branchHyspometryStruct(streamNum).BranchId))
        
    streamLineZ = streamProfileStruct(streamCount).Elevation;
    maxStreamElevation = max(streamLineZ); % Elevation of defined channel head
    minStreamElevation = min(streamLineZ);

   % Make a list of target pourpoints based on total height drop of 
    numPourPointSteps = (maxStreamElevation - minStreamElevation)/pourPointElevationStep;
    
    if(numPourPointSteps<1)
        
        highSubcatchmentPourPoint=max(streamLineZ)
    else
        
        targetPourPointElevationList = (linspace(minStreamElevation, maxStreamElevation,numPourPointSteps ))';
        targetPourPointElevationList = round(targetPourPointElevationList);
        
        highSubcatchmentPourPoint = max(targetPourPointElevationList);
    end
            
                    [closest, closestPourPointIndex] = ...
                        min(abs(streamLineZ - highSubcatchmentPourPoint));
    
                    subcatchmentOutlet =...
                        [streamProfileStruct(streamCount).XCoords(closestPourPointIndex),...
                        streamProfileStruct(streamCount).YCoords(closestPourPointIndex)]; 
                    
                    finalOutletElevation = round(streamLineZ(closestPourPointIndex));

                    %Extract DEM of drainage basins based list of evenly spaced pour points

                    finalDrainageBasinArray =...
                        ExtractSubcatchmentDEM(flowDirectionObj, demArray,subcatchmentOutlet(:,1), subcatchmentOutlet(:,2));

    highDumpName = ['Supercatchment',num2str(streamCount),'StreamNum',num2str(streamNum),...
                            'PourPointElevation',num2str(finalOutletElevation),'.tif'];
                        highDumpFilePath=fullfile(highSubcatchmentDumpFile, highDumpName);
                        SubcatchmentWrite(finalDrainageBasinArray, geospatialReferenceArray,demInfo,nanFlag,highDumpFilePath )


                [subcatchmentDemArray,subcatchmentDemGrid, subcatchmentdemInfo, subcatchmentgeospatialReferenceArray] =...
                    loadDEM(highSubcatchmentDumpFile, highDumpName, nanFlag);
                
               subcatchmentDEMVec = subcatchmentDemArray(~isnan(subcatchmentDemArray));
               superDemVec = demArray(~isnan(demArray));
               
               if(((length(subcatchmentDEMVec)/length(superDemVec))<0.1))
                
              
        
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
                minZAboveHead = finalOutletElevation;
                remainingSteps = (maxZAboveHead -  minZAboveHead)/pourPointElevationStep;
                remainingStepsNum = round(remainingSteps);

                remainingStepsList = (linspace(minZAboveHead, maxZAboveHead,remainingStepsNum ))';

                aboveHeadElevationList = round(remainingStepsList);
    
                % Avoid making subcatchment of only 1 pixel
                aboveHeadElevationListLength = length(aboveHeadElevationList);
                aboveHeadElevationList = aboveHeadElevationList(1:(aboveHeadElevationListLength-1));
             
                
                % Repeat progressive hypsometry routine for above channel head catchments
            
                for stepNum = 1: length(aboveHeadElevationList)

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

    
                        aboveHeadArray(stepNum,1) = outlet;
                        aboveHeadArray(stepNum,2) = hypsoPeak;

                end
                
                %%
                targetTrib = streamNetworkStruct(streamNum);

            % Find each downstream link for specified channel link
            for tribCount = 1:100
                if(~isnan(targetTrib.tribtoIX) && targetTrib.tribtoIX>0)
                targetStreamTribs(tribCount,1) = targetTrib.tribtoIX;
                targetTribNum = targetTrib.tribtoIX;
                targetTrib = streamNetworkStruct(targetTribNum);
                tribCount = tribCount +1;

                else
                break
                end
            end
%%
% Check that stream profile has been built, and write structure array with
% x,y for each channel pixel, organized by channel link

if(exist('targetStreamTribs','var'))
allTribFields = [streamNetworkStruct(streamNum).IX; targetStreamTribs(:,1)];
else 
    targetTribNum = targetTrib.IX;
         allTribFields = targetTribNum;
         
end
%%
for j = 1:length(allTribFields)
    
    if j==1
    branchPHStruct(j) = struct('Outlets',branchHyspometryStruct(allTribFields(j)).Outlets,...
        'Modes', branchHyspometryStruct(allTribFields(j)).Modes);
    
    %upstreamMinElevation = min(branchHyspometryStruct(allTribFields(j)).Outlets);
    else
             upstreamElevationList = vertcat(branchHyspometryStruct(allTribFields(1:j-1)).Outlets); 
             upstreamMinElevation = min(upstreamElevationList);
         
     
           
         branchIX= branchHyspometryStruct(allTribFields(j)).Outlets<upstreamMinElevation;
         branchPHStruct(j) = struct('Outlets',branchHyspometryStruct(allTribFields(j)).Outlets(branchIX),...
        'Modes',branchHyspometryStruct(allTribFields(j)).Modes(branchIX))
    end
 end
    

 %%   
    
     
        
    
if exist ('aboveHeadArray','var')
    
branchPHOutlets = vertcat(aboveHeadArray(:,1),branchPHStruct.Outlets) 
branchPHModes =  vertcat(aboveHeadArray(:,2),branchPHStruct.Modes) 

else
    
branchPHOutlets = vertcat(branchPHStruct.Outlets) 
branchPHModes =  vertcat(branchPHStruct.Modes) 
end
               
[allBranchOutlets, sortIX] = sort(branchPHOutlets, 'ascend')
allBranchModes = branchPHModes(sortIX);


modeOutletArray = [allBranchOutlets,allBranchModes];
[~, ia, ic]  = unique(modeOutletArray, 'rows');
modeOutletArrayUnique = modeOutletArray(ia,[1,2]);
nanIX = ~isnan(modeOutletArrayUnique(:,1));
modeOutletArrayUnique = modeOutletArrayUnique(nanIX,[1,2]);
clear nanIX

    
                % Name output .txt file by Supercatchment number and chain number
    
                outputFileName = ['Supercatchment',num2str(count),'StreamNum',num2str(streamNum), '_ModeOutletArray.txt'];
                fullOutputFile = fullfile(outputFilePath, outputFileName);
                dlmwrite(fullOutputFile,modeOutletArrayUnique)
                clear aboveHeadArray branchPHOutlets branchPHModes branchPHStruct
               
    end
end
end
end
               
    


 % End supercatchment loop


               