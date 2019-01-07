%% Supercatchment extraction
%
% Description
% 
% All catchments above a specified outlet elevation and above a specified
% peak elevation are extracted and written. 
%
% Output
%
% All catchment DEMs are written to specified folder as .tif

%% Set filepath and variables

DefinePHroot;
addpath(topoToolboxFilePath); 

targetOutletElevation = 540;
AreaThresh = 2777; % In pixels, minimum number of accumulating pixels for to define a channel head
supercatchmentPeakElevation = 1500; % Usually a value high enough to incoporate all catchments that reach the main divide
pourPointSensitivtiyThreshold = 30; % In meters, how close can pourpt be to target?

groupNameFilePath =  fullfile(phDataFilePath,groupArea,'DEM');
groupNameRegionDEM =  ['WesternAlps_MatterhornROI_1arcSecond_SRTM_patchedWithResampled90m_30mRes_Float32_nullSet.tif'];

demNanFlag = -32768; % Nan for Incoming DEM
subcatchmentNanFlag = -32768; % NaN for outgoing Subcatchments

outputFilePath=...
    fullfile(phDataFilePath, groupArea,'Supercatchments_540m');

mkdir(outputFilePath);


%% Bring in projected DEM of target region, write to GridObj

[demArray,demGrid, demInfo, geospatialReferenceArray] =...
    loadDEM(groupNameFilePath, groupNameRegionDEM, demNanFlag);

%% Fill sinks, delineate channel network

% DEM preprocessing, including pit (sink) filling, flow direction
% calculation. 

demSinksFilled = fillsinks(demGrid);
flowDirectionObj = FLOWobj(demSinksFilled);
flowAccumulationStruct  = flowacc(flowDirectionObj);

% Use threshold area to delineate channel network, write to StreamObj

streamAreaThreshold=flowAccumulationStruct>AreaThresh; 
streamNetworkStreamObj =...
     STREAMobj(flowDirectionObj,streamAreaThreshold);
 
streamNetworkStruct=STREAMobj2mapstruct(streamNetworkStreamObj);

%% Create a structure array of chains (channels) (x,y,z)

for streamCount= 1:length(streamNetworkStruct)
    
    % Use dummy first-order channel links to define all "channels" in DEM
    
    if(streamNetworkStruct(streamCount).streamorder==1)
        
        streamNum = streamNetworkStruct(streamCount).IX;
        
        % Get x,y,z for each channel pixel
        [streamX, streamY, streamZ,~,~] =...
            longProfileGenerator(streamNetworkStruct, streamNetworkStreamObj,...
            demGrid, streamNum, streamAreaThreshold);
     
        % Generate a structure array that contains x,y,z of each channel
        streamProfileStruct(streamCount) = ...
            struct('XCoords',streamX,'YCoords',streamY, 'Elevation', streamZ);
    end
end

%% Find outlets closest to target low elevation

for streamCountNum=1:length(streamProfileStruct)
    
    if(~isempty(streamProfileStruct(streamCountNum).Elevation))
    
    %For each channel, find elevation closest to target
    profileZ = streamProfileStruct(streamCountNum).Elevation;
    tempDiffVec = abs(profileZ-targetOutletElevation);
    [closestToTargetElevation, minimumIndex] = min(tempDiffVec); 
    
    % Check if channel elevation and tartget outlet are within imposed
    % threshold
    
        if(closestToTargetElevation<pourPointSensitivtiyThreshold)
        
        % Create structure array of low elevation channel pixels
        % These pixels will be used to extract catchments above target
        % outlet elevation     
        
        lowElevationBasinOutlets(streamCountNum,:)=...
        [streamProfileStruct(streamCountNum).XCoords(minimumIndex),...
        streamProfileStruct(streamCountNum).YCoords(minimumIndex)];
    
        end
    
   end
end

% Prune repeated low elevation outlets
lowElevationBasinOutlets = unique(lowElevationBasinOutlets,'rows');



%% Extract and write supercatchments above target elevation

supercatchmentCountNum = 0; % All catchment counter
supercatchmentNameNum = 0; % Extracted (above peak elevatin thresh only) counter

for outletCount = 1:length(lowElevationBasinOutlets)
     
    finalDrainageBasinArray = ExtractSubcatchmentDEM(flowDirectionObj, demArray,...
    lowElevationBasinOutlets(outletCount,1), lowElevationBasinOutlets(outletCount,2));

    maxDrainageBasinElevation = max(max(finalDrainageBasinArray));


    % Check if extracted catchment peak is above imposed peak elevation
    if(maxDrainageBasinElevation>supercatchmentPeakElevation)
        
        supercatchmentCountNum=supercatchmentCountNum+1;
        supercatchmentNameNum = supercatchmentNameNum+1;
    
        outputFileName=...
            [groupArea,'Supercatchment',num2str(supercatchmentNameNum),'.tif'];

        fullOutputFile = fullfile(outputFilePath, outputFileName);

        SubcatchmentWrite(finalDrainageBasinArray, geospatialReferenceArray,demInfo,...
            subcatchmentNanFlag, fullOutputFile)  
    
   
    else
        
     supercatchmentCountNum = supercatchmentCountNum+1;
     
    end
    
end


