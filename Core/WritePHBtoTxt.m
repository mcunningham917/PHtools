%% Generate list of PHBs and write to txt file
% 
% Description
% 
% Find steps in PH output, and writes a geotiff of elevation band around
% the PH mode. Also writes out a geotiff of entire catchment above h_change
% for the purpose of generating polygon layer.
% 
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

    %Input list of (outlet, mode) pairs
    progressivePourPointSubcatchmentFilePath = fullfile(fullfile(phAnalysisFilePath,groupArea,'Subcatchments','25mStep', num2str(supercatchmentFileName)));
    
    %Output file path for PHB layer
    allSupercatchmentPHBfilePath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchments');

    mkdir(allSupercatchmentPHBfilePath);



    
    %% Delineate channel network for supercatchment
    
    [supercatchmentDemArray, supercatchmentDemGrid, supercatchmentDemInfo, supercatchmentGeospatialReferenceArray] =...
        loadDEM(supercatchmentFilePath, supercatchmentDemName, nanFlag);

    demSinksFilled = fillsinks(supercatchmentDemGrid);
    flowDirectionObj = FLOWobj(demSinksFilled);
    flowAccumulationStruct  = flowacc(flowDirectionObj);
    streamNodeThreshNum=areaThreshPixelNum; 
    streamNodeThreshold = flowAccumulationStruct>streamNodeThreshNum;
    streamNetworkStreamObj = STREAMobj(flowDirectionObj,streamNodeThreshold);
    streamNetworkStruct=STREAMobj2mapstruct(streamNetworkStreamObj);
    

    %Get list of chains
clear realFirstOrderStreamList
    firstOrderStreamList = dir(progressivePourPointSubcatchmentFilePath);

        for i = 1:length(firstOrderStreamList)
    
            realFolder = strfind(firstOrderStreamList(i).name,'Super')
    
            if(~isempty(realFolder))
                realFirstOrderStreamList(i) = i;
            end
            
        end
    
realFirstOrderStreamList=realFirstOrderStreamList(realFirstOrderStreamList>0);

%%
    
%Find benches for each first order tributary in the supercatchment

        for highTributaryNum = 1:length(realFirstOrderStreamList)
            clear allSubcatchmentDataArray hypsoPeakList outletList
    
            chainNum = realFirstOrderStreamList(highTributaryNum);
            chainNumName = fullfile(progressivePourPointSubcatchmentFilePath, firstOrderStreamList(chainNum).name);
            allSubcatchmentDataArray= dlmread(chainNumName);
            
            streamNumFile = firstOrderStreamList(chainNum).name;
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
                
                if(hypsoPeakElevation<1000)
                            hypsoPeakName = num2str( hypsoPeakElevation, '%04d' );
                        else
                            hypsoPeakName = num2str(hypsoPeakElevation)
                end
          
                ModeOutletPair = [benchOutletElevation, hypsoPeakElevation]; 
             
                benchOutputFileName = ['HypsoPeak', hypsoPeakName,'PourPointElevation',num2str(benchOutletElevation),...
                    'Supercatchment',num2str(streamSupercatchment),'StreamNum',num2str(streamNum)];

               fullOutputFileForSuper = fullfile(allSupercatchmentPHBfilePath, benchOutputFileName);
               
                dlmwrite(fullOutputFileForSuper,ModeOutletPair);
               
               
        end
        end
        end
        end
end

 





