%% Generate list of PHBs and write to txt file
% 
% Description
% 
% Find steps in PH output, and writes a .txt containing (outlet, mode)
% pairs for each PHB in supercatchment. For a list of supercatchments, two 
% .txt's are written: one master file that contains all PHBs for all
% supercatchments passed, and one for each supercatchment.
% 
%% Set variables

minBenchLength; % In catchments, consecutive number of catchments with shared PH mode 
spillOverElevations = 25; % Elevation above and below PH modal elevation
areaThreshPixelNum = Ac; % In pixels, must be same as in PH
cuSumThresh = .02; % For PH bench separation
pixelLength = 30;
nanFlag = -32768;

allROIPHBcount=1;
%% 
for count = supercatchmentNum
    clear supercatchmentPHBArray
    PHBCount=1;
    streamSupercatchment = count;


    supercatchmentFilePath = fullfile(phDataFilePath,groupArea,'Supercatchments');
    supercatchmentFileName = ['Supercatchment', num2str(streamSupercatchment)];
    supercatchmentOutFileName = ['Supercatchment', num2str(streamSupercatchment),'PHBs'];
    supercatchmentDemName = [groupArea, 'Supercatchment', num2str(streamSupercatchment),'.tif']

    %Input list of (outlet, mode) pairs
    progressivePourPointSubcatchmentFilePath = fullfile(fullfile(phAnalysisFilePath,groupArea,'Subcatchments','25mStep', num2str(supercatchmentFileName)));
    
    %Output file path for PHB layer
    allSupercatchmentPHBfilePath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt');
    outputFileNameROI = [groupArea,'_allPHBs.txt'];
    
    SupercatchmentBenchFiles = fullfile(phAnalysisFilePath, groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt','Supercatchments',supercatchmentOutFileName);
    
    outputFileName = [supercatchmentFileName,'_allOutletModePairs.txt'];

    mkdir(allSupercatchmentPHBfilePath);
    mkdir(SupercatchmentBenchFiles);


    
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

                    superCount = PHBCount;
         
                    supercatchmentPHBArray(superCount,1) =benchOutletElevation;
                    supercatchmentPHBArray(superCount,2) =hypsoPeakElevation;
                    supercatchmentPHBArray(superCount,3) = count;
                    supercatchmentPHBArray(superCount,4) = streamNum;
                    
                    PHBCount = PHBCount+1;
                    
                    fullOutputFileForSuper = fullfile(SupercatchmentBenchFiles, outputFileName);

                    columnNames = {'Outlets','Modes','Supercatchment', 'StreamID'};
                    supercatchmentOutTable = array2table(supercatchmentPHBArray,'VariableNames',columnNames)
                    writetable(supercatchmentOutTable, fullOutputFileForSuper, 'Delimiter','\t');
               
                   
               
               
        end
        end
        end
        end
        
        
        superPHBListLength = length(supercatchmentPHBArray);
        allROIPHBcountNew = (allROIPHBcount+ superPHBListLength)-1;
        allROIArray(allROIPHBcount:allROIPHBcountNew,:) = supercatchmentPHBArray;
        
        fullOutputFileForSuper = fullfile(allSupercatchmentPHBfilePath, outputFileNameROI);

                    columnNames = {'Modes','Outlets','Supercatchment', 'StreamID'};
                    roiOutTable = array2table(allROIArray,'VariableNames',columnNames)
                    writetable(roiOutTable, fullOutputFileForSuper, 'Delimiter','\t');
        allROIPHBcount = allROIPHBcount+superPHBListLength;
               
end

 





