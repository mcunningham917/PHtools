%% PHBplot
% 
% Description
% 
% For a list of PHBs, plot (h_change, h_bench) pairs
% 
%
Colors;
count = supercatchmentNum;
numSupers = length(supercatchmentNum);

figureOutputFilePath = fullfile(phAnalysisFilePath,groupArea,'halfSqKmAc', 'Figures');
mkdir(figureOutputFilePath);

%% Plot options
peakElevation = peakElevationForOutputFig;
plotColor;
markerSize = 25;
baseLevelOutlet =0;

%% Input and output files

allFilepath = fullfile(phAnalysisFilePath,groupArea,'halfSqKmAc','PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt')
supercatchmentFilepath = fullfile(phAnalysisFilePath,groupArea,'halfSqKmAc','PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt','Supercatchments')
roiFileName = [groupName,'_allPHBs.txt']


allPHBs = dir(supercatchmentFilepath);
    j=0;
for i = 1:length(allPHBs)
    

    fileName = allPHBs(i).name;
    
    if(~isempty(strfind(fileName, 'Supercatchment')))
        
        j = j+1;
        allPHBList(j) = allPHBs(i);
    end
end


%%
targetList = 1:length(allPHBList)

for i = 1:length(targetList) 
    
    
    clear phbOutletArray modes outlets
    targetFolder = targetList(i);
    fileName = allPHBList(targetFolder).name;
    
   
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, 'PHBs');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    superNum = str2double(inBetweenTextSuper); 
    
    %if(length(count)>i)
    if(superNum == count)
    
    
    supercatchmentFigureOutputFilePath = fullfile(phAnalysisFilePath,groupArea,'halfSqKmAc', 'Figures','SupercatchmentPHBs');
    mkdir(supercatchmentFigureOutputFilePath);
    
    benchSuperFolderName = ['Supercatchment', num2str(superNum),'PHBs']
    benchSuperListName = ['Supercatchment', num2str(superNum),'_allOutletModePairs.txt'];
    benchSuperOutName = ['Supercatchment', num2str(superNum)];
    phbOutletArray = dlmread(fullfile(supercatchmentFilepath,benchSuperFolderName, benchSuperListName),'\t',2);

%%

     modeOutletArray = phbOutletArray(:,[1,2]);
    modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
    zeroIndices = find(modeOutletArray(:,2)==0);
    modeOutletArray(zeroIndices,:)=[];
    
    
    [~, ia, ic]  = unique(phbOutletArray, 'rows');
    modeOutletArrayUnique = phbOutletArray(ia,[1,2]);
    zeroIndices = find(modeOutletArrayUnique(:,2)==0);
    modeOutletArrayUnique(zeroIndices,:)=[];

%     modes = modeOutletArrayUnique(:,2)
%     outlets = modeOutletArrayUnique(:,1)
%     prunedDeltaH = modes-outlets;

    modes = modeOutletArray(:,2)
    outlets = modeOutletArray(:,1)
    deltaH = modes-outlets;
    [deltaHProb, deltaHBinCenters] = ksdensity(deltaH,'Support','positive','BoundaryCorrection','reflection');

    %modes = modeOutletArrayUnique(:,2)
    %outlets = modeOutletArrayUnique(:,1)
   

    

%% Make PHB plot

    figure(1)
    clf;
    hold on;
    grid on
    set(gca,'xtick',[0:1000:peakElevation], 'Linewidth', 1.5)
    set(gca,'ytick',[0:500:peakElevation], 'Linewidth', 1.5)
    set(gca,'fontname', 'Arial','FontSize',25, 'fontweight', 'bold')
    
    

    plot(outlets,modes,'.','color',plotColor,'MarkerSize', markerSize);

    ylim([0 peakElevation]);
    xlim([0 peakElevation]);


    ylabel('Modal elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'none')
    xlabel('Outlet elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')
    
    fig1Name = [benchSuperOutName,'_hBench_vs_hChange.png']
    fig1FilePath = fullfile(supercatchmentFigureOutputFilePath, fig1Name);
    fig1 = gcf;
    saveas(fig1,fig1FilePath,outputFigType)
    
    
    
    %% Plot deltaH PDF
    
    figure(2)
    clf;
    hold on;
    grid on;
    set(gca,'FontSize',18, 'fontweight', 'bold')

    plot(deltaHProb, deltaHBinCenters, 'color', plotColor, 'Linewidth', 2)

    ylabel('$${\Delta{h}}$$ [-]', 'Fontsize', 25, 'Fontweight', 'normal','Interpreter', 'latex')
    xlabel('p($${\Delta{h}}$$) ', 'Fontsize', 25, 'Fontweight', 'normal', 'Interpreter', 'latex')
    
    ylim([0 pdfHeight]);
    %xlim([0 peakElevation]);



    % Save figure



    fig2Name = [benchSuperOutName,'_deltaH_Distribution.png']
    fig2FilePath = fullfile(supercatchmentFigureOutputFilePath, fig2Name);
    fig2 = gcf;
    saveas(fig2,fig2FilePath,outputFigType)



 

    
    end
end




%% Plot all PHBs for ROI
% 
% phbOutletArray = dlmread(fullfile(allFilepath, roiFileName),'\t',2);
% modeOutletArray = phbOutletArray(:,[1,2]);
% modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
% zeroIndices = find(modeOutletArray(:,2)==0);
% modeOutletArray(zeroIndices,:)=[];
% 
% modes = modeOutletArray(:,2)
% outlets = modeOutletArray(:,1)
% 
% 
% 
% modes = modeOutletArray(:,2)
% outlets = modeOutletArray(:,1)
% 
% % Make all ROI PHB plot
% 
% figure(1)
% hold on
% grid on
% set(gca,'xtick',[0:1000:peakElevation], 'Linewidth', 1.5)
% set(gca,'ytick',[0:500:peakElevation], 'Linewidth', 1.5)
% set(gca,'fontname', 'Arial','FontSize',25, 'fontweight', 'bold')
% 
% plot(outlets,modes,'.','color',plotColor,'MarkerSize', markerSize);
% 
% ylim([0 peakElevation]);
% xlim([0 peakElevation]);
% 
% 
% ylabel('Modal elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'none')
% xlabel('Outlet elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')
% 
% yIntercept=500;
% upperBoundSlope =1;
% boundSeeds = linspace(0, peakElevation, length(phbOutletArray(:,1)));
% upperBoundLine = yIntercept+(boundSeeds*upperBoundSlope)
% plot(boundSeeds,upperBoundLine, 'r-', 'linewidth',1)
% 
% %%
% 
% figName = [groupArea,'_hBench_vs_hChange_allSupercatchments.png']
% figFilePath = fullfile(figureOutputFilePath, figName);
% fig1 = gcf;
% saveas(fig1,figFilePath,outputFigType)

