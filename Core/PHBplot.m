%% PHBplot
% 
% Description
% 
% For a list of PHBs, plot (h_change, h_bench) pairs
% 
%
Colors;
numSupers = length(supercatchmentNum);
supercatchmentNum = supercatchmentNum';

figureOutputFilePath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName, 'Figures');
mkdir(figureOutputFilePath);

%% Plot options
peakElevation = peakElevationForOutputFig;
plotColor;
markerSize = 25;
baseLevelOutlet =0;

%% Input and output files

%allFilepath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'PHBs','Cusum02_BenchLength3Steps','Tales')
%roiFileName = [groupName,'_allPHBs.txt']


supercatchmentFilepath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName,'PHBs','Cusum02_BenchLength3Steps','Tables')



allPHBs = dir(supercatchmentFilepath);
    j=0;
for i = 1:length(allPHBs)
    

    fileName = allPHBs(i).name;
    
    if(~isempty(strfind(fileName, 'Supercatchment')))
        
        j = j+1;
        allPHBList(j) = allPHBs(i);
    end
end


%% Gather list of supercatchment tables to plot from folder

targetList = 1:length(allPHBList)

for i = 1:length(targetList) 
    
    
    clear phbOutletArray modes outlets
    targetFolder = targetList(i);
    fileName = allPHBList(targetFolder).name;
    
   
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, '_all');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    superNumFileList(i) = str2double(inBetweenTextSuper); 
    
end

%% Plot list of specified supers

[~,targetSupersList] = find(superNumFileList==supercatchmentNum)
    
for j = 1:length(targetSupersList)
    
    superNum = superNumFileList(targetSupersList(j));
    supercatchmentFigureOutputFilePath = fullfile(phAnalysisFilePath,groupArea,AcSubFolderName, 'Figures','SupercatchmentPHBs');
    mkdir(supercatchmentFigureOutputFilePath);
    
    benchSuperFolderName = ['Supercatchment', num2str(superNum),'PHBs']
    benchSuperListName = ['Supercatchment', num2str(superNum),'_allOutletModePairs.txt'];
    benchSuperOutName = ['Supercatchment', num2str(superNum)];
    phbOutletArray = dlmread(fullfile(supercatchmentFilepath, benchSuperListName),'\t',2);

    modeOutletArray = phbOutletArray(:,[1,2]);
    modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
    zeroIndices = find(modeOutletArray(:,2)==0);
    modeOutletArray(zeroIndices,:)=[];
   
    [~, ia, ic]  = unique(phbOutletArray, 'rows');
    modeOutletArrayUnique = phbOutletArray(ia,[1,2]);
    zeroIndices = find(modeOutletArrayUnique(:,2)==0);
    modeOutletArrayUnique(zeroIndices,:)=[];

    %modes = modeOutletArrayUnique(:,2)
    %outlets = modeOutletArrayUnique(:,1)
    %prunedDeltaH = modes-outlets;

    modes = modeOutletArray(:,2)
    outlets = modeOutletArray(:,1)
    deltaH = modes-outlets;
    [deltaHProb, deltaHBinCenters,bw] = ksdensity(deltaH,'Support','positive','BoundaryCorrection','reflection');
    [deltaHProb, deltaHBinCenters] = ksdensity(deltaH,'Support','positive','BoundaryCorrection','reflection','Bandwidth',bw*.75);

    %modes = modeOutletArrayUnique(:,2)
    %outlets = modeOutletArrayUnique(:,1)
   

    

%% Make PHB plot

    figure(1)
    clf;
    hold on;
    %grid on
    set(gca,'xtick',[0:1000:peakElevation], 'Linewidth', 1.5)
    set(gca,'ytick',[0:500:peakElevation], 'Linewidth', 1.5)
    set(gca,'fontname', 'Arial','FontSize',25, 'fontweight', 'bold')
    
    

    plot(outlets,modes,'.','color',plotColor,'MarkerSize', markerSize);
    
    ylim([1800 peakElevation]);
    xlim([1800 peakElevation]);


    ylabel('Modal elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'none')
    xlabel('Outlet elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')
    
    reflineSeeds = linspace(0, peakElevation);
    reflinePoints = linspace(0, peakElevation);
    
%     upperSeeds = linspace(upperMode, peakElevation);
%     upperPoints = linspace(upperMode, peakElevation);
%     
%     lowerSeeds = linspace(0, peakElevation);
%     lowerPoints = linspace(0, peakElevation);
    
    plot(reflineSeeds, reflinePoints, 'Linewidth',1.5,'Color', black)
    %plot(reflineSeeds, reflinePoints+upperMode, 'Linewidth',1.5,'Color', redOrange)
    %plot(reflineSeeds, reflinePoints+lowerMode, 'Linewidth',1.5,'Color', blue)
   
    
    
    fig1Name = [benchSuperOutName,'_hBench_vs_hChange.png']
    fig1FilePath = fullfile(supercatchmentFigureOutputFilePath, fig1Name);
    fig1 = gcf;
    saveas(fig1,fig1FilePath,outputFigType)
    
%%

%% Make deltaH plot

    figure(2)
    clf;
    hold on;
    %grid on
    set(gca,'xtick',[0:1000:peakElevation], 'Linewidth', 1.5)
    set(gca,'ytick',[0:500:pdfHeight], 'Linewidth', 1.5)
    set(gca,'fontname', 'Arial','FontSize',25, 'fontweight', 'bold')
    
    

    plot(outlets,(modes-outlets),'.','color',plotColor,'MarkerSize', markerSize);
    
    ylim([0 pdfHeight]);
    xlim([0 peakElevation]);


    ylabel('$${\Delta{h}}$$ [-]', 'Fontsize', 25,'Interpreter', 'latex','Fontweight', 'bold')
    xlabel('Outlet elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')
    
    reflineSeeds = linspace(0, pdfHeight);
    reflinePoints = linspace(0, pdfHeight);
    
%     upperSeeds = linspace(upperMode, peakElevation);
%     upperPoints = linspace(upperMode, peakElevation);
%     
%     lowerSeeds = linspace(0, peakElevation);
%     lowerPoints = linspace(0, peakElevation);
    
    %plot(reflineSeeds, reflinePoints, 'Linewidth',1.5,'Color', black)
    %plot(reflineSeeds, reflinePoints+upperMode, 'Linewidth',1.5,'Color', redOrange)
    %plot(reflineSeeds, reflinePoints+lowerMode, 'Linewidth',1.5,'Color', blue)
   
    
    
    fig2Name = [benchSuperOutName,'_hBench_vs_hChange.png']
    fig2FilePath = fullfile(supercatchmentFigureOutputFilePath, fig2Name);
    fig2 = gcf;
    %saveas(fig1,fig1FilePath,outputFigType)
    
    
    %% Plot deltaH PDF
    
    figure(3)
    clf;
    hold on;
    %grid on;
    set(gca,'FontSize',25, 'fontweight', 'bold')

    plot(deltaHProb, deltaHBinCenters, 'color', plotColor, 'Linewidth', 3)

    ylabel('$${\Delta{h}}$$ [m]', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'latex')
    xlabel('p($${\Delta{h}}$$) $$[m^{-1}]$$ ', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'latex')
    
    ylim([0 pdfHeight]);
    %xlim([0 peakElevation]);



    %Save figure



    fig3Name = [benchSuperOutName,'_deltaH_Distribution.png']
    fig3FilePath = fullfile(supercatchmentFigureOutputFilePath, fig3Name);
    fig3 = gcf;
    saveas(fig3,fig3FilePath,outputFigType)



 

    
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

