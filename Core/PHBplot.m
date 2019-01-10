%% PHBplot
% 
% Description
% 
% For a list of PHBs, plot (h_change, h_bench) pairs
% 
%

count = supercatchmentNum;

figureOutputFilePath = fullfile(phAnalysisFilePath,groupArea, 'Figures');
mkdir(figureOutputFilePath);

%% Plot options
peakElevation = peakElevationForOutputFig;
plotColor = 'b'
markerSize = 25;
baseLevelOutlet =0;

%% Input and output files

filepath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt')
supercatchmentFigureFilepath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt','Supercatchments')
roiFileName = [groupArea, '_allPHBs.txt']


allPHBs = dir(filepath);
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
    
    clf;
    clear phbOutletArray
    targetFolder = targetList(i);
    fileName = allPHBList(targetFolder).name;
    
   
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, 'PHBs');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    superNum = str2double(inBetweenTextSuper); 
    
    
    supercatchmentFigureOutputFilePath = fullfile(phAnalysisFilePath,groupArea, 'Figures','SupercatchmentPHBs');
    mkdir(supercatchmentFigureOutputFilePath);
    
    benchSuperFolderName = ['Supercatchment', num2str(superNum),'PHBs']
    benchSuperListName = ['Supercatchment', num2str(superNum),'_allOutletModePairs.txt'];
    benchSuperOutName = ['Supercatchment', num2str(superNum)];
    phbOutletArray = dlmread(fullfile(filepath,benchSuperFolderName, benchSuperListName),'\t',2);

%%

    modeOutletArray = phbOutletArray(:,[1,2]);
    modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
    zeroIndices = find(modeOutletArray(:,2)==0);
    modeOutletArray(zeroIndices,:)=[];

    modes = modeOutletArray(:,2)
    outlets = modeOutletArray(:,1)

    modes = modeOutletArray(:,2)
    outlets = modeOutletArray(:,1)

%% Make PHB plot

    figure(superNum)
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



%% Save figure

figName = [benchSuperOutName,'_hBench_vs_hChange.png']
figFilePath = fullfile(supercatchmentFigureOutputFilePath, figName);
fig1 = gcf;
saveas(fig1,figFilePath,outputFigType)
end




%% Plot all PHBs for ROI

phbOutletArray = dlmread(fullfile(filepath, roiFileName),'\t',2);
modeOutletArray = phbOutletArray(:,[1,2]);
modeOutletArray(any(isnan(modeOutletArray), 2), :) = [];
zeroIndices = find(modeOutletArray(:,2)==0);
modeOutletArray(zeroIndices,:)=[];

modes = modeOutletArray(:,2)
outlets = modeOutletArray(:,1)



modes = modeOutletArray(:,2)
outlets = modeOutletArray(:,1)

%% Make all ROI PHB plot

figure(1)
hold on
grid on
set(gca,'xtick',[0:1000:peakElevation], 'Linewidth', 1.5)
set(gca,'ytick',[0:500:peakElevation], 'Linewidth', 1.5)
set(gca,'fontname', 'Arial','FontSize',25, 'fontweight', 'bold')

plot(outlets,modes,'.','color',plotColor,'MarkerSize', markerSize);

ylim([0 peakElevation]);
xlim([0 peakElevation]);


ylabel('Modal elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold','Interpreter', 'none')
xlabel('Outlet elevation [m]', 'Fontsize', 25, 'Fontweight', 'bold', 'Interpreter', 'none')

%%

figName = [groupArea,'_hBench_vs_hChange_allSupercatchments.png']
figFilePath = fullfile(figureOutputFilePath, figName);
fig1 = gcf;
saveas(fig1,figFilePath,outputFigType)

