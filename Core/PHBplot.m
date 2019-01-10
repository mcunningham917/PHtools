%% PHBplot
% 
% Description
% 
% For a list of PHBs, plot (h_change, h_bench) pairs
% 
% 

figureOutputFilePath = fullfile(phAnalysisFilePath,groupArea, 'Figures');
supercatchmentFigureOutputFilePath = fullfile(phAnalysisFilePath,groupArea, 'Figures','SupercatchmentPHBs');
mkdir(figureOutputFilePath);
mkdir(supercatchmentFigureOutputFilePath);


for count =supercatchmentNum
    clf;
    clear phbOutletArray
   
 
%% Plot options
peakElevation = peakElevationForOutputFig;
plotColor = 'b'
markerSize = 25;

baseLevelOutlet =0;


%% Input and output files

nanFlag = -32768;
supercatchmentBenchFolderName = ['Supercatchment',num2str(count)]
supercatchmentFileName = [groupArea,'Supercatchment', num2str(count),'.tif']
supercatchmentFolderName = [groupArea,'Supercatchment', num2str(count)]


filepath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchmentsTxt')
roiFileName = [groupArea, '_allPHBs.txt']


allPHBs = dir(filepath);
%%
targetList = 1:length(allPHBs)

for i = 1:length(targetList) 
    
    
    targetFolder = targetList(i);
    fileName = allPHBs(targetFolder).name;
    
    if(~isempty(strfind(fileName, 'Supercatchment')))
    
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, 'PHBs');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    superNum = str2double(inBetweenTextSuper); 
    
    benchSuperName = ['Supercatchment', num2str(superNum),'PHBs']
    supercatchmentPHBFileName = ['Supercatchment',num2str(superNum),'_allOutletModePairs.txt'];
    phbOutletArray = dlmread(fullfile(filepath, benchSuperName, supercatchmentPHBFileName),'\t',2);
    
    else 
        i=i+1;
    end
    
   
    
  
    
end


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

figure(supercatchmentNum)
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



%% Save figure

figName = [supercatchmentBenchFolderName,'_hBench_vs_hChange.png']
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
