%% PHBplot
% 
% Description
% 
% For a list of PHBs, plot (h_change, h_bench) pairs
% 
% 

Defaults;
addpath(topoToolboxFilePath); 
figureOutputFilePath = fullfile(phAnalysisFilePath,groupArea, 'Figures')
mkdir(figureOutputFilePath);


for count =supercatchmentNum
    clf;
    clear phbOutletArray
    
peakElevation = 4000;
plotColor = 'b'
markerSize = 25;
supercatchmentNum =count;

baseLevelOutlet =0;
lowELA=3400
highELA = 3800;

nanFlag = -32768


supercatchmentBenchFolderName = ['Supercatchment',num2str(supercatchmentNum)]
supercatchmentFileName = [groupArea,'Supercatchment', num2str(supercatchmentNum),'.tif']
supercatchmentFolderName = [groupArea,'Supercatchment', num2str(supercatchmentNum)]
benchSuperName = ['Supercatchment', num2str(supercatchmentNum)]

filepath = fullfile(phAnalysisFilePath,groupArea,'PHBs','Cusum02_BenchLength3Steps','AllSupercatchments')


allPHBs = dir(filepath);
%%
targetList = 1:length(allPHBs)

for i = 1:length(targetList) 
    
    
    targetPHB = targetList(i);
    fileName = allPHBs(targetPHB).name;
    
    superWordLocation1 = strfind(fileName, 'Supercatchment');
    superWordLocation2 = strfind(fileName, 'Stream');
    inBetweenTextSuper = fileName(superWordLocation1+14:superWordLocation2-1);
    supercatchmentNum = str2double(inBetweenTextSuper); 
   
    if(supercatchmentNum == count)
   
    
    
   
    word1Location = strfind(fileName, 'HypsoPeak');
    word2Location = strfind(fileName, 'Pour');
    inBetweenText = fileName(word1Location+9:word2Location-1);
    hypsoPeak = str2double(inBetweenText);
    
    word3Location = strfind(fileName, 'Elevation');
    word4Location = strfind(fileName, 'Super');
    inBetweenText = fileName(word3Location+9:word4Location-1);
    outlet = str2double(inBetweenText);
    
    deltaH = hypsoPeak-outlet;
   
 
    phbOutletArray(i,1) = outlet;
    phbOutletArray(i,2)=hypsoPeak;   

   
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
figFilePath = fullfile(figureOutputFilePath, figName);
fig1 = gcf;
saveas(fig1,figFilePath,outputFigType)
end


