%% Define root paths

groupArea = 'Model/FastScape_10Myr_1mmyr_100yrTimeStep_D_01';
groupName = 'Model'
Defaults;

addpath(phToolsPath)
addpath(topoToolboxFilePath)

%% Define job parameters, include Region, supercatchment

% SupercatchmentNum can be a single value or list, and contains the ID
% numbers for all supercatchments to be passed to PHRun.
Ac=111;
minBenchLength=3;
supercatchmentNum = [9];
outputFigType = 'png';
peakElevationForOutputFig = 2500;

Colors;
plotColor=orange;
%%
RunPH
