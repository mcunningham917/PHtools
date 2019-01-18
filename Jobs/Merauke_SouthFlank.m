%% Define root paths

groupArea = 'Merauke';
Defaults;

addpath(phToolsPath)
addpath(topoToolboxFilePath)

%% Define job parameters, include Region, supercatchment

% SupercatchmentNum can be a single value or list, and contains the ID
% numbers for all supercatchments to be passed to PHRun.
Ac=555;
minBenchLength=3;
supercatchmentNum = [26];
outputFigType = 'png';
peakElevationForOutputFig = 6000;
%%
RunPH
