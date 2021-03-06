%% Define root paths

Defaults;
addpath(phToolsPath)
addpath(topoToolboxFilePath)

%% Define job parameters, include Region, supercatchment

% SupercatchmentNum can be a single value or list, and contains the ID
% numbers for all supercatchments to be passed to PHRun.

groupArea = 'Merauke';
supercatchmentNum = [9,16];
outputFigType = 'png';
peakElevationForOutputFig = 4000;
%%
RunPH
