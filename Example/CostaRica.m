%% Define root paths

Defaults
addpath(phToolsPath)
addpath(topoToolboxFilePath)

%% Define job parameters, include Region, supercatchment

groupArea = 'CostaRica';
supercatchmentNum = [1:3];
outputFigType = 'png';
peakElevationForOutputFig = 4000;
%%
RunPH
