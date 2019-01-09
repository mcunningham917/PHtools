%% Define root paths

Defaults
addpath(phToolsPath)

%% Define job parameters, include Region, supercatchment

groupArea = 'CostaRica';
supercatchmentNum = [9];
outputFigType = 'png';
peakElevationForOutputFig = 4000;
%%
RunPH
