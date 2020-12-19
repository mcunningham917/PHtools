%% Define root paths

groupArea = 'Taiwan_NS';
groupName = 'Taiwan_NS';
Defaults;

addpath(phToolsPath)
addpath(topoToolboxFilePath)

%% Define job parameters, include Region, supercatchment

% SupercatchmentNum can be a single value or list, and contains the ID
% numbers for all supercatchments to be passed to PHRun.
Ac=555;
minBenchLength=3;
supercatchmentNum = [1:23];

Colors;
plotColor = blue;
outputFigType = 'png';
peakElevationForOutputFig = 4000;
%%
%RunPH
