%% ProminentHypsoID
% 
% Description
% 
% Identifies endpoints of steps in elevation of mode using CUSUM. 
% 
% Input
% 
% Vector of modes.
% 
% Output
% 
% 2xn array, where each row is the bounding x-position for each (outlet,
% mode) step.
%
%

function [X] = ProminentHypsoID(hypsoPeakList, cusumThresh, benchLengthNum)

%% Implement CUSUM

hypsoPeakListNormal = hypsoPeakList./(max(hypsoPeakList));
cumulativeHypsoPeaks = cumsum(hypsoPeakListNormal);
cuSumRange = max(cumulativeHypsoPeaks);
linearTrendLine = linspace(0, cuSumRange, length(hypsoPeakList));
linearTrendLine = linearTrendLine';
detrendedCuSum = cumulativeHypsoPeaks - linearTrendLine;

detrendedCuSumSmooth = medfilt1(detrendedCuSum,3);

firstDiffDetrendedCuSumSmooth = diff(detrendedCuSumSmooth); % First derivative of detrended HypsoPeak list
secondDiffDetrendedCuSumSmooth = diff(firstDiffDetrendedCuSumSmooth); % Second derivative of detrended HypsoPeak list

benchEdges = find(abs(secondDiffDetrendedCuSumSmooth)>cusumThresh);
benchEdges = benchEdges + 1;

if length(benchEdges) ==1
    benchEdges = [1, benchEdges];
end

if(length(benchEdges)==0)
    maxDiff = max(abs(secondDiffDetrendedCuSumSmooth))
    if(maxDiff<cusumThresh)
        benchEdges = [1, length(hypsoPeakList)]
    else
        benchEdges = 0;
    end
end
    

%% Define bench by number of modes at the same elevation

distanceBetweenBenchPoints = diff(benchEdges);

benchEndPointsCount = 0;

for benchTestNum = 1: length(benchEdges)
    
    if(benchTestNum<length(benchEdges))
        
            if(distanceBetweenBenchPoints(benchTestNum)>benchLengthNum)
                benchEndPointsCount = benchEndPointsCount+1;
                benchEndPoints(benchEndPointsCount,1) = (benchEdges(benchTestNum))+1;
                benchEndPoints(benchEndPointsCount,2) = (benchEdges(benchTestNum+1));
            end
    else
        if((length(hypsoPeakListNormal)-benchEdges(benchTestNum))>benchLengthNum)
            benchEndPointsCount = benchEndPointsCount+1;
            benchEndPoints(benchEndPointsCount,1) = (benchEdges(benchTestNum))+1;
            benchEndPoints(benchEndPointsCount,2) = length(hypsoPeakListNormal);
        end
        
    end
    
end

if exist('benchEndPoints','var')

X = benchEndPoints;

else
    
   X = 0;
end