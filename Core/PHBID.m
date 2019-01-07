%% PHBID
%
% Description
% 
% Generates a structure array of (outlet, mode) pairs.
% 
% Input
% 
% Vector of endpoints for segment of constant modal elevation, vector of 
% outlets, modes, a specified minimum PHB length.
% 
% Output
% 
% Structure array of (outlet, mode) pairs for each PHB step.
% 

function [benchStruct] = PHBID(benchEndPoints, outletList, hypsoPeakList, benchLength)
%% Apply slope test to points between bench edges

sizeBenchArray = size(benchEndPoints);
benchNum = sizeBenchArray(1);
benchStructCount = 0;
%%
for benchCount = 1:benchNum
    
    clear benchListOutlets benchListHypsoPeaks benchoutlets lowSlopeApparentBench benchListNanClear benchList
        
      benchMinOutlet = outletList(benchEndPoints(benchCount, 1));
      benchMaxOutlet = outletList(benchEndPoints(benchCount, 2));
      
      benchListOutlets = find(outletList<= benchMaxOutlet & outletList>= benchMinOutlet);
      benchListHypsoPeaks = hypsoPeakList(benchListOutlets);
      benchOutletElevations = outletList(benchListOutlets);
      %streamPixelPeaks = streamPixelList(benchListOutlets);
      
      benchListLength = length(benchListHypsoPeaks);
      
        %if((max(benchListHypsoPeaks)-min(benchListHypsoPeaks))>minBenchRelief)
            %benchListOutlets = [];
        %end
      
      
      if(~isempty(benchListOutlets) && length(benchListOutlets)>=benchLength)
          benchStructCount = benchStructCount + 1;
          verifiedBenchStruct(benchStructCount) = struct('BenchOutlets', benchOutletElevations,...
              'HypsoPeaks',benchListHypsoPeaks);
      end
      
     
      
end

if exist('verifiedBenchStruct','var')

benchStruct = verifiedBenchStruct;

else
    
   benchStruct = [];
end
