%% longProfileGenerator
% 
% Description
% 
% Generate the flow path for a specified channel. Returns x,y,z for
% profile.
% 
% Input
% 
% Stream network structure array generated with topotoolbox. In structure
% array, channels are indexed by their relationship to confluences, such
% that a new field is created in the structure array at each confluence.
% This function uses the stream network structure array to link together
% all channel links below a specified link. 
% 
% 
% 
% 



function [x,y,z,d,ix] = longProfileGenerator(streamNetworkStruct, streamNetworkStreamObj, demGrid, streamNum, AreaThresh)
% Bring in DEM
demArray = demGrid.Z;

fullDEMSize = size(demArray);
fullDEMX = fullDEMSize(:,1);
fullDEMY = fullDEMSize(:,2);
demCoordIndex = [1:(fullDEMX*fullDEMY)]';
[demX,demY]=ind2coord(demGrid, demCoordIndex);
demArrayX = reshape(demX, [fullDEMX,fullDEMY ]); 
demArrayY = reshape(demY, [fullDEMX,fullDEMY ]); 


% Use supplied stream network structure array to specify channel link
targetTrib = streamNetworkStruct(streamNum);

% Find each downstream link for specified channel link
for tribCount = 1:10
    if(~isnan(targetTrib.tribtoIX) && targetTrib.tribtoIX>0)
        targetStreamTribs(tribCount,1) = targetTrib.tribtoIX;
        targetTribNum = targetTrib.tribtoIX;
        targetTrib = streamNetworkStruct(targetTribNum);
        tribCount = tribCount +1;

    else
        break
   end
end

% Check that stream profile has been built, and write structure array with
% x,y for each channel pixel, organized by channel link

if(exist('targetStreamTribs','var'))
allTribFields = [streamNetworkStruct(streamNum).IX; targetStreamTribs(:,1)];
else 
    targetTribNum = targetTrib.IX;
         allTribFields = targetTribNum;
end

%Use each channel link to write x,y coordinates for entire channel

for allTribsCount = 1:length(allTribFields)
    
        clear tribX tribY tribIX tribD tribZ
        
        thisTrib = allTribFields(allTribsCount);
        tribX=streamNetworkStruct(thisTrib).X';
        tribY=streamNetworkStruct(thisTrib).Y';
        
        tribX = tribX(~isnan(tribX));
        tribY = tribY(~isnan(tribY));
        tribIX = coord2ind(demGrid,tribX,tribY);
        tribIX = tribIX(1:(length(tribIX)));
        [~,streamObjIX] = intersect(streamNetworkStreamObj.IXgrid, tribIX);
        tribD = streamNetworkStreamObj.distance(streamObjIX);
        distIX = streamNetworkStreamObj.IXgrid(streamObjIX);
        distArray = [distIX,tribD];
        tribDSort = sort(distArray);
        tribDFinal = tribDSort(:,2);
        [~, tribSort] = sort(tribIX);
        tribIX = tribIX(tribSort);
        
        tribX = tribX(tribSort);
        tribY = tribY(tribSort);
        tribZ = demArray(tribIX);
       
       if allTribsCount ==1
           streamProfileStruct(allTribsCount) = struct('XCoords',tribX,'YCoords',tribY, 'Elevation', tribZ,...
               'Distance', tribDFinal, 'Index', tribIX);
       else
          
           upstreamDist = vertcat(streamProfileStruct.Distance);
  
               upstreamTribMinDist = min(upstreamDist);
               excessDistance=tribD-upstreamTribMinDist;
               
               [upstreamNeighborTrib,~ ] = find(excessDistance>0);

               excessStreamTribs = struct('ExcessTribsIX', upstreamNeighborTrib);
               
                   tribX(upstreamNeighborTrib)=[];
                   tribY(upstreamNeighborTrib)=[];
                   tribZ(upstreamNeighborTrib)=[];
                   tribD(upstreamNeighborTrib)=[];
                   tribIX(upstreamNeighborTrib)=[];
               
                   streamProfileStruct(allTribsCount) = struct('XCoords',tribX,'YCoords',tribY, 'Elevation', tribZ,...
                       'Distance', tribD,'Index', tribIX);

           end

end

% Concatenate x,y,z for every link in channel chain

x = vertcat(streamProfileStruct.XCoords);
y = vertcat(streamProfileStruct.YCoords);
z = vertcat(streamProfileStruct.Elevation);
d = vertcat(streamProfileStruct.Distance);
fullStreamIX = vertcat(streamProfileStruct.Index);
[ix,~,~] = coord2ind(demArrayX,demArrayY,x,y);

%%
 
% d = zeros(length(x),1);
% for i = 2:(length(x)-1)
%     d(i) = sqrt(((x(i)-x(i+1))^2) + ((y(i)-y(i+1))^2));
% end

%d = cumsum(d);
%d = flipud(d);

end

       