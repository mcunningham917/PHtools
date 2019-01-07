%% TrackToDivide 
% 
% Description
% 
% Delineates a flow path by following maximum upstream accumlation area.
% Ouputs z and IX for each pixel along flow path
% 
% 


function [aboveHeadZ, upstreamNodes]= TrackToDivide(demGrid, flowDirectionObj, flowAccumulationStruct)

    givingNodes = flowDirectionObj.ix;
    receivingNodes = flowDirectionObj.ixc;

    givingNodes = flipud(givingNodes);
    receivingNodes = flipud(receivingNodes);

    flowAccumulationGiving = flowAccumulationStruct.Z(givingNodes);
    flowAccumulationReceiving = flowAccumulationStruct.Z(receivingNodes);

    
    j = 0;
    receiver = receivingNodes(1);

    for i = 1:length(receivingNodes)

        numGivers = find(receivingNodes==receiver);

        if(~isempty(numGivers))

       [~, upstreamIndex] = max(flowAccumulationGiving(numGivers));
       upstreamIndex = numGivers(upstreamIndex);
       j = j+1;
       upstreamAreaNode(j,1) = givingNodes(upstreamIndex);
       receiver = givingNodes(upstreamIndex);

       i = i+max(numGivers);

        else
        end
    end

    % Use the list of new max flow accumulations to extract upper catchments
    
    aboveHeadZ = demGrid.Z(upstreamAreaNode);
    upstreamNodes = upstreamAreaNode;
end