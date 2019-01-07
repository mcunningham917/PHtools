function  finalDrainageBasinGRID = ExtractSubcatchmentGRID(dem, flowDirectionObj, x, y)

    finalDrainageBasinStruct =...
        drainagebasins(flowDirectionObj,x, y);
    
    
     
    finalDrainageBasinArray=finalDrainageBasinStruct.Z;
    finalBasinMask=finalDrainageBasinStruct.Z;
    finalDrainageBasinArrayIndex = find(finalBasinMask==1);
    finalDrainageBasinArray = double(finalDrainageBasinArray);
    finalDrainageBasinArray(finalDrainageBasinArrayIndex) = dem(finalDrainageBasinArrayIndex);
    finalDrainageBasinArray(finalDrainageBasinArray==0)= NaN;
    
    
    fullDEMSize = size(finalDrainageBasinArray);
    fullDEMX = fullDEMSize(:,1);
    fullDEMY = fullDEMSize(:,2);
    demCoordIndex = [1:(fullDEMX*fullDEMY)]';
    [demX,demY]=ind2coord(finalDrainageBasinStruct, demCoordIndex);
    demArrayX = reshape(demX, [fullDEMX,fullDEMY ]); 
    demArrayY = reshape(demY, [fullDEMX,fullDEMY ]); 
    
    finalDrainageBasinGRID = GRIDobj(demArrayX,demArrayY, finalDrainageBasinArray);
end