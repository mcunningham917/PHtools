%% ExtractSubcatchmentDEM 
% 
% Description
% 
% Creates DEM array of the drainage basin that drains to a supplied pour
% point.
% 

function  finalDrainageBasinArray = ExtractSubcatchmentDEM(flowDirectionObj,dem, x, y)

    finalDrainageBasinStruct =...
        drainagebasins(flowDirectionObj,x, y);
     
    finalDrainageBasinArray=finalDrainageBasinStruct.Z;
    finalBasinMask=finalDrainageBasinStruct.Z;
    finalDrainageBasinArrayIndex = find(finalBasinMask==1);
    finalDrainageBasinArray = double(finalDrainageBasinArray);
    finalDrainageBasinArray(finalDrainageBasinArrayIndex) = dem(finalDrainageBasinArrayIndex);
    finalDrainageBasinArray(finalDrainageBasinArray==0)= NaN;
end