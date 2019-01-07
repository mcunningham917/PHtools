%% SubcatchmentWrite
% 
% Description
% 
% Given the array of a catchment, write out a geotiff.
% 
% 


function subcatchmentWrite(subcatchmentDEM,geospatialReferenceArray, geotiffInfo, nanFlag, outputFile)

    [drainageBasinRow, drainageBasinCol] = find(~isnan(subcatchmentDEM));
 
    topDEM=min(drainageBasinRow);
    bottomDEM=max(drainageBasinRow);
    leftDEM=min(drainageBasinCol);
    rightDEM=max(drainageBasinCol);
    
    subcatchmentDEM(isnan(subcatchmentDEM))= nanFlag;
    subcatchmentDEM=double(subcatchmentDEM);
    
    subImage = subcatchmentDEM(topDEM:bottomDEM, leftDEM:rightDEM, :);

    
    xi = [leftDEM - .5, rightDEM + .5];
    yi = [topDEM - .5, bottomDEM + .5];
    [xlimits, ylimits] = intrinsicToWorld(geospatialReferenceArray, xi, yi);
    subR = geospatialReferenceArray;
    subR.RasterSize = size(subImage);
    subR.XLimWorld = sort(xlimits);
    subR.YLimWorld = sort(ylimits);


    geotiffwrite(outputFile, subImage, subR, 'GeoKeyDirectoryTag', geotiffInfo.GeoTIFFTags.GeoKeyDirectoryTag)

end