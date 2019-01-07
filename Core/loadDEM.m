%% loadDEM
% 
% Description
% 
% Given filepath for DEM, create DEM array, Gridobj, and geospatial
% reference object
%
%
%
%%
function  [rawDEM,demGrid, demInfo, geospatialReferenceArray] = loadDEM(filePath, fileName, nanFlag)

[rawDEM,~]= geotiffread(fullfile(filePath,fileName));
  
   [~,geospatialReferenceArray] = geotiffread(fullfile(filePath,fileName));
     rawDEM(rawDEM<0)=NaN;       % define negatives as null
  rawDEM = double(rawDEM);
  rawDEM(rawDEM==nanFlag) = NaN;
  rawDEM(rawDEM==0)=NaN;
  
  demInfo = geotiffinfo(fullfile(filePath, fileName));
  
demGrid =...
    GRIDobj(fullfile(filePath, fileName));
fullDEMSize = size(rawDEM);
fullDEMX = fullDEMSize(:,1);
fullDEMY = fullDEMSize(:,2);
demCoordIndex = [1:(fullDEMX*fullDEMY)]';
[demX,demY]=ind2coord(demGrid, demCoordIndex);
demArrayX = reshape(demX, [fullDEMX,fullDEMY ]); 
demArrayY = reshape(demY, [fullDEMX,fullDEMY ]); 

    if(fullDEMX >1 && fullDEMY>1)
    % Create Grid with properly assigned NaNs

    demGrid = GRIDobj(demArrayX,demArrayY,rawDEM);
    else
    end
end