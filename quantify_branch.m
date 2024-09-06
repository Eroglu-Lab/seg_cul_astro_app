function quantify_branch(primaryB, layer2B, layer34B, terminalB, folderPath)
primaryB(primaryB < 0) = 0;
layer2B(layer2B < 0) = 0;
layer34B(layer34B < 0) = 0;
terminalB(terminalB < 0) = 0;
% output the area of the branches for each type 
primary_area_ROI = bwlabel(primaryB);
primary_area_ROI_idx = label2idx(primary_area_ROI);
area_all = cellfun(@length,primary_area_ROI_idx);
[lenx, leny] = size(primaryB);
id0 = find(primaryB(:) == 1);
[id0x, id0y] = ind2sub([lenx, leny], id0);
skel_B = bwskel(logical(primaryB),'MinBranchLength',15);
maskDist = bwdist(1 - primaryB);
maskDist = imgaussfilt(maskDist,2).*primaryB;
[centerx, centery] = find(maskDist == max(maskDist(:)));
radius = maskDist(centerx, centery);
% remove the soma region to obtain the branches (the branches should still be one connected component)
distAll = sqrt((id0x - centerx).^2 + (id0y - centery).^2);
idrm = id0(distAll <= radius);
primaryB(idrm) = 0;
skel_B(idrm) = 0;
lenAll_primary = sum(skel_B(:));
lenAll_cell = zeros(length(primary_area_ROI_idx),1);
for i = 1:length(lenAll_cell)
    lenAll_cell(i) = sum(skel_B(primary_area_ROI_idx{i}));
end

fileID = fopen([folderPath,'primary_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', length(area_all));
fprintf(fileID,'total area %f\n', sum(area_all));
fprintf(fileID,'total length %f\n', lenAll_primary);
fprintf(fileID, ['ID ' 'area ' 'length\n']);
for i = 1:length(area_all)
    fprintf(fileID, '%d %f %f\n', i,area_all(i), lenAll_cell(i));
end
fclose(fileID);

[totalArea, totalLength, areaAll, lengthAll] = genScore_branch(layer2B);

fileID = fopen([folderPath,'secondary_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', length(areaAll));
fprintf(fileID,'total area %f\n',totalArea);
fprintf(fileID,'total length %f\n', totalLength);
fprintf(fileID, ['ID ' 'area ' 'length\n']);
for i = 1:length(areaAll)
    fprintf(fileID, '%d %f %f\n', i,areaAll(i), lengthAll(i));
end
fclose(fileID);




[totalArea, totalLength, areaAll, lengthAll] = genScore_branch(layer34B);

fileID = fopen([folderPath,'fine_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', length(areaAll));
fprintf(fileID,'total area %f\n',totalArea);
fprintf(fileID,'total length %f\n', totalLength);
fprintf(fileID, ['ID ' 'area ' 'length\n']);
for i = 1:length(areaAll)
    fprintf(fileID, '%d %f %f\n', i,areaAll(i), lengthAll(i));
end
fclose(fileID);





[totalArea, totalLength, areaAll, lengthAll] = genScore_branch(terminalB);
fileID = fopen([folderPath,'terminal_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', length(areaAll));
fprintf(fileID,'total area %f\n',totalArea);
fprintf(fileID,'total length %f\n', totalLength);
fprintf(fileID, ['ID ' 'area ' 'length\n']);
for i = 1:length(areaAll)
    fprintf(fileID, '%d %f %f\n', i,areaAll(i), lengthAll(i));
end
fclose(fileID);
















end