function [totalArea, totalLength, areaAll, lengthAll] = genScore_branch(regionx)


regionxROI = bwlabel(regionx);
regionxROIidx = label2idx(regionxROI);
regionx_skel = bwskel(logical(regionx),'MinBranchLength',15);
totalArea = sum(regionx(:));

totalLength = sum(regionx_skel(:));
areaAll = zeros(length(regionxROIidx),1);
lengthAll = zeros(length(regionxROIidx),1);
for i = 1:length(regionxROIidx)
    areaAll(i) = sum(regionx(regionxROIidx{i}));
    lengthAll(i) = sum(regionx_skel(regionxROIidx{i}));

end

end