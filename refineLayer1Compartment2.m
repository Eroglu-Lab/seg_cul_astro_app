function [refinedlayer2] = refineLayer1Compartment2(layer1, layer2, layer34)
    [lenx, leny] = size(layer1);
    se = strel('disk', 1);
    layer2ROI = bwlabel(layer2);
    layer2ROIidx = label2idx(layer2ROI);
    layer2ROIidx = layer2ROIidx(:);
    layer1 = imdilate(layer1 , se);
    selectedID = layer2ROI(layer1(:) == 1);
    selectedID = unique(selectedID);
    selectedID(selectedID == 0) = [];
    refinedlayer2 = zeros(lenx, leny);
    refinedlayer2(cell2mat(layer2ROIidx(selectedID))) = 1;

end