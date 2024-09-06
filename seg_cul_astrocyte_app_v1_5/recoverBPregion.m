function branchPointRegion = recoverBPregion(maskBP2, maskDist_2, xxshift, yyshift)
    
    maskBP2ROI = bwlabel(maskBP2);
    maskBP2ROIidx = label2idx(maskBP2ROI);
    radiusx = cellfun(@(c) max(maskDist_2(c)), maskBP2ROIidx);
    [lenx, leny] =size(maskBP2);
    branchPointRegion = zeros(lenx, leny);
    for i = 1:length(maskBP2ROIidx)
        outputID = regionGrowxx(maskBP2ROIidx{i},round(radiusx(1)), lenx,leny,xxshift, yyshift);
        branchPointRegion(outputID) = 1;
    end


end