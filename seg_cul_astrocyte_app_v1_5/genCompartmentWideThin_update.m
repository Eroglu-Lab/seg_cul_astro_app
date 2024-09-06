function [widePart, thinPart] = genCompartmentWideThin_update(mask, maskskel,maskDist, thres0)
    [lenx, leny] = size(mask);
    terminalBranch = zeros(lenx, leny);
    maskDist = imgaussfilt(maskDist,2).*mask;
    maskSkelDist = maskskel .* maskDist;
    maskBP = bwmorph(maskskel,'branchpoints');
    se = strel('disk', 1);
    % some branch points are not real
    [branch_s, maskBP] = findRealBranch(maskskel,maskBP);
%     branch_s = double(maskskel - maskBP > 0);
    branch_s_ROI = bwlabel(branch_s);
    branch_s_idx = label2idx(branch_s_ROI);
    branch_s_idx = branch_s_idx(:);
    widthMean = cellfun(@(c) mean(maskSkelDist(c)), branch_s_idx);
    idThin = find(widthMean <= thres0);
    idWide = find(widthMean > thres0);
    maskWide = zeros(size(mask));
    maskThin = zeros(size(mask));
    maskWide(cell2mat(branch_s_idx(idWide))) = 1;
    maskThin(cell2mat(branch_s_idx(idThin))) = 1;
    maskDist2 = max(maskDist(:)) - maskDist;
    maskDist2 = imimposemin(maskDist2, (maskWide + maskThin));
    maskWS = double(watershed(maskDist2)).*mask;
    maskWSidx = label2idx(maskWS);
    maskWSidx = maskWSidx(:);
    maskWideLabel = unique(maskWS(cell2mat(branch_s_idx(idWide))));
    maskWideLabel(maskWideLabel == 0) = [];
    maskThinLabel = unique(maskWS(cell2mat(branch_s_idx(idThin))));
    maskThinLabel(maskThinLabel ==0) = [];
    widePart = zeros(size(mask));
    thinPart = zeros(size(mask));
    widePart(cell2mat(maskWSidx(maskWideLabel))) = 1;
    thinPart(cell2mat(maskWSidx(maskThinLabel))) = 1;
%     thinPartROI = bwlabel(thinPart);
    widePart = imclose(widePart, se);
    thinPart = imclose(thinPart, se);
end