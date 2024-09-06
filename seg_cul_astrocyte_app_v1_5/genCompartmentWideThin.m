function [widePart, thinPart, terminalBranchOut,branchPointRegion, branch_s_idx2, widthMean,maskWS2] = genCompartmentWideThin(mask, maskskel,maskDist, thres0,xxshift, yyshift)
    [lenx, leny] = size(mask);
    se = strel('disk', 1);
    terminalBranch = zeros(lenx, leny);
    maskDist = imgaussfilt(maskDist,2).*mask;
    maskSkelDist = maskskel .* maskDist;
    maskBP = bwmorph(maskskel,'branchpoints');
    [centerx, centery] = find(maskDist == max(maskDist(:)));
    maskCenter = zeros(lenx, leny);
    maskCenter(centerx, centery ) =1;
    dist2Center = bwdistgeodesic(logical(mask), logical(maskCenter));
    [branch_s, maskBP] = findRealBranch(maskskel,maskBP);
%     branch_s = double(maskskel - maskBP > 0);
    branch_s_ROI = bwlabel(branch_s);
    branch_s_idx = label2idx(branch_s_ROI);
    branch_s_idx = branch_s_idx(:);
    
    se2 = strel('disk', 3);
    mask_astro2 = imopen(mask, se2);
    masksspSkel2 = bwskel(logical(mask_astro2),'MinBranchLength',15);
    maskBP2 = bwmorph(masksspSkel2,'branchpoints');
    [branch_s_2, maskBP2] = findRealBranch(masksspSkel2,maskBP2);
    branch_s_ROI2 = bwlabel(branch_s_2);
    branch_s_idx2 = label2idx(branch_s_ROI2);
    branch_s_idx2 = branch_s_idx2(:);
    maskDist_2 = bwdist(1 - mask_astro2);
    branchPointRegion = recoverBPregion(maskBP2, maskDist_2, xxshift,yyshift);
    widthMean = cellfun(@(c) mean(maskDist_2(c)), branch_s_idx2);
    idWide = find(widthMean > thres0); % the ids ralated to the large branches
    idWide0 = unique(branch_s_ROI(cell2mat(branch_s_idx2(idWide))));
    idThin = setdiff([1:length(branch_s_idx)], idWide0);
    maskDist_3_r = max(maskDist_2(:)) - maskDist_2;
    maskDist_3_r = imimposemin(maskDist_3_r, branch_s_2);
    maskWS2 = double(watershed(maskDist_3_r)).*mask_astro2;
    maskWSidx2 = label2idx(maskWS2);
    maskWSidx2 = maskWSidx2(:);
    maskWideLabel = unique(maskWS2(cell2mat(branch_s_idx2(idWide))));
    maskWideLabel(maskWideLabel == 0) = [];
    widePart = zeros(size(mask));
    widePart(cell2mat(maskWSidx2(maskWideLabel))) = 1;
    
    
    
    maskDist2 = max(maskDist(:)) - maskDist;
    maskDist2 = imimposemin(maskDist2, branch_s);
    maskWS = double(watershed(maskDist2)).*mask;
    maskWSidx = label2idx(maskWS);
    maskWSidx = maskWSidx(:);
%     thinPartROI = bwlabel(thinPart);
    widePart = imclose(widePart, se);
    thinPart = mask - widePart;
    %terminal points should be shared by only one branch
    maskskel2 = bwskel(logical(mask), 'MinBranchLength',15);
    endPoints = bwmorph(maskskel2, 'endpoints');
    endPoints = endPoints.*(1 - maskBP);
    endPointsx = imdilate(endPoints, se);
    maskskel2ROI = maskskel2 - endPointsx;
    maskskel2ROI(maskskel2ROI < 0) = 0;
    maskskel2ROIx = bwlabel(maskskel2ROI);
    maskskel2ROIidx = label2idx(maskskel2ROIx);
    se = strel('disk',2);
    endPointsx = imdilate(endPointsx, se);
    endPointsxidx = label2idx(bwlabel(endPointsx));
    endPointsxidx = endPointsxidx(:);
    rmEndPoint = [];
    for i = 1:length(endPointsxidx)
        labeltmp = unique(maskskel2ROIx(endPointsxidx{i}));
        labeltmp(labeltmp ==0) = [];
        if(length(labeltmp) >1)
            rmEndPoint = [rmEndPoint;i];
        end
    end
    % remove end points that are shared by multiple branches
    endPoints(cell2mat(endPointsxidx(rmEndPoint))) = 0;
    endPointsROI = bwlabel(endPointsx);
    endPointsROIidx = label2idx(endPointsROI);
    countsOfendPoints = zeros(length(endPointsROIidx),1);
    
    for i = 1:length(branch_s_idx)
        tmpID = branch_s_idx{i};
        tmpIDGrowed = regionGrowxx(tmpID, 1, lenx, leny, xxshift, yyshift);
        labelEndPoints = unique(endPointsROI(tmpIDGrowed));
        labelEndPoints(labelEndPoints == 0) = [];
        if(~isempty(labelEndPoints))
            countsOfendPoints(labelEndPoints) = countsOfendPoints(labelEndPoints) +1 ;
        end
    end
    endPoints(cell2mat(endPointsxidx(countsOfendPoints > 1))) = 0;
    endPointsBranchLabel = unique(maskWS(endPoints(:) == 1));
    endPointsBranchLabel(endPointsBranchLabel == 0) = [];
    terminalBranch(cell2mat(maskWSidx(endPointsBranchLabel))) = 1;
    terminalBranchROI = bwlabel(terminalBranch);
    terminalBranchROIidx = label2idx(terminalBranchROI);
    terminalBranchROIidx = terminalBranchROIidx(:);
    rmID3 = [];
    % the distance to the center should be the farthest comparing to all
    % its neighbor regions
    for i = 1:length(terminalBranchROIidx)
        curMaxDist = max(dist2Center(terminalBranchROIidx{i}));
        growedID = regionGrowMatrix(terminalBranchROIidx{i}, 10, lenx, leny);
        growedID(terminalBranch(growedID) == 1) = [];
        surroundingMaxDist = dist2Center(growedID);
        surroundingMaxDist(isinf(surroundingMaxDist)) = [];
        if(max(surroundingMaxDist) > curMaxDist)
            rmID3  =[rmID3; i];
        end
    end
    terminalBranchROIidx(rmID3) = [];
%     terminalBranchOut = terminalBranch;
%     terminalBranchOut((cell2mat(terminalBranchROIidx(rmID3)))) = 0;
%     terminalBranchOutROI = bwlabel(terminalBranchOut);
%     terminalBranchOutROIidx = label2idx(terminalBranchOutROI);
    maxDist2center = cellfun(@(c) max(dist2Center(c)), terminalBranchROIidx);
    dist2centerPool = dist2Center(mask(:) == 1);
    mux = mean(dist2centerPool(~isinf(dist2centerPool)));
    sigmax = std(dist2centerPool(~isinf(dist2centerPool)));
    % adhoc remove certain terminal points that are too close to the center
    % also remove any terminal points that are at the outter part of the
    % region
    terminalBranchROIidx((maxDist2center - mux)/ sigmax < 0 | isinf(maxDist2center)) = []; 
    terminalBranchOut = zeros(lenx, leny);
    terminalBranchOut(cell2mat(terminalBranchROIidx)) = 1;
end