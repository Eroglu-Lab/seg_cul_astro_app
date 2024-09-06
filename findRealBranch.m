function [newBranch, newMaskBP] = findRealBranch(maskskel,maskBP)
% some branches are not real branch because removing the branch point from
% the branch will only output two branches rather than three
    se = strel('disk', 1);
    [lenx, leny] = size(maskskel);
    maskBP = imdilate(maskBP, se);
    branch_s = maskskel - maskBP;
    branch_s(branch_s(:) < 0) = 0;
    branch_s_ROI = bwlabel(branch_s);
    branch_s_idx = label2idx(branch_s_ROI);
%     branch_s_idx = branch_s_idx(:);
    se = strel('disk', 2);
    maskBP2 = imdilate(maskBP, se);
    maskBP2ROI = bwlabel(maskBP2);
    maskBP2ROIidx = label2idx(maskBP2ROI);
    countOfBranch = zeros(length(maskBP2ROIidx), 1);
    for i = 1:length(maskBP2ROIidx)
        branchesLabel = unique(branch_s_ROI(maskBP2ROIidx{i}));
        branchesLabel(branchesLabel == 0) = [];
        countOfBranch(i) = length(branchesLabel);
    end
    maskBP2ROIidx(countOfBranch < 3) = [];
    maskBP2ROIidx = maskBP2ROIidx(:);
    maskBP2 = zeros(lenx, leny);
    maskBP2(cell2mat(maskBP2ROIidx)) = 1;
    newMaskBP = maskBP2;
    newBranch = maskskel - newMaskBP;
    newBranch(newBranch(:) < 0) = 0;
end