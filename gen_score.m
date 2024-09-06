function [num_primary,size_primary_mean,circu_primary_mean, size_primary_all,circu_primary_all] = gen_score(mitox, primaryB, xxshift, yyshift)

    mito_primary = mitox.*primaryB;
    mito_primaryROI = bwlabel(mito_primary);
    mito_primaryROIidx = label2idx(mito_primaryROI);
    mito_primaryROIidx(cellfun(@length,mito_primaryROIidx)<5) = [];
    num_primary = length(mito_primaryROIidx);
    size_primary_mean = mean(cellfun(@length, mito_primaryROIidx));
    size_primary_all = cellfun(@length, mito_primaryROIidx);
    circu_primary_all = zeros(length(mito_primaryROIidx),1);
    [lenx, leny] = size(mitox);
    for i = 1:length(mito_primaryROIidx)
        tmpID=  mito_primaryROIidx{i};
        tmpIDss = regionGrowxx(tmpID, 1, lenx, leny,xxshift, yyshift);
        cc = length(setdiff(tmpIDss, tmpID));
        circu_primary_all(i) = 4*pi*size_primary_all(i)/(cc^2);
    end

    circu_primary_mean = mean(circu_primary_all(~isinf(circu_primary_all)));



end