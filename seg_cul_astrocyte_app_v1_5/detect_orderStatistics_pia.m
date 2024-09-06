function [newSynIdAllfin,zscoreMap] = detect_orderStatistics_pia(imregion1, imregion1G, svar,zthresca3)


[h,w] = size(imregion1);
iters1 = 1; %11x11x5
xxshift1 = zeros(2*iters1+1, 2*iters1+1);
yyshift1 = zeros(2*iters1+1, 2*iters1+1);
for i = -iters1:iters1
    for j = -iters1:iters1
        xxshift1(i+iters1+1,j+iters1+1) = i;
        yyshift1(i+iters1+1,j+iters1+1) = j;
    end
end
zscoreMap = zeros(h, w);
zscoreThres = zthresca3;
[lenx, leny] = size(imregion1G);
detectionRegion = zeros(size(imregion1G));
detectionZscoreRegion = zeros(size(imregion1G));
detectionRegionZthres = zeros(size(imregion1G));
detectionRegionZcombined = zeros(size(imregion1G));
% localVariance = [];
% count = 1;
for i = (floor(max(imregion1G(:)))+1):-1:1
    mask = double(imregion1G > i);
    maskroi = bwlabel(mask);
    maskroiIDx = label2idx(maskroi);
    lengthx = cellfun(@length, maskroiIDx);
    maskroiIDx(lengthx < 20) = [];
    for j = 1:length(maskroiIDx)
        %%%%%%%patch0
        idxtmp = maskroiIDx{j};
        idxtmpneiL1 = regionGrowxx(idxtmp,1, lenx,leny,xxshift1, yyshift1);
%         idxtmpneiL3 = regionGrowxx(idxtmp,3, lenx,leny,xxshift1, yyshift1);
        idxtmpnei = setdiff(idxtmpneiL1, idxtmp);
%         idxtmpnei3 = setdiff(idxtmpneiL3, idxtmp);
        idxtmpneiInner = regionGrowxx(idxtmpnei,1, lenx,leny,xxshift1, yyshift1);
        idxtmpInner = intersect(idxtmpneiInner, idxtmp);
        signal = imregion1(idxtmpInner);
        signalNei = imregion1(idxtmpnei);
        signal(isnan(signal)) = [];
        signalNei(isnan(signalNei)) = [];
%         signalNei3 = imregion1(idxtmpnei3);
%         if(length(signal) > 5 && sum(isnan(signalNei)) == 0)
        if(length(signal) > 5 && length(signalNei) > 5)
            [mutmp, sigmatmp] = ksegments_orderstatistics_fin(signal, signalNei);
            meanDiff = mean(signal) - mean(signalNei);
            zscorepatch0 = (meanDiff - mutmp.*svar)./(svar.*sigmatmp);
        else
            zscorepatch0 = nan;
        end
        %%%%%pick the largest from the two
        idxtmpNewLabel = idxtmp(detectionRegion(idxtmp) == 0);  
        
        detectionZscoreRegion(idxtmpNewLabel) = zscorepatch0;
        detectionRegion(idxtmp) = 1;
        idx_lenx = rem(idxtmp, lenx);
        idx_leny = ceil(idxtmp./lenx);
        ratiox = length(idxtmp)/((max(idx_lenx) - min(idx_lenx))*(max(idx_leny) - min(idx_leny)));
        if(~isnan(zscorepatch0) && zscorepatch0 > max(zscoreMap(idxtmp)) && length(idxtmp) > 5 && length(idxtmp) < 2000 && ratiox > 0.001)
            zscoreMap(idxtmp) = zscorepatch0;
        end
    end
% 
end
newSynIdAllfin = double(zscoreMap > zthresca3);
end