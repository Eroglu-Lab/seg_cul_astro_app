function newSynIdAllfin = mito_detect(imregion1, xxshift1, yyshift1, smoFactor, zthresca3)

if(smoFactor == 0)
imregion1G = imregion1;
else
imregion1G = imgaussfilt(imregion1,smoFactor);
end

localmean5 = imboxfilt(imregion1,[5,5]);
localmean3 = imboxfilt(imregion1,[3,3]);    
imregion2 = imregion1(3:end-2, 3:end-2);
localmean2 = 25.*localmean5(3:end-2, 3:end-2) - 9.*localmean3(3:end-2, 3:end-2);
localmean2 = 25.*localmean5 - 9.*localmean3;
localmean2 = localmean2./(25-9);
%
diff = (imregion1 - localmean3);
diff2 = (imregion1 - localmean2);
% diff2 = [diff(diff(:) <0); diff(diff(:) == 1); -diff(diff(:) < 0)];
svar = sqrt(var(diff(:)));
[h,w] = size(imregion1);
[lenx, leny] = size(imregion1);
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
%         signalNei3 = imregion1(idxtmpnei3);
        if(~isempty(signal) && ~isempty(signalNei))
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
        if(~isnan(zscorepatch0) && zscorepatch0 > max(zscoreMap(idxtmp)) && length(idxtmp) > 30)
            zscoreMap(idxtmp) = zscorepatch0;
        end
    end
% 
end
newSynIdAllfin = double(zscoreMap > zthresca3);
end