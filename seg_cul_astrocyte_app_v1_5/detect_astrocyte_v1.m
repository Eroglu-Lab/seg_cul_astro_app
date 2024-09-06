function mask_ssp = detect_astrocyte_v1(imregion1, xxshift1, yyshift1, smoFactor, zthresca3, newfolder)
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
for i = (floor(max(imregion1(:)))+1):-1:1
    mask = double(imregion1 > i);
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
        if(~isnan(zscorepatch0) && zscorepatch0 > max(zscoreMap(idxtmp)))
            zscoreMap(idxtmp) = zscorepatch0;
        end
    end
% 
end
newSynIdAllfin = double(zscoreMap > zthresca3);
imwrite(newSynIdAllfin, fullfile(newfolder, 'astro_detection_step1.tif'));

svar = sqrt(var(diff(zscoreMap(:) > zthresca3)));
%% detect scattered points/ regions in the remaining area
se = strel('disk', 1);
zscoreMapAll = zeros(lenx, leny);
zthresca3 = 3;
while(1)
    newSynIdAllfintmp = imdilate(newSynIdAllfin,se);
%     imregion1(newSynIdAllfintmp(:) == 1) = nan;
    imregion1G(newSynIdAllfintmp(:) == 1) = nan;
    imregion1(newSynIdAllfin(:) == 1) = nan;
%     imregion1G(newSynIdAllfin(:) == 1) = nan;
[remainedDetection,zscoreMaptmp] = detect_orderStatistics_pia(imregion1, imregion1G, svar,zthresca3);
zscoreMapAll = zscoreMapAll + zscoreMaptmp;
newSynIdAllfin = newSynIdAllfin + remainedDetection;
if(sum(remainedDetection(:)) == 0)
    break;
end
% imwrite(newSynIdAllfin + remainedDetection, 'F:\eroglu_lab\pia\0025_astrocyte_detection_test_2_round.tif')


end

imwrite(newSynIdAllfin, fullfile(newfolder, 'astro_detection_step2.tif'));


se0 = strel('disk', 1);
newSynIdAllfin = imclose(newSynIdAllfin, se0);
% newSynIdAllfin = im2double(imread('F:\eroglu_lab\pia\0025_astrocyte_detection_initial_components_v2.tif'));
%% examine the shortest path between each pair of the component
ccROI = bwlabel(newSynIdAllfin);
ccROIidx = label2idx(ccROI);
ccROIidxGrow = cell(length(ccROIidx),1);
for i = 1:length(ccROIidxGrow)
    
   ccROIidxGrow{i} = regionGrowMatrix(ccROIidx{i},20, lenx, leny); 
    
    
end
graphNode = [];
scoreMap = 1./(imregion1G.^2);
ccROIidx = ccROIidx(:);
lengthCEll = cellfun(@length, ccROIidx);
rootID = find(lengthCEll == max(lengthCEll));
for i = 1:(length(ccROIidxGrow) - 1)
    
    
    for j = (i+1):(length(ccROIidxGrow))
        if(~isempty(intersect(ccROIidxGrow{i}, ccROIidxGrow{j})))
            [idx, idy] = ind2sub([lenx, leny], cell2mat(ccROIidx([i;j])));
            boundingbox = [min(idx), max(idx), min(idy), max(idy)];
            bwLargeTmp1 = zeros(lenx, leny);
            bwLargeTmp1(ccROIidx{i}) = 1; 
            bwLargeTmp2 = zeros(lenx, leny);
            bwLargeTmp2(ccROIidx{j}) = 1;
            bwSmallTmp1 = bwLargeTmp1(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
            bwSmallTmp2 = bwLargeTmp2(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
            scoreMapTMP = scoreMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
            bwdistTwo = graydist(scoreMapTMP, [find(bwSmallTmp1(:) == 1)]) + graydist(scoreMapTMP, [find(bwSmallTmp2(:) == 1)]);
            bwdistTwo(isnan(bwdistTwo)) = inf;
            paths = double(bwdistTwo <= (min(bwdistTwo(:)) + 0.0001)); %this path is to find the 
            paths = bwmorph(paths, 'thin', inf);
            graphNode = [graphNode;[i,j,mean(scoreMapTMP(paths(:) == 1))]];
        end
        
    end
end
graphNode(isnan(graphNode(:,3)),:) = [];
G = graph(graphNode(:,1), graphNode(:,2), graphNode(:,3));

[T, pred] = minspantree(G,'Root',rootID);
ee = T.Edges.EndNodes;
% pathMap = zeros(lenx, leny);
pathMapCell = cell(size(ee,1),1);
boundingBoxCell = cell(size(ee,1),1);
parfor i = 1:size(ee,1)
%     disp(i)
    m = ee(i,1);
    n = ee(i,2);
    [idx, idy] = ind2sub([lenx, leny], cell2mat(ccROIidx([m;n])));
    boundingbox = [min(idx), max(idx), min(idy), max(idy)];
    bwLargeTmp1 = zeros(lenx, leny);
    bwLargeTmp1(ccROIidx{m}) = 1; 
    bwLargeTmp2 = zeros(lenx, leny);
    bwLargeTmp2(ccROIidx{n}) = 1;
    bwSmallTmp1 = bwLargeTmp1(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
    bwSmallTmp2 = bwLargeTmp2(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
    scoreMapTMP = scoreMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4));
    bwdistTwo = graydist(scoreMapTMP, [find(bwSmallTmp1(:) == 1)]) + graydist(scoreMapTMP, [find(bwSmallTmp2(:) == 1)]);
    bwdistTwo(isnan(bwdistTwo)) = inf;
    paths = double(bwdistTwo <= (min(bwdistTwo(:)) + 0.0001)); %this path is to find the 
    paths = bwmorph(paths, 'thin', inf);
    pathMapCell{i} = paths;
    boundingBoxCell{i} = boundingbox;
end
pathMap = zeros(lenx, leny);
largestRegionID = find(cellfun(@length,ccROIidx) == max(cellfun(@length,ccROIidx)));
newSynIdAllfinSkel = bwskel(logical(newSynIdAllfin));
newSynIdAllfinSkelROI = newSynIdAllfinSkel.*ccROI;
newSynIdAllfinSkelROIidx = label2idx(newSynIdAllfinSkelROI);
newSynIdAllfinDist = bwdist(1 - newSynIdAllfin);
widthAllRegion = cellfun(@(c) mean(newSynIdAllfinDist(c)), newSynIdAllfinSkelROIidx);
widthAllRegion = double(round(widthAllRegion));
for i = 1:size(ee,1)
    if(ee(i,1) == largestRegionID)
        se = strel('disk',widthAllRegion(ee(i,2)));
        boundingbox = boundingBoxCell{i};
        paths = imdilate(pathMapCell{i},se);
        pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) = pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) + paths;
    elseif(ee(i,2) == largestRegionID)
        se = strel('disk',widthAllRegion(ee(i,1)));
        boundingbox = boundingBoxCell{i};
        paths = imdilate(pathMapCell{i},se);
        pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) = pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) + paths;
    else
        se = strel('disk',round((widthAllRegion(ee(i,1)) + widthAllRegion(ee(i,2)))/2));
        boundingbox = boundingBoxCell{i};
        paths = imdilate(pathMapCell{i},se);
        pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) = pathMap(boundingbox(1):boundingbox(2), boundingbox(3):boundingbox(4)) + paths;
    end
end
imwrite(double( newSynIdAllfin + pathMap >0), fullfile(newfolder, 'astro_detection_ssp.tif'))

mask_ssp = double( newSynIdAllfin + pathMap >0);



end