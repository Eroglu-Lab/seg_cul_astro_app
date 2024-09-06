
path = 'F:\eroglu_lab\seg_cul_astrocyte_app\example\0025_astrocyte\';
file3 = 'astro_detection_ssp.tif';

iters1 = 1; %11x11x5
xxshift1 = zeros(2*iters1+1, 2*iters1+1);
yyshift1 = zeros(2*iters1+1, 2*iters1+1);
for i = -iters1:iters1
    for j = -iters1:iters1
        xxshift1(i+iters1+1,j+iters1+1) = i;
        yyshift1(i+iters1+1,j+iters1+1) = j;
    end
end

se = strel('disk', 3);
mask_astro2 = imopen(mask_astro, se);
masksspSkel2 = bwskel(logical(mask_astro2),'MinBranchLength',30);
figure; imagesc(mask_astro2 + masksspSkel2)


mask_astro = im2double(tiffreadVolume(fullfile(path,file3)));
masksspSkel = bwskel(logical(mask_astro),'MinBranchLength',15);
maskDist = bwdist(1 - mask_astro);
[widePart, thinPart, terminalBranchOut,widebranchPointRegion, widebranch_s_idx2, widthMean,widemaskWSidx2] = genCompartmentWideThin(mask_astro, masksspSkel,maskDist, 6,xxshift1, yyshift1);

[refinedWidePart] = refineWideCompartment(widePart, widebranchPointRegion);
refinedWidePart = refinedWidePart.*mask_astro;

thinBranch = mask_astro - refinedWidePart;
thinBranchskel = bwskel(logical(thinBranch),'MinBranchLength',15);
thinBranchDist = bwdist(1 - thinBranch);
[layer2, layer34] = genCompartmentWideThin_update(thinBranch, thinBranchskel,thinBranchDist, 2);
[refinedlayer2] = refineLayer1Compartment2(refinedWidePart, layer2, layer34);
layer34 = thinBranch - refinedlayer2;