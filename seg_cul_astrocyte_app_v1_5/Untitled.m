mask = double(imread('F:\eroglu_lab\pia\wrong_example\astro_detection_ssp.tif'));
mask(mask > 0) = 1;
maskskel =  bwskel(logical(mask),'MinBranchLength',15);
maskDist = bwdist(1 - logical(mask));
iters1 = 1;
xxshift= zeros(2*iters1+1, 2*iters1+1);
yyshift = zeros(2*iters1+1, 2*iters1+1);
for i = -iters1:iters1
    for j = -iters1:iters1
        xxshift(i+iters1+1,j+iters1+1) = i;
        yyshift(i+iters1+1,j+iters1+1) = j;
    end
end
thres0 = 7
[widePart, thinPart, terminalBranchOut, dist2center] = genCompartmentWideThin(mask, maskskel,maskDist, thres0,xxshift, yyshift);
%instead of connecting the widePart, we can connect the thin part in the
%new design
[refinedWidePart] = refineWideCompartment(maskDist, dist2center, widePart, thinPart, xxshift, yyshift);