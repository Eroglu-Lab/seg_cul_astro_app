function [refinedWidePart] = refineLayer1Compartment(layer1, layer2, layer34, thinBranchskel,thinBranchDist,dist2center, xxshift, yyshift)
% wide part should be all be connected though due to the low intensity the
% thin part might not always be connected
% several rules for the wide part
% basically, the small/scattered ones are not wide part
% two nearby wide part should be connected
    [lenx, leny] = size(layer2);
    refinedWidePart = zeros(lenx, leny);
    se = strel('disk',1);
    finePart = layer2 + layer34;
    layer2 = imclose(layer2, se);
    layer2ROI = bwlabel(layer2);
    layer2ID = [1:max(layer2ROI(:))];
    layer34ROI = bwlabel(layer34);
    layer34ID = [1:max(layer34ROI(:))];
    allMap = layer1 + (layer2ROI + double(layer2ROI > 0)) + (layer34ROI + double(layer34ROI > 0).*(1 + max(layer2ID)));
    layer2ID_new = 1 + layer2ID;
    layer34ID_new = 1 + max(layer2ID) + layer34ID;
    % build the graph of the branches (search the path from the root branches to each ) 
    maskBP = bwmorph(thinBranchskel,'branchpoints');
    [branch_s, maskBP] = findRealBranch(thinBranchskel,maskBP);
    branch_s_ROI = bwlabel(branch_s);
    branch_s_ROI_idx = label2idx(branch_s_ROI);
    maskBP2 = imdilate(maskBP, se);
    branch_s2 = imdilate(branch_s, se);
    branch_s2_ROI = bwlabel(branch_s2);
    maskBP2_ROI = bwlabel(maskBP2);
    maskBP2_ROI_idx = label2idx(maskBP2_ROI);
    layer1_ss = imdilate(layer1, se);
    source_bc = unique(branch_s2_ROI(layer1_ss(:) == 1));
    source_bc(source_bc == 0) = [];
    graphx = [];
    for i = 1:length(maskBP2_ROI_idx)
        bcLabel = unique(branch_s2_ROI(maskBP2_ROI_idx{i}));
        bcLabel(bcLabel ==0 ) = [];
        graphx = [graphx; [bcLabel(1:end-1), bcLabel(2:end), ones(length(bcLabel(1:end-1)),1)]];
    end
    source_n = max(graphx(:)) + 1;
    graphx = [graphx; [repmat(source_n, length(source_bc),1),source_bc(:), zeros(length(source_bc), 1)]];
    graph_s = graph(graphx(:,1), graphx(:,2), graphx(:,3));
    
    % extract the branch labels that are in the layer2 regions, for each
    % layer2 connected components choose the small branch that is close to
    % the 
    layer2ROI_idx = label2idx(layer2ROI);
    ccbranchlabel = zeros(length(layer2ROI_idx),1); % the terminal branch
    allBranchLabelLayer2 = cell(length(layer2ROI_idx),1); % all the branches id inside each of the layer2 region
    for i = 1:length(layer2ROI_idx)
        uid = unique(branch_s_ROI(layer2ROI_idx{i}));
        uid(uid == 0) = [];
        if(~isempty(uid))
            allBranchLabelLayer2{i} = uid;
            dist2centerBC = cellfun(@(c) min(dist2center(c)), branch_s_ROI_idx(uid));
            ccbranchlabel(i) = uid(dist2centerBC == min(dist2centerBC));
        end
    end
    layer2BranchMask = zeros(length(branch_s_ROI_idx),1);
    branch_s_ROI_idx = branch_s_ROI_idx(:);
    for i = 1:length(layer2ROI_idx)
        curID = ccbranchlabel(i);
        if(curID ~= 0)
            ssPath = shortestpath(graph_s, source_n, curID);
            if(~isempty(ssPath))
                ssPathx = ssPath;
                ssPathx(1) = [];
                ssPathx(end) = [];
                % if the path goes through other connected components of layer2
                % then remove it
                
                if(~isempty(ssPathx))
                    pathid = cell2mat(branch_s_ROI_idx(ssPathx));
                    ccLayer2 = unique(layer2ROI(pathid));
                    ccLayer2(ccLayer2 == 0) = [];
                    ccLayer2(ccLayer2 == i) = [];
                    if(~isempty(ccLayer2))
                        layer2BranchMask(allBranchLabelLayer2{i}) = 1;
                        layer2BranchMask(ssPathx) = 1;
                    end
                elseif(isempty(ssPathx))
                    layer2BranchMask(allBranchLabelLayer2{i}) = 1;

                end
            end
        end
    end
    maskDist2 = max(thinBranchDist(:)) - thinBranchDist;
    maskDist2 = imimposemin(maskDist2, (branch_s));
    maskWS = double(watershed(maskDist2)).*finePart;
    maskWSidx = label2idx(maskWS);
    maskWSidx = maskWSidx(:);
    maskWideLabel = unique(maskWS(cell2mat(branch_s_ROI_idx(logical(layer2BranchMask)))));
    maskWideLabel(maskWideLabel == 0) = [];    
    refinedWidePart(cell2mat(maskWSidx(maskWideLabel))) = 1;
    refinedWidePart
end