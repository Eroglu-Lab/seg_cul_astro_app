function [refinedWidePart] = refineWideCompartment(widePart, bpregion)
% simply pick the one large connected component
    [lenx, leny] = size(widePart);
    refinedWidePart = zeros(lenx, leny);
    se = strel('disk',1);
    widePart = widePart + bpregion;
    widePart = imclose(widePart, se);
    
    widePartROI = bwlabel(widePart);
    widePartROIidx = label2idx(widePartROI);
    widePartROIidx = widePartROIidx(:);
    lenWidepart = cellfun(@length, widePartROIidx);
    refinedWidePart(widePartROIidx{lenWidepart == max(lenWidepart)}) = 1;
    
    
    
    
%     refinedWidePart(cell2mat())
%     se = strel('disk',1);
%     wholeMask = double( widePart + thinPart > 0);
%     widePartROI = bwlabel(widePart);
%     widePartROIidx = label2idx(widePartROI);
%     lengthx = cellfun(@length, widePartROIidx);
%     centerLabel = find(lengthx == max(lengthx));
%     wholeMaskROI = widePartROI + bwlabel(thinPart) + max(widePartROI(:)).*thinPart;
%     wholeMaskROIidx = label2idx(wholeMaskROI);
%     wholeMaskROIidx = wholeMaskROIidx(:);
% 
%     nodex = [];
%     wholeMaskROIidxGrowed = cell(length(wholeMaskROIidx), 1);
%     for i = 1:length(wholeMaskROIidxGrowed)
%         tmpID = wholeMaskROIidx{i};
%         tmpIDGrowed = regionGrowxx(tmpID, 1, lenx, leny, xxshift, yyshift);
%         wholeMaskROIidxGrowed{i} = tmpIDGrowed;
%     end
%     for i = 1:length(wholeMaskROIidx)
%         tmpIDGrowed = wholeMaskROIidxGrowed{i};
%         tmpIDlabel = unique(wholeMaskROI(tmpIDGrowed));
%         tmpIDlabel(tmpIDlabel == 0) = [];
%         tmpIDlabel(tmpIDlabel == i) = [];
%         for k = 1:length(tmpIDlabel)
%             nodex = [nodex; [i, tmpIDlabel(k), 1./mean(maskDist(wholeMaskROIidx{i}))]];
%         end
%     end
%     if(~isempty(nodex))
%         G0 = graph(nodex(:,1), nodex(:,2), nodex(:,3));
%         for i =1:max(widePartROI(:))
%             if(i ~= centerLabel)
%                 pathx = shortestpath(G0,centerLabel, i);
%                 refinedWidePart(cell2mat(wholeMaskROIidx(pathx))) = 1;
%             end
%         end
%     else
%         refinedWidePart = widePart;
%     end
end