function [sIdx, eIdx] = find_runs(mask)
    d = diff([false, mask, false]);
    sIdx = find(d==1);
    eIdx = find(d==-1)-1;
end
