function [perfIdx] = perfusionIdx(x,x_dtr,len,pad)
%%% Perfusion index for PPG
%%% Assume x and x_dtr are 1-column signals
    segs = buffer(x,len,(len-pad))';
    segs = segs(len/pad:end,:);
    
    segs_dtr = buffer(x_dtr,len,(len-pad))';
    segs_dtr = segs_dtr(len/pad:end,:);
    
    perfIdx = (max(segs_dtr')' - min(segs_dtr')')./(mean(segs,2)+0.0001)*100;
    
end