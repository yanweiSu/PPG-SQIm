function [ssqi] = sSQI(x,len,pad)
%%% Skewness SQI for PPG
%%%Assume x is row data
    segs = buffer(x,len,(len-pad))';
    segs = segs(len/pad:end,:);
    segskew = skewness(segs,0,2); % signal/de-sample bias/dim
    ssqi = segskew;
%     ssqi = repelem(segskew',w);
%     ssqi = ssqi(1:length(x));
end