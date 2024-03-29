function [recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, c, Band, coeff) 

% Band = absolute frequency band usually is chosen to be 0.1Hz
% tfrsqtic = in the scale of 1Hz %(Hao) No it's actually in the scale
% betwenn 0~0.5 as the ratio of sampling freuqency
% c is the extracted curve.

%coeff=0.1152; % computed in the numerical example 'genFig1Fig2_test_coeff_SST.m'
%(Hao) Not quite sure what this is

alpha = tfrsqtic(2)-tfrsqtic(1);
RR = round((Band/Hz)/alpha);


recon = [] ;

C = 2 * alpha / coeff ;

for kk = 1: length(c)
	idx = max(1,c(kk)-RR): min(length(tfrsqtic),c(kk)+RR); % Extract relative magnitude of each frequency
	recon(kk) = C * sum(tfrsq(idx,kk),1) ;
end

