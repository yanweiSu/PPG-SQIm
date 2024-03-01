function [sig, tfr, tfrtic, phi_value,PsiInv] = iterWarping(sig, basicTF, phi, iter)
% This program iteratively warps on signal's fundamental, 1st harmonics, 2nd
% harmonics, ... and so on. The number of iterations is set by [I].
% [sig] is the signal
% [phi] is the pre-extracted fundamental phase
% [basicTF] is as following:
fs = basicTF.fs;
win = basicTF.win;
fr = basicTF.fr;
HighFreq = 8/fs;
LowFreq = basicTF.LowFreq/fs;
hop = 1;
% Amplitude (for reconstruction)
[h, ~, ~] = hermf(win,1,5);
h0 = h(floor(size(h,2)/2)+1);

I = length(iter);
phi_value = cell(I,1);
PsiInv = cell(I,1);
for i = 1:I
    fprintf("warping %d(%d)\n", i, iter(i));
    %% construct "inverse of phi(psiInv)"
    M = ceil(phi(end)-phi(1))*fs;  % sampling on the fundamental phase
    % tau = 1/fs;
    val = linspace(phi(1), phi(end), M);
    psiInv = zeros(1, floor(M/iter(i)));
    for k = 1:length(psiInv)
        [~, psiInv(k)] = min(abs(phi-val(iter(i)*k)));
    end
    PsiInv{i} = psiInv;
    %% Prepare for unwarping back (i.e. phi_hat in the note)
    phi_val = zeros(size(phi));
    val2 = val(iter(i):iter(i):iter(i)*floor(M/iter(i)));
    for k = 1:length(phi_val)
        [~, phi_val(k)] = min(abs(val2-phi(k))); % phi's value index
    end
    phi_value{i} = phi_val;

    %% Warping and SST on the warped signal
    sig = sig(psiInv);
    sig = sig - mean(sig);
    [~, ~, tfr, ~, tfrtic] = ConceFT_sqSTFT_C(sig, ...
        LowFreq, HighFreq, fr/fs, hop, win, 1, 5, 1, 1, 0);
    % TFR plot (For checking)
%     if i <= I
%         figure;
%         set(gcf,'Position',[100 50 1000 700]);
%         imageSQ((0:size(tfr,2)-1)./fs, tfrtic*fs, abs(tfr), 0.99);
%         axis xy; colormap(1-gray); colorbar
%         xlabel('time(sec)','FontSize',20);
%         ylabel('frequency(Hz)','FontSize',20);
%         ax = gca;
%         ax.FontSize = 20;
%     end

    %% extract the fundamental phase phi of the warped signal
    if i < I
    % The range of fundamental IF. Need to modify
    idx0 = find(tfrtic*fs>(iter(i+1)-0.4) & tfrtic*fs<(iter(i+1)+0.4));
    [fund] = CurveExt(abs(tfr(idx0,:))', 3.0);
    fund = fund + idx0(1) - 1;
    tmp = Recon_sqSTFT_v2(tfr, tfrtic, fs/hop, fund, 0.2, h0);
    phi = unwrap(angle(tmp))/2/pi;  % This is the phase
    end
end