function [features] = trouble_handling(signal)

sig = signal(:,ep)-mean(signal(:,ep));
% First get the fundamental phase
% basic parameters for STFT
sampling_rate = 100; %(need to change!)
basicTF.win = sampling_rate*10+1; %window length(sec)(need to change!)
basicTF.hop = 1; %(need to change!)([hop] sample)
basicTF.fs = sampling_rate;
basicTF.fr = 0.02; %frequency resolution %(need to change!)
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection;
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 3.15/basicTF.fs; % highest frequency/sampling freq(need to change!)
advTF.LowFreq = 0.01/basicTF.fs; % lowest frequency/sampling freq(need to change!)
advTF.lpc = 0;
[~, ~, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(sig, advTF.LowFreq-0.02/basicTF.fs, ...
    advTF.HighFreq, basicTF.fr/basicTF.fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0); 

% Curve extraction & reconstruction
[h, ~, ~] = hermf(100*10+1,1,5); h0 = h(floor(size(h,2)/2)+1);
% fundamental frequency
idx0 = find(tfrsqtic*sampling_rate>0.5 & tfrsqtic*sampling_rate<1.9);
[fund] = CurveExt_M(abs(tfrsq(idx0,:))', 5.0);
fund = fund + idx0(1) - 1;
% plot((0:length(fund)-1)/sampling_rate, sampling_rate*tfrsqtic(fund), 'r');
tmp = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 100, fund, 0.2, h0);
phi_fund = unwrap(angle(tmp))/2/pi; % This the phase
%%
% Decomposite the new function
val = linspace(phi_fund(1),phi_fund(end),length(phi_fund));
psiInv = zeros(size(val));
for k = 1:length(psiInv) % sampling in codomain
    [~, psiInv(k)] = min(abs(phi_fund - val(k)));
end
% figure; plot(psiInv);
% figure; plot((0:length(phi_fund)-1)./100,phi_fund);
sig_prime = sig(psiInv);
% figure; plot((0:length(sig_prime)-1)./100, real(sig_prime));
%%
% basic parameters for STFT
sampling_rate = 100; %(need to change!)
basicTF.win = sampling_rate*10+1; %window length(sec)(need to change!)
basicTF.hop = 1; %(need to change!)([hop] sample)
basicTF.fs = sampling_rate;
basicTF.fr = 0.02; %frequency resolution %(need to change!)
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection;
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 8.15/basicTF.fs; % highest frequency/sampling freq(need to change!)
advTF.LowFreq = 0.01/basicTF.fs; % lowest frequency/sampling freq(need to change!)
advTF.lpc = 0;
[~, ~, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(sig_prime, advTF.LowFreq-0.02/basicTF.fs, ...
    advTF.HighFreq, basicTF.fr/basicTF.fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);
[h, ~, ~] = hermf(100*10+1,1,5); h1 = h(floor(size(h,2)/2)+1);

% figure;
% set(gcf,'Position',[100 50 1000 700]);
% imageSQ((0:size(tfrsq,2)-1)/100, tfrsqtic*100, abs(tfrsq), 0.99);
% axis xy; colormap(1-gray); colorbar
% title('TF representation win14');
% xlabel('time(sec)','FontSize',20);
% ylabel('frequency(Hz)','FontSize',20);
% ax = gca;
% ax.FontSize = 20;
% hold on;

% Curve extraction
lambda = 5.0; gamma = 9000.0;
harmonics = zeros(6,size(tfrsq,2));
% fundamental frequency
idx0 = find(tfrsqtic*sampling_rate>0.5 & tfrsqtic*sampling_rate<1.9);
[f] = CurveExt_M(abs(tfrsq(idx0,:))', lambda);
f = f + idx0(1) - 1;
% plot((0:length(f)-1)/sampling_rate, sampling_rate*tfrsqtic(f), 'r');
harmonics(1,:) = f;

% the 2nd harmonic
frng_low = 0.5; frng_high = 0.5;
rng_low = ceil(frng_low*size(tfrsq,1)/8.14);
rng_high = ceil(frng_high*size(tfrsq,1)/8.14);
tfrsq1 = tfrsq;
for t = 1:length(f)
    tfrsq1(f(t)+floor(prctile(f,90))+rng_high+1:end, t) = 0;
    tfrsq1(1:f(t)+floor(prctile(f,10))-rng_low, t) = 0;
end
idx = find(tfrsqtic*sampling_rate >= min(2*tfrsqtic(harmonics(1,:))*sampling_rate)-frng_low);
[c] = CurveExt_M_v2(abs(tfrsq1(idx,:))', lambda, 2, idx(1), tfrsqtic, f, gamma, tfrsqtic);
c = c + idx(1) - 1;
% plot((0:length(c)-1)/100, 100*tfrsqtic(c), 'r');
harmonics(2,:) = c;

for k = 3:6
    frng = 0.5;
    tfrsq1 = tfrsq;
    rng = ceil(frng*size(tfrsq1,1)/8.14);
    for t = 1:length(fund)
        tfrsq1(c(t)+floor(prctile(f,75))+rng+1:end, t) = 0;
        tfrsq1(1:c(t)+floor(prctile(f,25))-rng, t) = 0;
    end
    idx = find(tfrsqtic*sampling_rate >= min(tfrsqtic(f)+tfrsqtic(c))*sampling_rate-frng);
    [c] = CurveExt_M_v2(abs(tfrsq1(idx,:))', lambda, k/(k-1), idx(1), tfrsqtic, c, gamma, tfrsqtic);
    c = c + idx(1) - 1;
%     plot((0:length(c)-1)/sampling_rate, sampling_rate*tfrsqtic(c), 'r');
    harmonics(k,:) = c;
end
% hold off;
%% Wrapping back
% phi_fund ~ val(phi_val)
phi_val = zeros(size(phi_fund));
for k = 1:length(phi_val)
    [~,phi_val(k)] = min(abs(val-phi_fund(k)));
end

% disp(num2str(norm(val(phi_val)-phi_fund)/norm(phi_fund)));
% figure;
% plot((0:length(sig)-1)/100,sig); hold on;
% plot((0:length(sig)-1)/100,sig_prime(phi_val)); hold off;
% disp(num2str(norm(sig-sig_prime(phi_val))/norm(sig)));
%%
recon = zeros(size(harmonics)); super = zeros(size(sig));
for h = 1:6
    tmp = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 100, harmonics(h,:), 0.2, h1);
    recon(h,:) = tmp(phi_val);
    if h == 1
        super = recon(1,:);
    else
        super = super + recon(h,:);
    end
end
% figure;
% plot(sig); hold on; plot(real(super)); hold off;
% disp(num2str((norm(sig-real(super)'))/norm(sig)));
% figure;

features = zeros(1,15*5);
for h = 2:6
    tmp = (recon(1,:).^h).*conj(recon(h,:));
    tmp = unwrap(angle(tmp))/2/pi;
    tmp = tmp(1,3*100+1:end-3*100);
    % plot(unwrap(angle(tmp))/2/pi);
    
    % features
    tmpMoment = zeros(1,11);
    for mm = 1:size(tmpMoment,2)
        tmpMoment(1,mm) = moment(tmp, mm+2)/std(tmp,1)^(mm+2);
    end
    features(1,(h-2)*15+1:(h-1)*15) = [
        prctile(tmp,25) median(tmp) prctile(tmp,75) tmpMoment mad(tmp)
    ];
end

end