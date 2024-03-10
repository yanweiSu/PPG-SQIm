function sqi = SQI_M(x, fs)

% Make the signal to be a column vector
x = reshape(x, [], 1);

addpath('./lib')
addpath('./lib/Algorithm') % The code for extending signal
addpath('./lib/TF_analysis')

%% TF parameters for PPG
% Frequency bins for SST
fr = 0.02; fr_first = 0.025;    % for RD

% Time resolution for SST
hop = 5;    % samples

% Frequency range for SST
HighFreq = 18.5; LowFreq = 0.5; % Hz

% Window function:
window1 = 10;   % seconds
supp = 5;

% Amplitude (for reconstruction)
[h, ~, ~] = hermf(fs*window1+1, 1, supp);
h0 = h(floor(size(h,2)/2)+1);

%% Pre-processing: bandpass and signal extension
HOP = 1;
extSEC = 5 ; % the extension is of extSEC second
L = round( extSEC*fs );
extM = round(1.5*L);    % dimension of embedding / signals length
extK = round(2.5*extM); % number of points to estimate A / size of datasets
method.name = 'edmd';
method.param = 100;

% ----------- Bandpass
[b, a] = butter(6, [0.5, 20]/(fs/2)); % bandpass filter setting.
sig_dtr = filtfilt(b, a, x);% bandpass filter on x

% Reduce the signal by extending extSEC on both side
xxsig = SigExtension(sig_dtr, fs, HOP, extK, extM, extSEC, method);


%% ----------- multi-RD to get the fundmental IF ------------
% ---------- Faster RD: curve 1 & 2 [c_tmp] ------------
% TFR with hop = [fs] samples
[~, ~, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(xxsig, LowFreq/fs, HighFreq/fs, ...
    fr_first/fs, fs, window1*fs+1, 1, supp, 1, 1, 0);
% Truncate the extended part
tfrsq = tfrsq(:, extSEC+1:end-extSEC);

% multi-RD: 2-curves version
ver = 2;
nth_harmonics = 2;
tt = [0 size(tfrsq,2)];
% lambda = [0.1 0.1];
% mu = .005;
lambda = [1/fs 1/fs];
mu = .07;
fundRng = [.7 3];
[c_tmp] = MultiCurveExt(tfrsq, tfrsqtic*fs, lambda, mu, ...
    nth_harmonics, fundRng, 1.0, ver, tt, [-.5 .5]);

% TFR with hop = 5 samples
[~, ~, tfrsq5, ~, tfrsqtic5] = ConceFT_sqSTFT_C(xxsig, LowFreq/fs, HighFreq/fs, ...
    fr/fs, hop, window1*fs+1, 1, supp, 1, 1, 0);
tfrsq5 = tfrsq5(:, extSEC*fs/hop+1:end-extSEC*fs/hop);

% Interplolate the IF curve to be of fs/5 sampling rate
IF_fs1hz = tfrsqtic(c_tmp)*fs;
IF_fs1hz_interp = interp1(1:(fs/hop):size(tfrsq5,2), IF_fs1hz, ...
    1:size(tfrsq5,2), "linear");

%% Harmonics modes with fs/5 Hz
cALL = [];
tfr_curveExt = zeros(size(tfrsq5));
alpha = tfrsqtic5(2)-tfrsqtic5(1);
harmodes = zeros(6, size(tfrsq5,2));
for l = 1:6
    if l <= 2
        upperHz = IF_fs1hz_interp(:,l)*1.2;
        lowerHz = IF_fs1hz_interp(:,l)*0.8;
    else
        upperHz = (l+0.3)*tfrsqtic5(cALL(:,1))*fs;
        lowerHz = (l-0.3)*tfrsqtic5(cALL(:,1))*fs;
    end

    for m = 1:size(tfrsq5,2)
        % Lidx = max( 1 , round( (reflower(m)-LowFreq)/fr ) );
        % Hidx = min( length(tfrsqtic5) , round( (refupper(m)-LowFreq)/fr ) );
        Lidx = max( 1 , ceil(1+(lowerHz/fs-tfrsqtic5(1))/alpha) );
        Hidx = min( length(tfrsqtic5) , floor(1+(upperHz/fs-tfrsqtic5(1))/alpha) );
        tfr_curveExt(Lidx:Hidx, m) = tfrsq5(Lidx:Hidx, m);
    end
    cALL = [cALL CurveExt(abs(tfr_curveExt).', 1.0)];
    tfr_curveExt = zeros(size(tfrsq5));
    harmodes(l,:) = real(Recon_sqSTFT_v2(tfrsq5, tfrsqtic5, fs, cALL(:,l), .2, h0));
end

% The reconstructed ANHM signal
recon = sum(harmodes);

%% dws=1 modes
% [~, ~, tfrsq1, ~, tfrsqtic1] = ConceFT_sqSTFT_C(xxsig, ...
%     basicTF.LowFreq/fs, basicTF.HighFreq/fs, ...
%     fr/fs, 1, basicTF.win, 1, 5, 1, 1, 0);
% % Truncate the TF of extended part
% tfrsq1 = tfrsq1(:, extSEC*fs+1:end-extSEC*fs);  % fs = fs
% 
% % c is a dws=5 curve IF curve
% c = cALL(:,1);
% c = repmat(c, 1, hop);
% c = reshape(c', 1, [])';
% 
% tfr_curveExt = tfrsq1;
% ref = tfrsqtic(c)*fs;
% for m = 1:size(tfrsq1,2)
%     % ----- trim the TFR for CurveExt -----
%     lowbdd = max(1, round( (ref(m)*0.75-basicTF.LowFreq)/fr ));
%     upbdd = min(size(tfrsq1,1), round( (ref(m)*1.25-basicTF.LowFreq)/fr ));
%     tfr_curveExt([1:lowbdd upbdd:end], m) = 0;
% end
% cALL1 = CurveExt(abs(tfr_curveExt)', 1.0);

%% ----- SQI from hop=1 + warped modes ------
% % warping
% phi1_est = unwrap(angle(fund_est))/2/pi;
% [~, tfr1_warped, tfrtic1_warped, phi_value, ~] = ...
%     iterWarping(xxxsig, basicTF, phi1_est, 1);
% harms = [];
% for i = 1:6
%     % IF
%     idx0 = find(tfrtic1_warped*fs>(i-0.35) & tfrtic1_warped*fs<(i+0.35));
%     [c] = CurveExt(abs(tfr1_warped(idx0,:))', 1.0);
%     c = c + idx0(1) - 1;
% 
%     % harmonics
%     est = Recon_sqSTFT_v2(tfr1_warped, tfrtic1_warped, fs, c, .2, h0);
%     est = est(phi_value{1});
%     est = resample(est, 1, 5);  
%     harms = [harms; est];
% end
% recon_warped = sum(harms);
% 
% % ------- SQI
% PVP_dtr = sig_dtr;
% fs = fs;
% % resample the reconstructed cardiac signal into the same sampling rate as
% % the PVP signal.
% recon_cardiac = interp1( (1:hop:length(PVP_dtr))/fs, ...
%     real(recon_warped), (1:length(PVP_dtr))/fs, 'pchip')';
% reinterp_sig = interp1( (1:hop:length(PVP_dtr))/fs, ...
%     real(harms).', (1:length(PVP_dtr))/fs, 'pchip')';
% reinterp_sig = reshape(reinterp_sig, 6, []);
% 
% PVP_noise = PVP_dtr - recon_cardiac;
% SQI_tmp = zeros(length(recon_warped),1);
% for m = 1:length(recon_warped)-1    % -1???
%     SQI_tmp(m) = norm(harms(:,m),2)^2 / ( norm(harms(:,m),2)^2 + ...
%         sum(PVP_noise( (m-1)*hop+1 : m*hop ).^2)/hop);
% end
% 
% epoch_num = floor((length(PVP_dtr)-fs*epoch_len)./(shift_n*hop));
% % epoch_num = floor( (length(sig_dtr)-5*fs) / 30 ) = 25*fs / 30
% % shift_n = 30/5
% 
% SQI = zeros(epoch_num,1);
% avgAHI = zeros(epoch_num,6);
% 
% % initial
% s = 0;
% t = fs/hop*epoch_len;
% % t = fs
% 
% for i = 1:epoch_num
%     SQI(i)= median(SQI_tmp(s+1:t));
%     base = max(norm(PVP_dtr(s*hop+1:t*hop)), 0.001); % previent nan and 0
%     avgAHI(i,:) = max(sqrt(sum(reinterp_sig(:,s*hop+1:t*hop).^2, 2))/base, -100);
%     s = s + shift_n;    % s + 6 samples
%     t = t + shift_n;    % fs*step_len;
% end
% tmp_SQI = SQI;
% avgAHIWarped_list{k} = avgAHI;
% SQIWarped_table(k,:) = tmp_SQI;

%% SQI_M
PVP_dtr = sig_dtr;
% Reconstruct the noise with the original sampling rate
recon_cardiac = interp1((1:hop:length(PVP_dtr))/fs, recon, ...
    (1:length(PVP_dtr))/fs, 'pchip')';
PVP_noise = PVP_dtr - recon_cardiac;

% SQI_M definition
sqi = zeros(size(recon));
for i = 1:length(recon)
    sqi(i) = norm(harmodes(:,i),2)^2 / ...
        (norm(harmodes(:,i),2)^2 + sum(PVP_noise((i-1)*hop+1:i*hop).^2)/hop);
end


end



% %% For ML
% shift_n = round(0.5*(fs/hop));
% % PPG epoch is set to be 5 seconds for SQI
% epoch_len = 5;
% epoch_num = floor((length(PVP_dtr)-fs*epoch_len)./(shift_n*hop));
% SQI = zeros(epoch_num,1);
% 
% % initial
% s = 0;
% t = fs/hop*epoch_len;
% for i = 1:epoch_num
%     SQI(i)= median(SQI_tmp(s+1:t));
%     s = s + shift_n;    %fs/hop*step_len;
%     t = t + shift_n;    %fs*step_len;
% end
% 
% tmp_SQI = SQI.';    % Into a row vector
% 
% %% Parameters
% % segment length for sSQI, eSQI, and perfIdx
% pad = 0.5*fs;
% % PPG_M shift.
% PVP_SQI_step = round((pad/fs)*(fs/hop));
% 
% %%%% !!!NEED ATTENTION IF YOU CHANGE WINDOW LENGTH or PAD!!!
% t0 = (1:30*fs)./fs;    % 0s to 30s, fs=64Hz
% % The sampling period of PPG-SQI and AHI are PVP_SQI_step*(hop/fs)
% t2 = 2.5 + (0:52)*PVP_SQI_step*(hop/fs);
% 
% %%
% tmp_SQI = interp1(t2, tmp_SQI, t0, 'pchip');
% tmp_SQI(1:(t2(1)*fs)-1) = tmp_SQI(t2(1)*fs);
% tmp_SQI(t2(end)*fs+1:end) = tmp_SQI(t2(end)*fs);
% 
% % Temp: Downsample to 2Hz
% downsampling_rate = 2;
% tmp_SQI = tmp_SQI(1:fs/downsampling_rate:end);
% 
% % Truncation to 5~25s
% trim_time = 5;
% tmp_SQI = tmp_SQI(downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);
% sqi = tmp_SQI';
% 
