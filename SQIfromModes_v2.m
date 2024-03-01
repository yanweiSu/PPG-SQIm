clear;
close all;

addpath('./lib')
addpath('./lib/Algorithm') % The code for extending signal
addpath('./lib/TF_analysis')

PLOT = 0;

db = "TROIKA";
% db = "WESAD";
% db = "DaLiA_test";
% db = "DaLiA_train";
sig = csvread(strcat("./DATASETS/", db, "/", db, "_seg.csv")).';
n_seg = size(sig,2);

%% TF settings for PPG
fs = 64;
fr = 0.02;
fr_first = 0.025;
hop = 5;
HighFreq = 18.5;
LowFreq = 0.5;
window1 = 10;

basicTF.hop = hop;
basicTF.fs = fs;
basicTF.fr = fr; % frequency resolution%(need to change!)
basicTF.win = fs*window1+1;
basicTF.HighFreq = 18.5;
basicTF.LowFreq = 0.5;
f_component_size = (basicTF.HighFreq-basicTF.LowFreq)/basicTF.fr;
% Amplitude (for reconstruction)
[h, ~, ~] = hermf(basicTF.win,1,5);
h0 = h(floor(size(h,2)/2)+1);

% Extending signal
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets
method.name = 'edmd';
method.param = 100 ;

epoch_len = 5;
shift_n = round(0.5*(fs/hop));

%% For PPG
%%%%% indices saving
epoch_num = floor((length(sig(:,1))-fs*epoch_len)./30);

SQI_table = zeros([n_seg,epoch_num]);
SQIWarped_table = zeros([n_seg,epoch_num]);
avgAHI_list = cell([n_seg,1]);
avgAHIWarped_list = cell([n_seg,1]);
modeNum_list = zeros(n_seg,1);
ref_SQI_table = zeros([n_seg, epoch_num]);


%% Added at 2023/12/31: for T_examples 50, 108
start = 1;
allModes = load(strcat("./SQI_results_v2/modes_", db, ".mat"));
allModes_dws5 = allModes.modes_dws5_list;
allModes_fund_dws1 = allModes.modes_fund_list;
allIF = allModes.IF_list;

for k = 1:n_seg
% for k = [50 108]
    tic
    fprintf("seg %d/%d\n", k, n_seg);
    % ----------- Bandpass
    [b, a] = butter(6, [0.5, 20] / (fs/2));    % bandpass filter setting
    sig_dtr = filtfilt(b, a, sig(:,k));    % bandpass filter

    % Reduce the signal by extending extSEC on both side
    xxsig = SigExtension(sig_dtr, fs, HOP, extK, extM, extSEC, method);
    xxxsig = xxsig(extSEC*fs+1 : end-extSEC*fs);

    %% load the reconstructed harmonics
    % 6 harmoncis, dws = 5
    rec_sig = allModes_dws5{k};
    recon = sum(rec_sig);

    % 1 fundamental, dws = 1
    harms = allModes_fund_dws1{k};

    % IF, dws = 1
    IFcurves_fsHz = allIF{k}.dws1;

    %% --------- estimate number of modes ---------
    fund_est = harms(1,:);
    phi_est = unwrap(angle(fund_est));
    A_est = abs(fund_est);
    % figure; plot(real(fund_est));

    R_max = floor( (fs/2) / max(IFcurves_fsHz) );
    selection = zeros(R_max, 1);
    for rr = 1:R_max
        N = length(fund_est);
        C = cosenos(A_est, phi_est, rr);
        coef_est = ((sig_dtr.')*C')/(C*C');
        x_R = coef_est*C;   % this is a simple reconstruction
        
        mse = norm(sig_dtr.'-x_R)^2/N;
        selection(rr) = N*N*mse/((N-2*rr-1)^2);   % G
    end      
    [~,R] = min(selection(:,1));    % use GCV criterion for now
    R = min([R 6]);   % maximum is 6
    modeNum_list(k) = R;


    %%
    %%%%%% This is for "parfor"
    basicTF = [];
    basicTF.hop = 5;
    basicTF.fr = fr;
    basicTF.fs = fs;
    basicTF.win = fs*window1+1;
    basicTF.LowFreq = 0.5;
    basicTF.HighFreq = round( (max(IFcurves_fsHz) + .5)*7 );
    [h, ~, ~] = hermf(basicTF.win,1,5);
    h0 = h(floor(size(h,2)/2)+1);
    %%%%%%

    %% ----- SQI from hop=1 + warped modes ------
    % warping
    phi1_est = unwrap(angle(fund_est))/2/pi;
    [~, tfr1_warped, tfrtic1_warped, phi_value, ~] = ...
        iterWarping(xxxsig, basicTF, phi1_est, 1);
    harms = [];
    for i = 1:6
        % IF
        idx0 = find(tfrtic1_warped*fs>(i-0.35) & tfrtic1_warped*fs<(i+0.35));
        [c] = CurveExt(abs(tfr1_warped(idx0,:))', 1.0);
        c = c + idx0(1) - 1;

        % harmonics
        est = Recon_sqSTFT_v2(tfr1_warped, tfrtic1_warped, fs, c, .2, h0);
        est = est(phi_value{1});
        est = resample(est, 1, 5);  
        harms = [harms; est];
    end
    recon_warped = sum(harms);

    % ------- SQI
    PVP_dtr = sig_dtr;
    sampling_rate = fs;
    % resample the reconstructed cardiac signal into the same sampling rate as
    % the PVP signal.
    recon_cardiac = interp1( (1:basicTF.hop:length(PVP_dtr))/sampling_rate, ...
        real(recon_warped), (1:length(PVP_dtr))/sampling_rate, 'pchip')';
    reinterp_sig = interp1( (1:basicTF.hop:length(PVP_dtr))/sampling_rate, ...
        real(harms).', (1:length(PVP_dtr))/sampling_rate, 'pchip')';
    reinterp_sig = reshape(reinterp_sig, 6, []);

    PVP_noise = PVP_dtr - recon_cardiac;
    SQI_tmp = zeros(length(recon_warped),1);
    for m = 1:length(recon_warped)-1    % -1???
        SQI_tmp(m) = norm(harms(:,m),2)^2 / ( norm(harms(:,m),2)^2 + ...
            sum(PVP_noise( (m-1)*basicTF.hop+1 : m*basicTF.hop ).^2)/basicTF.hop);
    end

    epoch_num = floor((length(PVP_dtr)-sampling_rate*epoch_len)./(shift_n*basicTF.hop));
    % epoch_num = floor( (length(sig_dtr)-5*fs) / 30 ) = 25*fs / 30
    % shift_n = 30/5

    SQI = zeros(epoch_num,1);
    avgAHI = zeros(epoch_num,6);
    
    % initial
    s = 0;
    t = sampling_rate/basicTF.hop*epoch_len;
    % t = fs

    for i = 1:epoch_num
        SQI(i)= median(SQI_tmp(s+1:t));
        base = max(norm(PVP_dtr(s*basicTF.hop+1:t*basicTF.hop)), 0.001); % previent nan and 0
        avgAHI(i,:) = max(sqrt(sum(reinterp_sig(:,s*basicTF.hop+1:t*basicTF.hop).^2, 2))/base, -100);
        s = s + shift_n;    % s + 6 samples
        t = t + shift_n;    % sampling_rate*step_len;
    end
    tmp_SQI = SQI;
    avgAHIWarped_list{k} = avgAHI;
    SQIWarped_table(k,:) = tmp_SQI;


    %% ---------- SQI from hop=5 modes ---------
    PVP_dtr = sig_dtr;
    sampling_rate = fs;
    % resample the reconstructed cardiac signal into the samme sampling rate as
    % the PVP signal.
    recon_cardiac = interp1( (1:basicTF.hop:length(PVP_dtr))/sampling_rate, ...
        real(recon), (1:length(PVP_dtr))/sampling_rate, 'pchip')';
    reinterp_sig = interp1( (1:basicTF.hop:length(PVP_dtr))/sampling_rate, ...
        real(rec_sig).', (1:length(PVP_dtr))/sampling_rate, 'pchip')';
    reinterp_sig = reshape(reinterp_sig, 6, []);

    PVP_noise = PVP_dtr - recon_cardiac;
    for m = 1:length(recon)-1
        SQI_tmp(m) = norm(harms(:,m),2)^2 / ( norm(harms(:,m),2)^2 + ...
            sum(PVP_noise( (m-1)*basicTF.hop+1 : m*basicTF.hop ).^2)/basicTF.hop);
    end

    epoch_num = floor((length(PVP_dtr)-sampling_rate*epoch_len)./(shift_n*basicTF.hop));   
    SQI = zeros(epoch_num,1);
    avgAHI = zeros(epoch_num,6);

    % initial
    s = 0;
    t = sampling_rate/basicTF.hop*epoch_len;
    for i = 1:epoch_num
        SQI(i)= median(SQI_tmp(s+1:t));
        base = max(norm(PVP_dtr(s*basicTF.hop+1:t*basicTF.hop)),0.001); %previent nan and 0
        avgAHI(i,:) = max(sqrt(sum(reinterp_sig(:,s*basicTF.hop+1:t*basicTF.hop).^2,2))/base,-100);
        s = s + shift_n;    %sampling_rate/basicTF.hop*step_len;
        t = t + shift_n;    %sampling_rate*step_len;
    end
    tmp_SQI = SQI;
    avgAHI_list{k} = avgAHI;
    SQI_table(k,:) = tmp_SQI.';

    toc
end

%% save
clear reinterp_sig tfrtic1_warped tfr1_warped

save(strcat("./SQI_results_v2/", "SQIwarped_", db, ".mat"), ...
    'SQIWarped_table', 'avgAHIWarped_list', 'modeNum_list');

save(strcat("./SQI_results_v2/", "SQI_", db, ".mat"), ...
    'SQI_table','avgAHI_list', 'modeNum_list');



