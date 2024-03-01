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
fr_first = 0.025;   % for RD
hop = 5;    % samples
HighFreq = 18.5;   % Hz
LowFreq = 0.5;   % Hz
window1 = 10;   % seconds

% Amplitude (for reconstruction)
[h, ~, ~] = hermf(fs*window1+1, 1, 5);
h0 = h(floor(size(h,2)/2)+1);

% signal extension
HOP = 1;
extSEC = 5 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets
method.name = 'edmd' ;
method.param = 100 ;

% PPG epoch is set to be 5 seconds for SQI and Cn
epoch_len = 5;

%% save the reconstructed modes
IF_list = cell([n_seg, 1]);    % IFcurves
modes_fund_list = cell([n_seg, 1]);    % one fund_est with fs = fs
modes_dws5_list = cell([n_seg, 1]);    % all modes_est with fs = fs/5

%%
for k = 1:n_seg
% k = 50;
    tic
    fprintf("seg %d/%d\n", k, n_seg);
    % ----------- Bandpass
    [b, a] = butter(6, [0.5, 20]/(fs/2)); % bandpass filter setting.
    sig_dtr = filtfilt(b, a, sig(:,k));% bandpass filter
    
    % Reduce the signal by extending extSEC on both side
    xxsig = SigExtension(sig_dtr, fs, HOP, extK, extM, extSEC, method);

    %% ----------- multi-RD to get the fundmental IF ------------
    basicTF.win = fs*window1+1;
    basicTF.hop = hop;
    basicTF.fr = fr;
    basicTF.HighFreq = HighFreq;
    basicTF.LowFreq = LowFreq;

    % ---------- RD: curve 1 & 2 [c_tmp] ------------
    % TFR with hop = [fs] samples
    [~, ~, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(xxsig, ...
        basicTF.LowFreq/fs, basicTF.HighFreq/fs, ...
        fr_first/fs, fs, basicTF.win, 1, 5, 1, 1, 0);
    % Truncate the TF of extended part
    tfrsq = tfrsq(:, extSEC+1:end-extSEC);

    ver = 2;
    tt = [0 size(tfrsq,2)];
    % lambda = [0.1 0.1];
    % mu = .005;
    lambda = [1/fs 1/fs];
    mu = .07;
    fundRng = [.7 3];
    [c_tmp] = MultiCurveExt(tfrsq, tfrsqtic*fs, lambda, mu, ...
        2, fundRng, 1.0, ver, tt, [-.5 .5]);

    % ----------
    if PLOT
        figure; hold on;
        imageSQ((0:size(tfrsq,2)-1)./fs*hop, tfrsqtic*fs, abs(tfrsq), 0.99);
        axis xy; colormap(1-gray); %colorbar
        xlabel('time(sec)'); ylabel('frequency(Hz)'); set(gca, 'FontSize', 20);
        for l = 1:2
            plot((0:length(c_tmp(:,1))-1)./fs*hop, tfrsqtic(c_tmp(:,l))*fs, 'r-');
        end
        hold off;
        xlim([0 30])
        %figure; plot((0:length(sig_dtr)-1)./fs, sig_dtr);
    end

    % TFR with hop = 5 samples
    [~, ~, tfrsq5, ~, tfrsqtic5] = ConceFT_sqSTFT_C(xxsig, ...
        basicTF.LowFreq/fs, basicTF.HighFreq/fs, ...
        basicTF.fr/fs, hop, basicTF.win, 1, 5, 1, 1, 0);
    tfrsq5 = tfrsq5(:, extSEC*fs/hop+1:end-extSEC*fs/hop);

    IF_fs1hz = tfrsqtic(c_tmp)*fs;
    IF_fs1hz_interp = interp1(1:(fs/hop):size(tfrsq5,2), IF_fs1hz, 1:size(tfrsq5,2), "linear");

    % ----------
    if PLOT
        figure; hold on;
        imageSQ((0:size(tfrsq5,2)-1)./fs*hop, tfrsqtic5*fs, abs(tfrsq5), 0.99);
        axis xy; colormap(1-gray); %colorbar
        xlabel('time(sec)'); ylabel('frequency(Hz)'); set(gca, 'FontSize', 20);
        for l = 1:2
            plot((0:size(IF_fs1hz_interp,1)-1)./fs*hop, IF_fs1hz_interp(:,l), 'r-');
            plot((0:size(IF_fs1hz_interp,1)-1)./fs*hop, IF_fs1hz_interp(:,l)*1.3, 'b-');
            plot((0:size(IF_fs1hz_interp,1)-1)./fs*hop, IF_fs1hz_interp(:,l)*0.7, 'b-');
        end
        hold off;
        xlim([0 30])
        %figure; plot((0:length(sig_dtr)-1)./fs, sig_dtr);
    end

    %% dws=5 modes
    cALL = [];
    tfr_curveExt = zeros(size(tfrsq5));
    for l = 1:6
        if l <= 2
            refupper = IF_fs1hz_interp(:,l)*1.2;
            reflower = IF_fs1hz_interp(:,l)*0.8;
        else
            refupper = (l+0.3)*tfrsqtic5(cALL(:,1))*fs;
            reflower = (l-0.3)*tfrsqtic5(cALL(:,1))*fs;
        end

        for m = 1:size(tfrsq5,2)
            lowerbdd = max( 1 , round( (reflower(m)-basicTF.LowFreq)/fr ) );
            upperbdd = min( length(tfrsqtic5) , round( (refupper(m)-basicTF.LowFreq)/fr ) );
            tfr_curveExt(lowerbdd:upperbdd, m) = tfrsq5(lowerbdd:upperbdd, m);
        end
        cALL = [cALL CurveExt(abs(tfr_curveExt).', 1.0)];
        tfr_curveExt = zeros(size(tfrsq5));
    end

    % IF curves [Hz]
    IF_dws5 = tfrsqtic5(cALL)*fs;
    % harmonics modes
    rec_sig = zeros(6, size(tfrsq5,2));
    for l = 1:6
        rec_sig(l,:) = Recon_sqSTFT_v2(tfrsq5, tfrsqtic5, fs, cALL(:,l), .2, h0);
    end
    % recon = sum(rec_sig);

    % ----------
    if PLOT
        figure; hold on;
        imageSQ((0:size(tfrsq5,2)-1)./fs*hop, tfrsqtic5*fs, abs(tfr_curveExt), 0.99);
        axis xy; colormap(1-gray); %colorbar
        xlabel('time(sec)'); ylabel('frequency(Hz)'); set(gca, 'FontSize', 20);
        for l = 1:2
            plot((0:size(IF_fs1hz_interp,1)-1)./fs*hop, tfrsqtic5(cALL(:,l))*fs, 'r-');
        end
        hold off;
        xlim([0 30])
        %figure; plot((0:length(sig_dtr)-1)./fs, sig_dtr);
    end

    IF_list{k}.fs1hz = IF_fs1hz;
    IF_list{k}.dws5 = IF_dws5;

    modes_dws5_list{k} = rec_sig;


    %% dws=1 modes
    [~, ~, tfrsq1, ~, tfrsqtic1] = ConceFT_sqSTFT_C(xxsig, ...
        basicTF.LowFreq/fs, basicTF.HighFreq/fs, ...
        fr/fs, 1, basicTF.win, 1, 5, 1, 1, 0);
    % Truncate the TF of extended part
    tfrsq1 = tfrsq1(:, extSEC*fs+1:end-extSEC*fs);  % fs = fs

    % c is a dws=5 curve IF curve
    c = cALL(:,1);
    c = repmat(c, 1, hop);
    c = reshape(c', 1, [])';

    tfr_curveExt = tfrsq1;
    ref = tfrsqtic(c)*fs;
    for m = 1:size(tfrsq1,2)
        % ----- trim the TFR for CurveExt -----
        lowbdd = max(1, round( (ref(m)*0.75-basicTF.LowFreq)/fr ));
        upbdd = min(size(tfrsq1,1), round( (ref(m)*1.25-basicTF.LowFreq)/fr ));
        tfr_curveExt([1:lowbdd upbdd:end], m) = 0;
    end
    cALL1 = CurveExt(abs(tfr_curveExt)', 1.0);

    IF_list{k}.dws1 = tfrsqtic1(cALL1)*fs;

    %clear tfr_curveExt lowbdd upbdd

    % --------- estimate number of modes ---------
    fund_est = Recon_sqSTFT_v2(tfrsq1, tfrsqtic1, fs, cALL1, 0.2, h0);
    modes_fund_list{k} = fund_est;

    % ---------
%     if PLOT
%         figure; hold on;
%         imageSQ((0:size(tfrsq1,2)-1)./fs, tfrsqtic1*fs, sqrt(abs(tfrsq1)), 0.99);
%         axis xy; colormap(1-gray); %colorbar
%         xlabel('time(sec)'); ylabel('frequency(Hz)'); set(gca, 'FontSize', 20);
%         plot((0:length(cALL)-1)./fs, tfrsqtic1(cALL)*fs, 'r-');
%         hold off;
%         title('hop = 1 curve');
%         xlim([0 30])
%     end
    toc
end

clear tfrsq tfrsq5 tfrsq1 tfrsqtic tfrsqtic5 tfrsqtic1 tfr_curveExt

save(strcat("./SQI_results_v2/modes_", db, ".mat"), ...
    'IF_list', 'modes_fund_list', 'modes_dws5_list');
