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
labels = csvread(strcat("./DATASETS/", db, "/", db, "_label.csv")).';
y = labels;
n_seg = size(y,2);


%% Parameters
storePath = "."; %"E:/2023spring/Liver_S_BP/PTT/";

% For PPG calculation
fs = 64;
hop = 5;

% segment length for sSQI, eSQI, and perfIdx
pad = 0.5*fs;

% PPG-SQI and AHI shift.
PVP_SQI_step = round((pad/fs)*(fs/hop));

%% 2023/12/31: T_examples 50, 108
% 3 AHI, 1 PVP_SQI
x1 = load(strcat('./SQI_results_v2/SQI_', db, '.mat'));
% 1 sSQI, 1 eSQI
x2 = load(strcat('./SQI_results_v2/classical_SQI/cSQI_', db, '.mat'));

% features
X = [];
% label
Y = [];

for seg = 1:n_seg
fprintf("seg %d/%d\n", seg, n_seg);

% Interpolate to 64Hz
tmp_sSQI = x2.sSQI_list{seg}';
tmp_eSQI = x2.eSQI_list{seg}';
% tmp_perfIdx = x2.perfIdx_list{seg}';

tmp_AHI = x1.avgAHI_list{seg}';
tmp_numMode = x1.modeNum_list(seg);
tmp_SQI = x1.SQI_table(seg,:);
tmp_y = y(:,seg)';

%%%% !!!NEED ATTENTION IF YOU CHANGE WINDOW LENGTH or PAD!!!
t0 = (1:length(tmp_y))./fs;    % 0s to 30s, fs=64Hz
% The sampling period of skewness and entropy are 0.5s
t1 = 2 + (0:52)*(pad/fs);    % 2s to 28s, fs=2Hz
% The sampling period of PPG-SQI and AHI are PVP_SQI_step*(hop/fs)
t2 = 2.5 + (0:52)*PVP_SQI_step*(hop/fs);
% t3 = (1+[0:56]*(pad/fs)); % For perfIdx

%%
tmp_sSQI = interp1(t1, tmp_sSQI, t0, 'pchip');
tmp_sSQI(1:(t1(1)*fs)-1) = tmp_sSQI(t1(1)*fs);
tmp_sSQI(t1(end)*fs+1:end) = tmp_sSQI(t1(end)*fs);

tmp_eSQI = interp1(t1, tmp_eSQI, t0, 'pchip');
tmp_eSQI(1:(t1(1)*fs)-1) = tmp_eSQI(t1(1)*fs);
tmp_eSQI(t1(end)*fs+1:end) = tmp_eSQI(t1(end)*fs);

% tmp_perfIdx = interp1(t3,tmp_perfIdx,t0,'pchip');
% tmp_perfIdx(1:(t3(1)*fs)-1) = tmp_perfIdx(t3(1)*fs);
% tmp_perfIdx(t3(end)*fs+1:end) = tmp_perfIdx(t3(end)*fs);

tmp_AHI = interp1(t2, tmp_AHI.', t0, 'pchip').';
tmp_AHI(:,1:(t2(1)*fs)-1) = tmp_AHI(:,t2(1)*fs)*ones([1,t2(1)*fs-1]);
tmp_AHI(:,(t2(end)*fs+1):end) = tmp_AHI(:,t2(end)*fs)*ones([1,size(tmp_AHI,2)-t2(end)*fs]);

tmp_SQI = interp1(t2, tmp_SQI, t0, 'pchip');
tmp_SQI(1:(t2(1)*fs)-1) = tmp_SQI(t2(1)*fs);
tmp_SQI(t2(end)*fs+1:end) = tmp_SQI(t2(end)*fs);

% Temp: Downsample to 2Hz
downsampling_rate = 2;
tmp_sSQI = tmp_sSQI(1:fs/downsampling_rate:end);
tmp_eSQI = tmp_eSQI(1:fs/downsampling_rate:end);
% tmp_perfIdx = tmp_perfIdx(1:fs/downsampling_rate:end);
tmp_AHI =  tmp_AHI(:,1:fs/downsampling_rate:end);
tmp_SQI = tmp_SQI(1:fs/downsampling_rate:end);
tmp_y = floor(median(buffer(tmp_y, fs/downsampling_rate),1));

% Truncation to 5~25s
trim_time = 5;
tmp_sSQI = tmp_sSQI(downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);
tmp_eSQI = tmp_eSQI(downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);
% tmp_perfIdx = tmp_perfIdx(downsampling_rate*trim_time+1:end - downsampling_rate*trim_time);
tmp_AHI = tmp_AHI(:,downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);
tmp_SQI = tmp_SQI(downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);
tmp_y = tmp_y(downsampling_rate*trim_time+1 : end-downsampling_rate*trim_time);

% Appending
% X = [X;tmp_SQI',tmp_AHI',tmp_sSQI',tmp_eSQI',tmp_perfIdx'];
% X = [X;tmp_SQI',tmp_AHI',tmp_sSQI',tmp_eSQI'];
X = [X; tmp_SQI', tmp_AHI', tmp_sSQI', tmp_eSQI', tmp_numMode*ones(size(tmp_sSQI))'];
% X = [X; tmp_SQI', tmp_sSQI', tmp_eSQI', tmp_numMode*ones(size(tmp_sSQI))'];
Y = [Y;tmp_y'];

% Save 
if mod(seg,100)==0
    save(strcat('./params_v2/features&label_', db, '.mat'), 'X', 'Y');
    csvwrite(strcat('./params_v2/features_', db, '.csv'), X);
    csvwrite(strcat('./params_v2/labels_', db, '.csv'), Y);
end

end
% Save 
save(strcat('./params_v2/features&label_', db, '.mat'), 'X', 'Y');
csvwrite(strcat('./params_v2/features_', db, '.csv'), X);
csvwrite(strcat('./params_v2/labels_', db, '.csv'), Y);


