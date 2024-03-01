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
y = csvread(strcat("./DATASETS/", db, "/", db, "_label.csv")).';
n_seg = size(sig,2);

%% Parameters
% For PPG calculation
fs = 64;
% segment length for sSQI and eSQI
len = 4*fs;
pad = 0.5*fs;
len_perf = 2*fs;

%% indices
%skewness
sSQI_list = cell([n_seg, 1]);

%entropy
eSQI_list = cell([n_seg, 1]);

% % perfusion index
% perfIdx_list = cell([size(x.data,1),1]);

%%
for seg = 1:n_seg
% seg = 100;
    fprintf("seg %d/%d", seg, n_seg);
    [b, a] = butter(6, [0.5, 15]/(fs/2));% bandpass filter setting.
    PVP_dtr = filtfilt(b, a, sig(:,seg));% bandpass filter

    rng = prctile(PVP_dtr, 97.5) - prctile(PVP_dtr, 2.5);
    PVP_dtr = PVP_dtr .* (2/rng);
    PVP_dtr = PVP_dtr - prctile(PVP_dtr, 2.5) - 1;
    fprintf(" %.1f, %.1f\n", prctile(PVP_dtr, 97.5), prctile(PVP_dtr, 2.5));

    % figure; plot((1:length(PVP_dtr))./fs, sig(:,seg));
    tmp = buffer(PVP_dtr, len, len-pad)';
    tmp = tmp(len/pad:end, :);
    entropy_list = [];%zeros([size(tmp,1),1]);
    for i = 1:size(tmp,1) % [15 34]
%         figure; plot(tmp(i,:), 'LineWidth', 1.3); set(gca, 'FontSize', 20); ylim([0 1]); xlabel('samples')
%         figure; histogram(tmp(i,:), 'BinWidth', 0.05); set(gca, 'FontSize', 20); xlabel('au'); ylabel('counts'); xlim([0 1])
%         figure; histogram(tmp(i,:), 20); set(gca, 'FontSize', 20); xlabel('au'); ylabel('counts'); xlim([0 1])
%         
%         rng = prctile(tmp(i,:), 97.5) - prctile(tmp(i,:), 2.5);
%         tmp(i,:) = tmp(i,:) .* (2/rng);
%         tmp(i,:) = tmp(i,:) - prctile(tmp(i,:), 2.5) - 1;

        [counts, ~] = histcounts(tmp(i,:), 'BinWidth', 0.05);
        % [counts2, edges2] = histcounts(tmp(i,:), 20);
        % [counts,edges] = histcounts(tmp(i,:)/norm(detrend(tmp(i,:)),2), 'BinWidth', BinWidth);
        counts = counts/length(tmp(i,:));    % pmf
        % counts2 = counts2 / size(tmp,2);
        seq = counts(counts > 0);
        % seq2 = counts2(counts2 > 0);
        % entropy_list(i) = median(-sum((counts+eps).*log(counts+eps), 1)); % prevent Inf
        % entropy_list(i,:) = [-sum(seq.*log(seq)) -sum(seq2.*log(seq2))];
        entropy_list(i) = -sum(seq.*log(seq));
%         entropy_list(i) = median(-sum((counts+1e-10).*log(counts+1e-10)));
    end
    eSQI_list{seg} = entropy_list;

    
%     % 2s to 28s, hop 0.5s: 53pts
%     tt_ent = 2 + (0:52)*0.5;
%     tt_y = (0:1919)*(1/fs);
%     figure; hold on;
%     plot(tt_y, PVP_dtr, 'LineWidth', 1.2);
%     plot(tt_ent, entropy_list-min(min(entropy_list))+1, 'LineWidth', 1.2);
% %     plot(tt_ent, entropy_list(:,2)-min(min(entropy_list))+1, 'LineWidth', 1.2);
% %     plot(tt_ent, entropy_list(:,1), 'LineWidth', 1.2);
% %     plot(tt_ent, entropy_list(:,2), 'LineWidth', 1.2);
%     plot(tt_y, y(:,seg), '-m', 'LineWidth', 1.2);
%     legend("PPG", "entropy", "labels")
%     set(gca, 'FontSize', 16);
%     xlabel('time(sec)')
%     title(strcat(db, " seg ", string(seg)));
%     ylim([-3 3])

    
    %%
    % Calculate SQIs
    % [~, ~, eSQI_list{seg}] = SQI_eval(PVP_dtr, len, pad);
    [sSQI_list{seg}] = sSQI(PVP_dtr, len, pad);
    % [perfIdx_list{seg}] = perfusionIdx(sig(:,seg),PVP_dtr,len_perf,pad);
    
    if mod(seg,100)==0
    % save(strcat(storePath,'/SQI_results/classical_SQI/','classical_SQI_',type,'_',string(id),'.mat'),'sSQI_list','eSQI_list','perfIdx_list');
    save(strcat('./SQI_results_v2/classical_SQI/cSQI_', db, '.mat'), ...
        'sSQI_list', 'eSQI_list');
    end
end
% save(strcat(storePath,'/SQI_results/classical_SQI/','classical_SQI_',type,'_',string(id),'.mat'),'sSQI_list','eSQI_list','perfIdx_list');
save(strcat('./SQI_results_v2/classical_SQI/cSQI_', db, '.mat'), ...
    'sSQI_list', 'eSQI_list');

