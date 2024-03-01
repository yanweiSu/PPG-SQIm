function [PAT,ECG_idx] = PAT_calc(R_ECG,R_PPG,Fs)
    ECG_idx = [];
    %PPG_idx = [];
    PAT = [];
    for i = 1:length(R_ECG)
        tmp = min(find(R_PPG > R_ECG(i)));

        if isempty(tmp)
            break
        end

        if (R_PPG(tmp) - R_ECG(i))/Fs > .35 || (R_PPG(tmp) - R_ECG(i))/Fs < 0.12
            continue
            %(R_PPG(tmp) - R_ECG(i))/Fs;
        else
            (R_PPG(tmp) - R_ECG(i))/Fs;
            PAT = [PAT,(R_PPG(tmp) - R_ECG(i))/Fs];
            ECG_idx = [ECG_idx, i];
            %PPG_idx = [PPG_idx, R_PPG(tmp)];
        end

        R_PPG = R_PPG(tmp:end); 
    end
    
end