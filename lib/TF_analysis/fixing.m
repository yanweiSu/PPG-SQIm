function [tmp] = fixing(tmp,fs)

if nargin<2; fs=100; end

data = buffer(tmp, 6*fs, 3*fs);
qual = range(data(:,2:end));
for seg = find(qual > 0.4)
    tmp(1,(seg-1)*3*fs+1:(seg+1)*3*fs) = NaN;
%     tmp(1,(seg-1)*3*fs+1:seg*3*fs) = NaN;
end

ok = [NaN tmp NaN];
ini = 2;
while ini <= length(tmp)+1
    if ~isnan(ok(ini)) && isnan(ok(ini-1))
        fin = ini;
        while fin <= length(tmp)+1 % find the last index of this segment
            if ~isnan(ok(fin)) && isnan(ok(fin+1))
                break;
            end
            fin = fin+1;
        end
        med = mod(median(ok(1,ini:fin)),1);
        if med < 0.5
            ok(1,ini:fin) = ok(1,ini:fin) - (median(ok(1,ini:fin))-med);
        else
            ok(1,ini:fin) = ok(1,ini:fin) - (median(ok(1,ini:fin))-med+1);
        end
        ini = fin+1;
    else
        ini = ini+1;
    end
end

tmp = ok(1,2:end-1);
% tmp = fillmissing(tmp, 'movmedian', 70*100);

end