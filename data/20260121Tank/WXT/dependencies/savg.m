function [mday_avg, sval] = savg(mday,rawvar,minutes)
%   Ben Greenwood - script to average SCALAR over user specified time period (minutes)    
%   data is binned using histc() and averaged using accumarray().
%   sampling_period = raw.mday(2)-raw.mday(1);
%   avg_period = minutes/24/60;
%   bin_size = floor(avg_period/sampling_period);
%   bin_count = floor(length(raw.mday)/bin_size);

    bin_size = minutes/24/60;
    tmin = min(mday);
    tmax = max(mday);
    t0 = tmin-(minutes/24/60/2);
    tf = tmax+(minutes/24/60);
    bin_edges = t0:minutes/24/60:tf;
    
    [b,idx] = histc(mday,bin_edges);
    if size(rawvar,1) == 1 % transpose row to column if necessary
        rawvar=rawvar';
    end
    sval = accumarray(idx,rawvar,[],@mean);
    mday_avg = (tmin:(minutes/24/60):tmax)';
    if length(sval)>length(mday_avg) % sometimes one more bin than variable length
        sval(end)=[];
    end
end
