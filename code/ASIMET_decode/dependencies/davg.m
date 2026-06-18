function [mday_avg, dval] = davg(mday,rawvar,minutes)
%   Ben Greenwood - script to average DIRECTION variable (degrees) over user specified time period (minutes)
%   Variable is broken into x/y components, each components is binned using histc(),
%   and averaged using accumarray(). Finally, components are vector combined.
%   sampling_period = raw.mday(2)-raw.mday(1);
%   avg_period = minutes/24/60;
%   bin_size = floor(avg_period/sampling_period);
%   bin_count = floor(length(raw.mday)/bin_size);
    N = cos(rawvar*pi/180)';
    E = sin(rawvar*pi/180)';
    
    bin_size = minutes/24/60;
    tmin = min(mday);
    tmax = max(mday);
    t0 = tmin-(minutes/24/60/2);
    tf = tmax+(minutes/24/60);
    bin_edges = t0:minutes/24/60:tf;
    
    [b,idx] = histc(mday,bin_edges);
    if size(N,1) == 1 % transpose row to column if necessary
        N=N';
        E=E';
    end
    N_avg = accumarray(idx,N,[],@mean);
    E_avg = accumarray(idx,E,[],@mean);
    mday_avg = (tmin:(minutes/24/60):tmax)';
    if length(N_avg)>length(mday_avg) % sometimes one more bin than variable length
        N_avg(end)=[];
        E_avg(end)=[];
    end
    dval = mod(atan2(E_avg,N_avg)*180/pi,360);
end