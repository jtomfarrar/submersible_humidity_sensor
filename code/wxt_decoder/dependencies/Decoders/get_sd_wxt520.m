function [stat] = get_sd_wxt520(infile, rawfile, outfile, meta, p, average_winds)
% Load binary wxt520 file from SD card, decode, and save to .mat file
% Upon successful processing of a file, get_sd_wxt520 will return stat=1

% Required fields:
%  infile, outfile, meta (struct)
%  rawfile optionally

%% input fields:
%--------------------------------------------------------------------------
%
% infile - raw dat file eg 'VWX003_dock.dat'
% rawfile - intermediate file for matlab array 'VWX_3arr.mat'; %'VWXT520_arr.mat';
% outfile - destination matfile  eg 'VWXT003_v2.mat';
% meta - whatever metadata you like.
%     Recommended MINIMUM fields:
%     meta.site = 'WHOTS8 burnin'; %'NTAS burnin';
%     meta.deployment = '8';
%     meta.start_year = 2011; %2010;
% p - user specified parameters
%
%--------------------------------------------------------------------------
%
%%  Modifications
% 2020/08/27 B.Greenwood when bad records are removed, ensure indices don't reference removed records
% 2020/06/05 B.Greenwood remove records with bad timestamps. Identified as records out of sequence.
% 2019/04/03 B.Greenwood Fix problem with wdir. Previously, compass measurements
%   were not added to wdir. Also update meta description of wind variables
% 2019/04/03 B.Greenwood Include wind averaging ability
% 2018/05/14 B.Greenwood remove code which rotates winds by compass. UOP group decided to leave winds unrotated
%   for all decoders with the rational being that decoders should report as close to raw values as possible.
% 2018/02/01 B.Greenwood modify variable names to match UOP plot_1mn_burnin.m script
% 2018/01/29 B.Greenwood new 272 byte records, nest all freads inside try block, add EOF and CRC exceptions.
%   create throughput report and add it to meta.history. Move existing CRC lookup code into
%   locate_record_beginning() function. Use fseek to place cursor at beginning of record rather than skip.
%   insert firmware strings into meta-data
% 2017/09/13 N.Galbraith keep trying, expect trash at start of file
% 2016/04/20 N.Galbraith read SD card - new format including new time struct
% 2015/12/03 N.Galbraith ADD COMPASS AND "DIRECTION" before generating
%   wind vectors, per Seb's discussion
% 2014/09/16 N.Galbraith pass max wind to this function as
%   meta.history.decode.windspeedlimits, to window based on speed
% 2014/09/11 N.Galbraith make funct, avg winds, add some metadata, fix version date
% 2014/08/28 N.Galbraith slight clean up, (comments)
% 2014/08/27 N.Galbraith fix indexing problem that was overwriting some files
% 2014/08/14 N.Galbraith fix up metadata per UOP standards and re-name bpr_pa to bpr,
%   since it is apparently in mbars and make mday npts,1 (like the rest of the variables)
% 2013/05/14 N.Galbraith use daterange (matlab mdays) for start and end dates
%   will truncate if this variable exists

    fprintf('Reading WXT520 raw file: %s\n',infile);
    meta.decoder_version='2018/05/14';
    meta.decoder_author={'Nan Galbraith ngalbraith@whoi.edu','Ben Greenwood bgreenwood@whoi.edu'};
    stat=0;
    spdlims=[-999 999];
    used = 165;
    unused = 0;

    %% Open the file
    % 'ieee-be' or 'b' IEEE floating point with big-endian byte ordering
    % 'ieee-le' or 'l' IEEE floating point with little-endian byte ordering
    infid=fopen(infile,'rb','l');
    nread=0;

    locate_record_beginning(infid,used);
    
    % bgreenwood estimate # records and pre-allocate space for data
    file = dir(infile);
    nrecs = ceil(file.bytes / 272);
    fprintf('Using file-size (%d bytes) to estimate the number of records: %d\n',file.bytes,nrecs);
    B = ones(nrecs,48)*NaN; % Pre-allocate space for data
    S = ones(nrecs,64)*NaN;
    
    for irec=1:nrecs
        try 
            time = fread(infid,6,'uchar');  % First 8 bytes contain date-stamp
            year = fread(infid,1,'uint16'); % sec(uchar),min(uchar),hour(uchar),day(uchar),dow(not used),month(uchar),year(uint16)
            mday(irec) = datenum(year,time(6),time(5),time(3),time(2),time(1));
            flag(irec) = 0; % flag used to identify bad timestamps

            fread(infid,6,'uchar');    % record byte-size = 272 (ASCII) 
            fread(infid,1,'uint16');   % record byte-size = 272 (uint16)

            B(irec,:) = fread(infid,48,'float');    % read 48 floats of data
            S(irec,:) = fread(infid,64,'uchar');    % read 64 bytes of characters
            if S(irec,63) ~= used
                throw(MException('MATLAB:mycode','CRC_checksum_fail'));
            end

        catch ME
            if feof(infid)
                fprintf('EOF after record # %d, byte %d\n',irec-1,ftell(infid));
                break;
            elseif strcmp(ME.message,'CRC_checksum_fail')
                fprintf('Bad checksum; record #%d, byte %d\n',irec,ftell(infid)-272);
                locate_record_beginning(infid,used);
            else
                fprintf('Unknown exception, error reading record #%d, byte %d\n',irec,ftell(infid));
                fprintf('line #%d: %s\n',ME.stack(1).line,ME.message);
                keyboard;
            end
        end
    end
    fclose(infid);

	% 2020/08/27 BG mday array may be shorter than nrecs at this point
	% 2020/06/05 BG Identify bad timestamps for removal
    for i = 2:length(mday)
        if mday(i) < mday(i-1) % bad records contain stale data including timestamp in past
            flag(i)=1;
        end
    end    
    
    save(rawfile,'mday','B','S'); % save raw data file
    
    % Truncate data to date-range specified in do_load_wxt520.m
    iok = find(mday >= p.datemin & mday <= p.datemax & ~isnan(mday) & ~flag);
    
    % Check CRC
    crc = S(iok,63:64); % set to 0xA5A5 upon record write
    crc_bad = find( crc(:,1) ~= 165 | crc(:,2) ~= 165);
    if crc_bad
        iok = setdiff(iok,crc_bad); % remove records with bad CRC
        for i = crc_bad
            fprintf('bad CRC; record # %d\n',i);
        end
    end
    
    % unpack WXT520 data
    mday = mday(iok)';
    wdir11 = B(iok,1:11);
    wspd11 = B(iok,12:22);
    wspd_min = B(iok,23);
    wspd_max = B(iok,24);
    compass11 = B(iok,25:35);
    tilt_x = B(iok,36);
    tilt_y = B(iok,37);
    atmp = B(iok,38);
    hrh = B(iok,39);
    bpr = B(iok,40);
    precip = B(iok,41);
    rain_duration = B(iok,42);
    rain_intensity = B(iok,43);
    hail_accumulation = B(iok,44);
    hail_duration = B(iok,45);
    hail_intensity = B(iok,46);
    rain_peak_intensity = B(iok,47);
    hail_peak_intensity = B(iok,48);
    version = S(iok,1:20);
    brdversion = S(iok,21:36);
    modser = S(iok,37:40);
    sensor = S(iok,41:48);
    spare = S(iok,49:58);
    samp_count = S(iok,59);
    wndflag = S(iok,60);
    rhtpflag = S(iok,61);
    prcflag = S(iok,62);
    
    fprintf('%s calculating winds\n', datestr(now));
    wnde=ones(length(mday),1)*NaN; % pre-assign, fil with NaNs
    wndn=wnde; wspd=wnde; wdir=wnde; wnumav=wnde;
    ngdwnds=0;
    for ii=1:length(mday)
        tspd=wspd11(ii,:);
        % 2019/04/03 Ben Greenwood WXT raw winds use meteorological convention
        % 2018/05/14 BG UOP group decided to leave wind direction unrotated for all UOP decoders.
        % note that WXT will be off by 180 with respect to ASIMET (which uses oceanographic convention)
        % tdir=wdir11(ii,:)+compass11(ii,:); % 20151203 N.Galbraith add compass!!! to wind dir raw
        % tdir=rem(tdir,360);
        tdir=mod(wdir11(ii,:)+compass11(ii,:),360);
        gv=find(tspd < spdlims(2) );
        if gv
            [wnde11,wndn11] = md2vect(tspd(gv),tdir(gv));
            gv=find(wnde11 < spdlims(2) & wndn11 < spdlims(2));
            wnumav(ii)=length(gv);
            if gv
                wnde(ii) = mean(wnde11(gv));
                wndn(ii) = mean(wndn11(gv));
                [wspd(ii),wdir(ii)] = vect2md(wnde(ii),wndn(ii));
                ngdwnds=ngdwnds+1;
            end
        end
        % average compass11
        [cmpu11,cmpv11]=md2vect(ones(1,11),compass11(ii,:));
        cmpu(ii)=mean(cmpu11);
        cmpv(ii)=mean(cmpv11);
        [~,compass(ii)]=vect2md(cmpu(ii),cmpv(ii));
    end
    compass=compass';
    fprintf('%s done with wind, %d good 1minute records\n', datestr(now),ngdwnds);

    %% Conditionally average winds
    averaged_vars='';
    if (average_winds)
        [mday_avg,wnde_avg]=savg(mday,wnde,average_winds);
        [mday_avg,wndn_avg]=savg(mday,wndn,average_winds);
        [mday_avg,wspd_avg]=savg(mday,wspd,average_winds);
        [mday_avg,compass_avg]=davg(mday,compass,average_winds);
        [mday_avg,wdir_avg]=davg(mday,wdir,average_winds);
        averaged_vars='mday_avg wnde_avg wndn_avg wspd_avg compass_avg wdir_avg ';
    end
    
    fprintf('Throughput report\n');
    interval=[ mday(1) ];
    dctemp.record_intervals=[];
    for i = 2:length(mday)
        if mday(i)-mday(i-1) > (5/24/60) % 2020/06/05 BG increase threshold to 5 minutes
            fprintf('Interval %s to %s\n',datestr(interval(end)),datestr(mday(i-1)));
            dctemp.record_intervals=[ dctemp.record_intervals ; [ interval(end) mday(i-1) ] ];
            interval = [interval mday(i)];
        end
    end
    fprintf('Interval %s to %s\n',datestr(interval(end)),datestr(mday(i)));
    dctemp.record_intervals=[ dctemp.record_intervals ; [ interval(end) mday(i-1) ] ];
            
    % this is wordy; I'm trying to get a well-formed history.decode
    % where one field may have been set in the calling program
    dctemp.date=datestr(now,'yyyy-mm-dd');
    dctemp.infile  = [pwd '/' infile];
    dctemp.outfile = [pwd '/' outfile];
    dctemp.process='unpack WXT520 binary, assign variables and dates';
    if isfield(meta,'history')
        if isfield(meta.history,'decode')
            if isfield(meta.history.decode,'datelimits')
                dctemp.datelimits=datestr(meta.history.decode.datelimits);
                daterange = meta.history.decode.datelimits;
            end
            if isfield(meta.history.decode,'windspeedlimits')
                dctemp.windspeedlimits=meta.history.decode.windspeedlimits;
                spdlims=dctemp.windspeedlimits;

            end
        end
    end
    meta.history=dctemp;
    
    % instrument
    meta.instrument.manufacturer='Vaisala';
    meta.instrument.model='WXT520';
    meta.instrument.model_vers='SD';
    meta.instrument.software_vers=char(version(1,:));       % Firmware version string extracted in raw binary
    meta.instrument.software_brd=char(brdversion(1,:));     % Board version string extracted from raw binary
    meta.instrument.model_sn=char(modser(1,:));             % Module serial number
    meta.instrument.sensor_sn=char(sensor(1,:));            % Sensor serial number

    meta.data.wdir.long_name='Wind Direction';
    meta.data.wdir.methodology='Vector average of 11 five-second raw wind direction measurements; meteorological convention. Adjusted using compass';
    meta.data.wdir.units='degree';
    meta.data.wspd.long_name='Wind Speed';
    meta.data.wspd.methodology='Average of 11 five-second scalar wind speed measurements';
    meta.data.wspd.units='m/s';
    meta.data.wnde.long_name='Wind East';
    meta.data.wnde.methodology='Derived from raw wdir and wspd measurements; meteorological convention';
    meta.data.wnde.units='m/s';
    meta.data.wndn.long_name='Wind North';
    meta.data.wndn.methodology='Derived from raw wdir and wspd measurements; meteorological convention';
    meta.data.wndn.units='m/s';
    if isfield(meta.history,'windspeedlimits')
        meta.data.wndu.minmaxvalues=meta.history.windspeedlimits;
        meta.data.wndv.minmaxvalues=meta.history.windspeedlimits;
    end

    meta.data.wspd_min.long_name='the minimum Wind Speed for the minute';
    meta.data.wspd_min.units='m/s';
    meta.data.wspd_max.long_name='the maximum Wind Speed for the minute';
    meta.data.wspd_max.units='m/s';
    meta.data.compass.long_name='vector average of 11 five-second compass measurements';
    meta.data.compass.units='degree';
    meta.data.compass11.long_name='11 five-second samples (one minute) of compass data';
    meta.data.compass11.units='degree';
    meta.data.tilt_x.long_name='the average (of 11 samples) Tilt X for the minute';
    meta.data.tilt_x.units='degree';
    meta.data.tilt_y.long_name='the average (of 11 samples) Tilt Y for the minute';
    meta.data.tilt_y.units='degree';
    meta.data.atmp.long_name='Air Temperature';
    meta.data.atmp.units='Degree C';
    meta.data.hrh.long_name='Relative Humidity';
    meta.data.hrh.units='%';
    meta.data.bpr.long_name='Air Pressure';
    meta.data.bpr.units='mbar';
    meta.data.precip.long_name='Rain Accumulation';
    meta.data.precip.units='mm';

    meta.data.wdir11.long_name='11 five-second samples (one minute) of Wind Direction data, not corrected for compass';
    meta.data.wdir11.units='degree';
    meta.data.wspd11.long_name='11 five-second samples (one minute) of Wind Speed data';
    meta.data.wspd11.units='m/s';

    meta.data.rain_duration.long_name='Rain Duration';
    meta.data.rain_duration.units='seconds';
    meta.data.rain_intensity.long_name='Rain Intensity';
    meta.data.rain_intensity.units='mm/h)';
    meta.data.hail_accumulation.long_name='Hail Accumulation';
    meta.data.hail_accumulation.units='hits/cm2';
    meta.data.hail_duration.long_name='Hail Duration';
    meta.data.hail_duration.units='sec';
    meta.data.hail_intensity.long_name='Hail Intensity';
    meta.data.hail_intensity.units='hits/cm2hr';
    meta.data.rain_peak_intensity.long_name='Rain Peak Intensity';
    meta.data.rain_peak_intensity.units='mm/hr';
    meta.data.hail_peak_intensity.long_name='Hail Peak Intensity';
    meta.data.hail_peak_intensity.units='hits/cm2hr';
    
    %meta.mday='Time stamp using datenum in Matlab';

    myvars=' meta mday wnde wndn wdir wspd  wdir11 wspd11 wspd_min wspd_max compass compass11 tilt_x tilt_y ';
    myvars=[myvars 'atmp hrh bpr precip rain_duration '];
    myvars=[myvars 'rain_intensity hail_accumulation hail_duration hail_intensity rain_peak_intensity hail_peak_intensity'];
    myvars=[myvars ' version brdversion modser sensor wndflag rhtpflag prcflag ' averaged_vars];
    try
        savecmd=sprintf('save %s %s',outfile,myvars);
        eval (savecmd);
    catch saveerr
        disp(['savecmd is ' savecmd]);
        rethrow(saveerr);
        keyboard
    end
    stat=1;
end

function locate_record_beginning(infid,used)
    for nbskp=1:100000
        % not unused 'count' variable message suppressed, to keep this code
        % viable for older versions of Matlab NRG 20140911
        [A1,count]=fread(infid,1,'uchar'); %#ok<*NASGU>
        if  A1 == used,
            %fprintf('A1 at %d\n',ftell(infid)-1);
            [A2,count]=fread(infid,1,'uchar');
            if  A2 == used,
                %fprintf('A2 at %d\n',ftell(infid)-1);
                break;
            end
        end
    end
    fseek(infid,-272,'cof'); % bgreenwood; move file position to beginning of record
    fprintf('Next record begins at byte %d \n', ftell(infid));
end
        
        
function [vx,vy] = md2vect(mag,dir)
    % converts speed and heading to vectors east and north (vx,vy)
    jk=find(dir > 360);
    dir(jk)=dir(jk) - 360;
    jk=find(dir < 0);
    if ~isempty(jk)
        dir(jk)=(dir(jk))+360;
    end
    dir=dir*pi/180; % to radians
    vx = mag .* sin(dir);
    vy = mag .* cos(dir);
end

function [mag,dir] = vect2md(vx,vy)
    %function [mag,dir] = vect2md(vx,vy);
    % converts vectors east and north (vx,vy) to speed and heading
    mag = sqrt(vx.^2 + vy.^2);
    dir = atan2(vx,vy)*180/pi;
    jk=find(dir<0);
    if ~isempty(jk)
        dir(jk)=360+(dir(jk));
    end
end
