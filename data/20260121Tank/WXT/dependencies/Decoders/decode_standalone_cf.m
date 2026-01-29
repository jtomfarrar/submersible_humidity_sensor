function decode_standalone_cf( sa_file, DATADIR )
% 2018/10/18 bgreenwood; integrate into uop_plot, assume big-endian (CF)
% 2016/12/19 minor change to meta structure
% 2010/08/12 Seb Bigorre's idea, vtypes/vnames are now cells
% 2010/07/14 fixed sonic wind tilt (was uchar, is char)
% 2003/08 original version ngalbraith

infile = [sa_file '.DAT'];
outfile= [sa_file '.mat'];

meta.decoder.program = mfilename('fullpath'); 
meta.decoder.runtime = datestr(now);
meta.decoder.version = '2018/10/18';
meta.decoder.authors = {'Nan Galbraith, ngalbraith@whoi.edu','Ben Greenwood, bgreenwood@whoi.edu'};
meta.UOP_website = 'http://uop.whoi.edu';
meta.instrument.detailed_information = 'http://uop.whoi.edu/UOPinstruments/frodo/asimet';

% Infer instrument type from filename
if strncmp(lower(infile),'swnd',4)
    inst_type = 'swnd';
else
    inst_type = lower(infile(1:3));
end

nvars=0;vnames=[];
switch lower(inst_type)
    case  {'bpr'}
        nvars=1; vnames= {'bpr' };
        sfact =  1 ;
        vtypes={'float'};
        % 8 bytes end: 4 bytes unused, 2 used, 2 crc
        nbytes  = 8 + (60 * nvars * 4) + 8;
        numunused=4;
        meta.instrument.manufacturer='Heise';
        meta.instrument.model='DXD';
        meta.instrument.range ='850 to 1050 mbar';
        meta.instrument.accuracy='UOP lab calibration, 0.2mbar, manufacturer, 0.02% F.S.';
        meta.data.bpr.long_name='Barometric Pressure';
        meta.data.bpr.unit='mbar';
    case  {'hrh'}
        nvars=2; vnames= {'hrh ' ;'atmp'};
        sfact=[ 1  1 ];
        vtypes={'float'; 'float'};
        % 24  bytes end:  20 unused, 2 used, 2 crc
        nbytes  = 8 + (60 * nvars * 4) + 24;
        numunused=20;
        meta.instrument.manufacturer='Rotronic';
        meta.instrument.model='MP-101A';
        meta.instrument.range='0 to 100 %RH -40 to +60 degC';
        meta.instrument.accuracy='UOP lab calibration 1 %RH, 0.05degC';
        meta.data.hrh.long_name='Relative humidity';
        meta.data.hrh.unit='%RH';
        meta.data.atmp.long_name='Air Temperature';
        meta.data.atmp.unit='degC';
    case  {'lwr'}
        nvars=4; vnames = {'dome '; 'body '; 'tpile'; 'lwr  '};
        sfact=[.01 .01 1 .1];
        vtypes={'short'; 'short'; 'float'; 'short'};
        % short int,   short int,     float,   short
        %      2         2               4               2  ==  10
        % 4 bytes end: no unused, 2 used, 2 crc
        nbytes  = 8 + (60 * 10) + 4;
        numunused=0;
        meta.instrument.manufacturer='Eppley Precision';
        meta.instrument.model='Precision Infrared Radiometer (PIR)';
        meta.instrument.range='0 to 700 W/m^2';
        meta.instrument.resolution='0.1 W/m^2';
        meta.instrument.accuracy='2 W/m^2 (Colbo & Weller,2008)';
        meta.data.lwr.long_name='longwave radiation flux';
        meta.data.lwr.unit='W/m^2';                
        meta.data.dome.long_name='Dome Temperature';
        meta.data.dome.unit='K';
        meta.data.body.long_name='Body Temperature';
        meta.data.body.unit='K';
        meta.data.tpile.long_name='Thermopile voltage';
        meta.data.tpile.unit='volts';
    case  {'prc'}
        nvars=1; vnames= {'prc' };
        sfact= 1 ;
        vtypes={'float'};
        % 8 bytes end:  4 unused, 2 used, 2 crc
        nbytes  = 8 + (60 * nvars * 4) + 8;
        numunused=4;
        meta.instrument.manufacturer='RM Young';
        meta.instrument.model='50202 Self-siphoning rain gauge';
        meta.instrument.range='0 to 50mm';
        meta.instrument.resolution='0.1mm';
        meta.instrument.accuracy='UOP lab calibration 1mm/hr(Serra et al.,2001)';
        meta.data.prc.long_name='Rain Accumulation';
        meta.data.prc.unit='mm';
    case  {'spn'}
        nvars=2; vnames= {'swrt'; 'sdif' };
        sfact=[ 1  1 ];
        vtypes={'float';'float'};
        % 8 bytes end: 20 unused, 2 used, 2 crc
        numunused=20;
        nbytes  = 8 + (60 * nvars * 4) + numunused + 4;
    case  {'swr'}
        nvars=1; vnames= {'swr' };
        sfact= 1 ;
        vtypes={'float'};
        % 8 bytes end: 4 unused, 2 used, 2 crc
        nbytes  = 8 + (60 * nvars * 4) + 8;
        numunused=4;
        meta.instrument.manufacturer='Eppley';
        meta.instrument.model='Precision Spectral Pyranometer (PSP)';
        meta.instrument.range='0 to 2800 W/m^2';
        meta.instrument.resolution='0.1 W/m^2';
        meta.instrument.accuracy='2 W/m^2(Colbo & Weller,2008)';
        meta.data.swr.long_name='short wave radiation flux';
        meta.data.swr.unit='W/m^2';
    case  {'wnd'}
        nvars=8; vnames = {'wnde   '; 'wndn   '; 'wspd   '; 'spdmax '; ...
                'vane   '; 'compass'; 'tiltx  '; 'tilty  '  };
        sfact=[ .01 .01  .2  .2  .1 .1  .2 .2];
        vtypes={'short';'short';'uchar';'uchar';'short';'short';'schar';'schar'};
        % data bytes: 2 2  1   1    2  2  1   1
        % 12 bytes end: 8 unused, 2 used, 2 crc
        nbytes  = 8 + (60 * 12) + 12;
        fprintf('wnd %d bytes\n',nbytes);
        numunused=8;
        meta.instrument.manufacturer='RM Young';
        meta.instrument.model='5103 Wind monitor';
        meta.instrument.range='0 to 60m/s(wind speed), 0 to 360 deg(direction)';
        meta.instrument.resolution='0.01m/s, 0.1deg';
        meta.instrument.accuracy='UOP lab calibration 1%(wind speed), 3 degrees(direction)';
        meta.data.wspd.long_name='Wind Speed';
        meta.data.wspd.methodology='Scalar averaged wind speed over one minute, 8bit resolution';
        meta.data.wspd.units='m/s';
        meta.data.wnde.long_name='Wind TO the East (Oceanographic convention)';
        meta.data.wnde.units='m/s';
        meta.data.wndn.long_name='Wind TO the North (Oceanographic convention)';
        meta.data.wndn.units='m/s';
        meta.data.spdmax.long_name='the maximum 5 second averaged Wind Speed for the minute, 8bit resolution';
        meta.data.spdmax.units='m/s';
        meta.data.compass.long_name='vector average of 11 (5 second average) compass measurements for each minute';
        meta.data.compass.units='degree';
        meta.data.tilt_x.long_name='the average (of 11 samples) Tilt X for the minute';
        meta.data.tilt_x.units='degree';
        meta.data.tilt_y.long_name='the average (of 11 samples) Tilt Y for the minute';
        meta.data.tilt_y.units='degree';
        meta.data.vane.long_name='vane direction';
        meta.data.vane.methodology='last vane direction for the minute';
        meta.data.vane.units='degree';
    case  {'swnd'}
        nvars=10; vnames = {'wnde   '; 'wndn   '; 'wspd   '; 'spdmax '; ...
                'lxydir '; 'compass'; 'tiltx  '; 'tilty  ';'svel   '; 'temp   '};
        sfact=[ .01 .01  .2  .2  .1 .1  .2 .2 1. 1.];
        vtypes={'short';'short';'uchar';'uchar';'short';'short';'schar';'schar';'float';'float'};
        % data bytes: 2 2  1   1    2  2  1   1 4 4
        % 4 bytes end: 0 unused, 2 used, 2 crc
        nbytes  = 8 + (60 * 20) + 4;
        numunused=0;
        meta.instrument.manufacturer='Gill Instruments';
        meta.instrument.model='WindObserver II Ultrasonic Anemometer';
        meta.instrument.range='0 to 65 m/s (wind speed), 0 to 360 deg (direction)';
        meta.instrument.resolution='0.01 m/s, 0.1 deg';
        meta.instrument.accuracy='Manufacturer spec, 2%, 2 degrees';
        meta.data.wnde.long_name='Wind East';
        meta.data.wnde.units='m/s';
        meta.data.wndn.long_name='Wind North';
        meta.data.wndn.units='m/s';
        meta.data.wspd.long_name='Wind Speed';
        meta.data.wspd.units='m/s';
        meta.data.wspd.methodology='Scalar averaged wind speed over one minute';
        meta.data.spdmax.long_name='Maximum Wind Speed';
        meta.data.spdmax.methodology='Maximum 5 second wind speed (m/s). Measurements sampled at 40Hz, 195 samples in a 5s interval';
        meta.data.lxydir.long_name='Last XY Direction (like vane)';
        meta.data.lxydir.units='';
        meta.data.compass.long_name='Last Compass Direction';
        meta.data.compass.units='degree';
        meta.data.tilt_x.long_name='the average (of 11 samples) Tilt X for the minute';
        meta.data.tilt_x.units='degree';
        meta.data.tilt_y.long_name='the average (of 11 samples) Tilt Y for the minute';
        meta.data.tilt_y.units='degree';
        meta.data.svel.long_name='Speed of Sound';
        meta.data.svel.units='';
        meta.data.temp.long_name='Gill Temperature';
        meta.data.temp.units='Celsius';
    otherwise
        disp('Standalone (Compact Flash) Unknown instrument type.');
        disp('Fatal error. type ''dbquit'' and try again');
        keyboard
end

usedloc=nbytes-3;  % this byte indicates the record is "used"
% Open the file
% 'ieee-be' or 'b' IEEE floating point with big-endian byte ordering
% 'ieee-le' or 'l' IEEE floating point with little-endian byte ordering
% open all files as big endian; just change read of floats in cflash
infid=fopen(infile,'rb','b');

nread=0;   nhours=0;
used = 165;  unused = 0;
% isused1 = used;  isused2 = used;
% first, loop through to get array sizes:
while(1),
    [A,count]=fread(infid,nbytes,'uchar');
    if count == nbytes,
        if A(usedloc) ~= used,
            break;
        end
    else            % short record, or empty
        break;
    end
    nhours = nhours + 1;
end
fclose (infid);
fprintf('Processing standalone %s; read %d hourly records beginning ', infile, nhours);
if nhours > 0,
    % re-open, this time to read it for real
    infid=fopen(infile,'rb','b');
    % set up arrays
    nmins=nhours*60;
    mday=ones(nmins,1);
    for iiivv=1:nvars,
        eval([ deblank(vnames{iiivv}) ' =  mday; ' ]);  % cheat w/ empty mday array
    end
    % read the file, in a loop
    for irec=1:nhours,
        % read 6 bytes of date: hh,mm,ss,day,dow,mon
        [A,count]=fread(infid,6,'uchar');
        if ~isempty(A),
            nread=nread + 1;
            hr = A(1);
            minit = A(2);
            day = A(4);
            mon = A(6);
            % read 2 bytes year
            [year,count]=fread(infid,1,'uint16');
            for imin=0:59,
                oindx=(irec-1)*60  + imin + 1;
                mday(oindx)= datenum(year,mon,day,hr,imin,0);
            end
            if irec == 1
                fprintf('%d/%d/%d %d:%d.\n', year, mon,day,hr,minit); end
            % where will this data go in output array?
            irange= ((irec-1)*60 +1) : (irec-1)*60 + 60;
            % read 60 minutes of data, length depends on variable
            for iiivv=1:nvars,
                [B,count]=fread(infid,60,vtypes{iiivv});
                if count == 60,
                    B = B * sfact(iiivv);
                    cmd=sprintf('%s(%d:%d) = B;',...
                        deblank(vnames{iiivv}), irange(1),irange(end));
                    eval(cmd);
                else
                    fprintf('unexpected end reading data at hour %d var %d\n',...
                        irec, iiivv);
                    break;
                end
            end
            % read unused bytes, if any for this inst_type
            if numunused > 0,
                [dummy,count]=fread(infid,numunused,'uchar');
            end
            % read  trailer unsigned bytes. sp2 and sp3 should be A5 hex (165)
            [C,count]=fread(infid,4,'uchar');
            % isused1 = C(1);  %don't need these anymore
            % isused2 = C(2);
        end  % end isempty test
    end      % end nhours loop
    fclose(infid);
    % generate yearday, using year of first record
    if length(mday) > 1
      [stdv]=datevec(mday(1));
      styr=stdv(1);
      yday = mday - datenum(styr,1,0);
    end
    % truncate the data, if styday & endyday set
    if exist('styday','var') && exist('endyday','var') ,
        if endyday > styday,
            yindx = find( yday > styday & yday  < endyday );
            fprintf('using %d values from %d total\n',...
                length(yindx),length(yday));
            fprintf('using %.2f to %.2f from %.2f to %.2f\n',...
                yday(yindx(1)), yday(yindx(length(yindx))),...
                yday(1),yday(length(yday)));
            for iiivv=1:nvars,
                tcmd=sprintf('%s = %s(yindx); ', vnames{iiivv}, vnames{iiivv});
                % disp(tcmd)
                eval(tcmd);
            end
            mday=mday(yindx);
            yday=yday(yindx);
        end
    end
else
    fprintf('NO DATA!\n');
%         clear mday
    nvars=0;
end
vnames=strrep(vnames,' ','');
save([DATADIR 'processed/' outfile],vnames{:},'mday','meta');
copyfile(infile,DATADIR);
end