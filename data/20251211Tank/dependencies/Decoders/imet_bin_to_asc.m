% imet_bin_to_asc
% read an imet binary file, write out ascii
% requires:
%    dat_file, asc_file
% input files are big-endian, from logger processor
% Dallas Semiconductor DS87C530 microcontroller
% vers 29 Dec  2010, 
% updated 2011/July 19 to add version info to metadata
% 2013/05/03 use  better filename terms (but no version change)
% and just process 1 input file to 1 output file
infid=fopen(dat_file,'rb','b');
ofid=fopen(asc_file, 'wt');

disp(['processing ' dat_file ' to ' asc_file ])
meta.history.decode_program = mfilename('fullpath');
meta.history.decode.program_version_date = '2010/12/29';
meta.history.decode_date = datestr(now);

if (infid)
    % 'ieee-be' or 'b'   IEEE fp with big-endian byte ordering
    % 'ieee-le' or 'l'   IEEE fp with little-endian byte ordering

    % first are 7 time bytes, first five are char
    % next, mux_parm byte
    % now shorts: Ve Vn rotor1 rotor2 LastCompass
    % byte TiltX Tilty
    % short stemp
    % 4 byte float: thermistor resistance, muxed parameter
    % 2 "used" bytes - they should be A5A5
    %   recdefault = struct(...
    %     'hr',0,'min',0,'sec',0,'day',0,'mon',0, 'year',0,'muxpar',0, ...
    %     'wnde',0,'wndn',0,'wsavg',0,'wsmax',0,'wsmin',0, ...
    %     'vane',0, 'compass',0, 'bpr',0, 'hrh',0, 'atmp',0,...
    %     'swr',0,'dome',0,'body',0,'tpile',0,'lwr',0, ...
    %     'prlev',0, 'sst', 0,  'cond', 0,'mday',0 );
    % 64 byte records
    nrecs=1;
    used = 165;
    nbad=0;
    isused(1) = used;  isused(2) = used;
    while (isused(1) == used && isused(2) == used),
        % set these up so they always print (they should, anyway, but ... )
        hr=NaN;minu=NaN; day=NaN; mon=NaN; year=NaN; rec_no=NaN; muxpar =NaN;
        wnde=NaN; wndn=NaN; wsavg=NaN; wsmax=NaN; wsmin=NaN; vane=NaN; compass =NaN;
        bpr=NaN; hrh=NaN; atmp=NaN; swr=NaN;
        dome=NaN; body=NaN; tpile=NaN; lwr=NaN;
        prlev=NaN; sst=NaN; cond=NaN;
        if rem(nrecs,20000) == 0, fprintf(' %d', nrecs); end
        if rem(nrecs,200000) == 0, fprintf('\n'); end
        % read 5 bytes of date
        [UC5,count]=fread(infid,5,'uchar');
        if UC5 == -1 
            disp('eof')
            break
        end
        if count == 5,
            if isempty(UC5), disp('end data'); break; end
           
            hr = UC5(1);
            minu = UC5(2);
            day = UC5(3);
            mon = UC5(4);
            %if ( hr == 255  && minu == 255) ,
            %   disp('unused '); break; end;
            year = UC5(5) + 2000;
            
            % next, rec_num, short int
            [UI1,count]=fread(infid,1,'uint16');
            if count ~= 1 || ...
               isempty(UI1), disp('end data'); break; end

            rec_no  = UI1(1);

            % next, mux_parm byte
            [UC1,count]=fread(infid,1,'uchar');
            if count ~= 1, break; end
            muxpar = UC1(1);

            % next, 2 shorts wnde,wndn
            [I2,count]=fread(infid,2,'int16');
            if count ~= 2, break; end
            wnde  =  I2(1) / 100;
            wndn =  I2(2)  / 100;
            
            % next 3 unsigned shorts, wsavg,wsmax,wsmin
            [UI3,count]=fread(infid,3,'uint16');
            if count ~= 3, break; end
            wsavg  =  UI3(1)/ 100;
            wsmax =  UI3(2)/ 100;
            wsmin =  UI3(3)/ 100;
            
            % next 2 shorts, short vane, compass
            [I2,count]=fread(infid,2,'int16');
            if count ~= 2, break; end
            vane  =  I2(1) / 10;
            compass =  I2(2)  / 10;

            % next unsigned short bpr
            [UI1,count]=fread(infid,1,'uint16');
            if count ~= 1, break; end
            bpr  =  (UI1(1)/ 100) + 900; 
            
            % next short rh;
            [I1,count]=fread(infid,1,'int16');
            if count ~= 1, break; end
            hrh  =  I1(1)  / 100.0 ;
            
            % next unsigned short atmp 
            [UI1,count]=fread(infid,1,'uint16');
            if count ~= 1, break; end
            atmp  =  (UI1(1) / 1000) - 20 ;
            
            % next swr - it is signed!
            [I1,count]=fread(infid,1,'int16');            
            if count ~= 1, break; end
            swr  =  I1(1)/10 ;
            
            % next unsigned  dome,body;
            [UI2,count]=fread(infid,2,'uint16');
            if count ~= 2, break; end
            dome  =  UI2(1)/100 ;
            body  =  UI2(2)/100 ;
            
            % next short tpile lwr prlev
            [I3,count]=fread(infid,3,'int16');
            if count ~= 3, break; end
            tpile  =  I3(1) /10;
            lwr  =  I3(2) /10;
            prlev  =  I3(3)/100 ;

            % next unsigned short sst cnd
            [UI2,count]=fread(infid,2,'uint16');
            if count ~= 2, break; end
            sst  = ( UI2(1) / 1000) -5;
            cond  =  UI2(2) / 10000 ;

            % read and discard 4 shorts 1 long 3 shorts
            % 3rd (last) short is 2 "used " bytes
            [I4,count]=fread(infid,4,'int16');  % 3 bat shorts
            if count ~= 4, break; end
            % opt parm
            [L1 ,count] = fread(infid,1,'int32');  % opt parm, long
            if count ~= 1, break; end
            [I3 ,count] = fread(infid,3,'int16');  % 3 spare shorts
            if count ~= 3, break; end

            % 2 "used" bytes - they are A5A5 in a good record
            [UC2,count]=fread(infid,2,'uchar');
            isused1 = UC2(1);
            isused2 = UC2(2);

            fprintf(ofid,...
                '%2d  %2d  %2d  %2d  %4d  %5d ', ...
                hr, minu, day, mon, year, rec_no, muxpar );
            fprintf(ofid,...
                '%7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %6.2f  %6.2f', ...
                wnde, wndn, wsavg, wsmax, wsmin, vane, compass );
            fprintf(ofid,...
                '%8.2f  %7.3f  %7.3f  %6.1f ', ...
                bpr, hrh, atmp, swr );
            fprintf(ofid,...
                ' %7.3f  %7.3f  %7.1f  %6.1f  %7.2f  %7.3f  %7.4f\n',...
                dome, body, tpile, lwr, prlev, sst, cond);
            nrecs=nrecs + 1;
        else
            nbad=nbad +1;
            if nbad > 100
                break
            end
        end  % end rec length check (only first 5 fields)
    end      % end processing loop
    fprintf('\ndone\n');
    fclose(infid);
    fclose (ofid);
    % decrement counter  (?)
    nrecs=nrecs - 1;
    fprintf('read %d records\n',nrecs);
end  % end if infid

