%%  History
% 2019/10/02 Ben Greenwood; modify search criteria to identify first line of Rotronic data
% 2018/04/30 Ben Greenwood; create script

meta.decoder_version = '2019/10/02';
meta.decoder_author = 'Ben Greenwood; bgreenwood@whoi.edu';
meta.decoder_runtime = datestr(now);
meta.instrument = 'Rotronic HC2A';
meta.variables = 'Temperature(degC) and RH(%)';
meta.comment = 'Note: this is raw data';
    
%% Open the file
fin=fopen(infile,'r');
if (fin ~= -1)
    fprintf('Rotronic file %s opened\n',infile);
else
    fprintf('Unable to open file: %s\n',infile);
end

try
    % identify beginning of Rotronic.csv data section
    while(1)
        line = fgetl(fin);
        if strncmp(line,'[#D]',4)
            break;
        end
        if feof(fin)
            fprintf('Error, data not found in %s\n',infile);
            return;
        end
    end

    n=1;
    while(1)
        data = fscanf(fin,'%d/%d/%d,%d:%d:%d%*c%cM,%f,%f',9);
        data(3)=data(3)+2000;
        data(4)=mod(data(4),12); % convert times 12:00 to 0:00 
        if data(7) == 80 % ASCII 'p' found, convert from 12 hour to 24:00 format
            data(4) = data(4) + 12;
        end
        mday(n) = datenum(data(3),data(1),data(2),data(4),data(5),data(6));
        hrh(n)=data(8);
        atmp(n)=data(9);
        n = n + 1;
    end
catch ME
    if feof(fin)
        fprintf('%d Records found\n',n);
        if n > 0
            stat=1;
        end
    else
        fprintf('Unknown exception, error reading record #%d, byte %d\n',n,ftell(fin));
        fprintf('line #%d: %s\n',ME.stack(1).line,ME.message);
        keyboard;
    end
end
fclose(fin);

