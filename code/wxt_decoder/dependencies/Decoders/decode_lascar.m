% 2018/07/19 Ben Greenwood; integrate into uop_plot.m
% 2016/07/16 s.bigorre revises code to loadlascar_v2.m
%            Code revised to accomodate new date format (yyyy-mm-dd)
% 2010/01/01 n.galbraith issues loadlascar.m

infid=fopen(infile,'rt'); % defined in uop_decode.m
sdata=textscan(infid,...
    '%d %.0f-%.0f-%.0f %.0f:%.0f:%.0f %f %f %f %s %s',...
    'Headerlines',1,'Delimiter',',');
fclose(infid);

meta.decoder.version = '2016/07/14';
meta.decoder.author = 'Sebastien Bigorre';
meta.decoder.runtime = datestr(now);
meta.decoder.program = mfilename('fullpath');
meta.instrument = 'Lascar';
meta.instrument_SN = sdata{12}(1);
meta.variables = 'Temperature(degC),RH(%),Dew Point(degC)';
meta.comment = 'Note: this is raw data';

% SPURS2,Time,Celsius(?C),Humidity(%rh),Dew Point(?C),Serial Number
% 1,2016-07-11 17:00:00,24.0,50.5,13.1,010032184

year=sdata{2}; 
month=sdata{3};
day=sdata{4};
hour=sdata{5};
minute=sdata{6};
second=sdata{7};

Lascar.mday=datenum(year,month,day,hour,minute,second);
Lascar.atmp=sdata{8};
Lascar.hrh=sdata{9};
Lascar.dewpt=sdata{10};



