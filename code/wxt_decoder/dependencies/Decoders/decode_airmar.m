% 2022/06/01 Benjamin Greenwood decode WHOTS18 Campbell Airmar internal ASCII data
% files expected to be in format Airmarxxx.dat numerically increasing

BASE = '.';

files = dir( [ BASE '/*.dat']);
[~,findx] = sort( str2double( regexp( {files.name}, '\d+', 'match', 'once')));
files = files(findx);


first = 1;

for f = 1:length(files)
  fprintf(' * processing %s\n', files(f).name);

  % BG use readtable() to parse comma delimited data - appears to correctly interpret column names
  data = readtable(files(f).name);

  % BG Initialize Airmar cell array if this is the first file processed
  if first == 1
    Airmar.mday = [];
    for i = 2:21
      Airmar.(data.Properties.VariableNames{i}) = [];
    end
    first = 0;
  end

  % BG convert column 1 timestamp > datetime > datenum and append to Airmar cell array
  dt = datetime(data{:,1},'format','yyyy-MM-dd HH:mm:ss');
  Airmar.mday = [ Airmar.mday ; datenum(dt) ];

  % BG append additional columns
  for i = 2:21
    cname = data.Properties.VariableNames{i}; % column name
    Airmar.(cname) = [ Airmar.(cname) ; data{:,i} ];
  end
end

meta.decoder.program = mfilename('fullpath');
meta.decoder.version = '2022/06/01';
meta.decoder.runtime = datestr(now);
meta.decoder.author = 'Benjamin Greenwood, bgreenwood@whoi.edu';
meta.UOP_group_website = 'uop.whoi.edu';
meta.instrument.manufacturer = 'Airmar';
meta.instrument.model = 'WX200';
meta.units = 'Winds reported in m/s, Air temp in degC, pressure in mbar, SOG in km/hr, magvar/heading/lat/lon/pitch/roll in degrees';
