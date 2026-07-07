%add meta.data fields for variables in an asimet data file
% atmp  compass    hrh        lwr        precip     swr        wnde       wndn_corr  wspd            
% body       cond       latitude   mday       sal        tpile      wnde_corr  wsmax      year       
% bpr        dome       longitude  meta       sst        vane       wndn       wsmin  

% vers Aug 2014, removing the rotated winds
% ngalbraith

meta.data.mday.long_name='Matlab days with UTC time';
meta.data.mday.units='days';

meta.data.yday.long_name='Year day with UTC time';
meta.data.mday.units='days';
meta.data.yday.convention='Naval convention: Jan 1 at 12:00 is yday 1.5';

meta.data.latitude.long_name='latitude of measurement, degrees north';
meta.data.latitude.units='degrees';

meta.data.longitude.long_name='longitude of measurement, degrees east';
meta.data.longitude.units='degrees';

meta.data.year.long_name= 'start year'; 
  
meta.data.atmp.long_name= 'air temperature';
meta.data.atmp.units= 'degree_C';    
meta.data.bpr.long_name= 'barometric pressure';
meta.data.bpr.units= 'millibar';
meta.data.hrh.long_name= 'relative humidity';
meta.data.hrh.units= 'percent';
meta.data.swr.long_name= 'shortwave radiation';
meta.data.swr.units= 'W m-2';
meta.data.lwr.long_name= 'longwave radiation';
meta.data.lwr.units= 'W m-2';
meta.data.precip.long_name= 'precipitation level';
meta.data.precip.units= 'mm';

meta.data.sst.long_name= 'sea temperature';
meta.data.sst.units= 'degree_C';
meta.data.sal.long_name= 'sea water practical salinity';
meta.data.sal.units='1'; 
meta.data.cond.long_name= 'sea water conductivity';
meta.data.cond.units= 'S m-1';

meta.data.wnde.long_name= 'wind velocity east, uncorrected';
meta.data.wnde.units=  'meters/second';
meta.data.wndn.long_name= 'wind velocity north, uncorrected';
meta.data.wndn.units=  'meters/second';

% meta.data.wnde_corr.long_name= 'wind velocity east with magnetic correction';
% meta.data.wnde_corr.units=  'meters/second';
% meta.data.wndn_corr.long_name= 'wind velocity north with magnetic correction';
% meta.data.wndn_corr.units=  'meters/second';


meta.data.wspd.long_name= 'wind velocity';
meta.data.wspd.units=  'meters/second';
meta.data.wsmax.long_name=  'maximum observed wind velocity';
meta.data.wsmax.units='meters/second';
meta.data.wsmin.long_name= 'minimum observed wind velocity';
meta.data.wsmin.units= 'meters/second';
meta.data.vane.long_name= 'last anemometer vane reading';
meta.data.vane.units= 'degrees';
meta.data.compass.long_name= 'last compass reading';
meta.data.compass.units= 'degrees';

meta.data.dome.long_name= 'LWR done reading';
meta.data.dome.units= 'degree_K';
meta.data.body.long_name=  'LWR case reading';
meta.data.body.units= 'degree_K';
meta.data.tpile.long_name= 'thermopile voltage';
meta.data.tpile.units= 'volts';

% meta.data.wnde_corr.long_name= 'wind velocity east, corrected for magnetic declination';
% meta.data.wnde_corr.units= 'meters/second';
% meta.data.wndn_corr.long_name= 'wind velocity north, corrected for magnetic declination';
% meta.data.wndn_corr.units= 'meters/second';
