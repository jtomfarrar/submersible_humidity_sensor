% imet_asc_to_mat
% load matlab-generated ascii file, save a matfile
% 2019/04/03 Ben Greenwood wind direction not transmitted by logger; derive wind direction from vector components wnde and wndn
% 2019/04/03 Ben Greenwood conditionally calculate wind averages
% vers 2016/11/02
%     Use Gibbs Seawater Toolbox (GSW), gsw_matlab_v3_05_5 for pracitcal salinity
% version 2012/12/20:
%      apply magvar if available;
% updated july 19 2011 to find start_year better
% July 25 2013 add more metadata re magvar applied, if applied
% Aug 3 2016 call add_meta_data to get long names of variables
% requires:
%    asc_file - ascii (input) file name
%    mat_file   - output (mat) file name
%    start_year, or meta.start_year, or year
%    PATH to seawater toolbox must be set
%    latitude: used in salinity calc
% if available:
%  sst_depth_m - otherwise, uses .5
%  minmd, maxmd - used to remove records with wild dates
%  magvar - in degrees - if it exists, vectors are corrected
%  (note raw vectors are kept also)
% clear all derived variables
meta.history.decode.program.name = mfilename('fullpath');
meta.history.decode.program.version_date = '2017/01/18';
meta.history.decode.outputfile=mat_file;

wnde_corr=[];  wndn_corr=[]; sal=[];
disp(['loading ' asc_file ]);
idat = load(asc_file);

fprintf('will save dates from %s to %s\n', datestr(minmd),...
    datestr(maxmd));

%  hr, minu, day, mon, year, rec_no, muxpar 1 2 3 4 5 6 7
hour=idat(:,1);
minute=idat(:,2);
day=idat(:,3);
month=idat(:,4);
year= idat(:,5);
mday = datenum(year, month,day,hour,minute,0);
year=year(1);
% you may want to specify time limits; do that
% in the routine that calls this one
ginx= find(mday > minmd & mday < maxmd);
fprintf('input length %d output %d\n',...
    length(mday),length(ginx));

% wnde, wndn, wsavg, wsmax, wsmin, vane, compass
%  8     9     10       11    12    13    14
wnde    = idat(ginx,8);
wndn    = idat(ginx,9);
wspd    = idat(ginx,10);
wsmax   = idat(ginx,11);
wsmin   = idat(ginx,12);
vane    = idat(ginx,13);
compass = idat(ginx,14);

meta.data.magnetic_variation.applied_to='none';
% Skip the magnetic correction unless magvar exists and ~= 0
% Used in Radians, stored in degrees
if exist('magvar','var')
    if magvar ~= 0
        magvarrad = magvar * ((2*pi)/360);
        wnde_corr = (wnde * cos(magvarrad)) + (wndn * sin(magvarrad));
        wndn_corr = (-1 * wnde * sin(magvarrad)) + (wndn * cos(magvarrad));
        meta.data.magnetic_variation.applied_to='wnde_corr,wndn_corr';
        meta.data.magnetic_variation.applied=magvar;
    else
        disp('not rotating wind, magvar = 0');
    end
else
    disp('not rotating wind, magvar does not exist');
end
% bpr, hrh, atmp, swr
%  15  16   17    18
bpr     = idat(ginx,15);
hrh     = idat(ginx,16);
atmp    = idat(ginx,17);
swr     = idat(ginx,18);
% dome, body, tpile, lwr, prlev, sst, cond
dome    = idat(ginx,19);
body    = idat(ginx,20);
tpile   = idat(ginx,21);
lwr     = idat(ginx,22);
precip   = idat(ginx,23);
sst     = idat(ginx,24);
cond    = idat(ginx,25);
mday    = mday(ginx);

%% calculate practical salinity
mysal=which('gsw_SP_from_C');
if ~isempty(mysal)
    if ~exist('sst_depth_m','var')
        sst_depth_m =.5;
    end
    if ~exist('latitude','var')
        latitude=41.5242;
    end
    if ~exist('longitude','var')
        longitude=-70.6712;
    end
    fprintf('calculating salinity with sensor depth %.2f m\n', sst_depth_m);
    try
        % note GSW uses height (m) not depth
        gswp=gsw_p_from_z(-sst_depth_m, latitude);
        sal=gsw_SP_from_C(cond*10,sst,gswp);
        meta.data.sal.algorithm='Gibbs Seawater Toolbox (GSW), gsw_matlab_v3_05_5';
        meta.data.sal.latitude_used=latitude;
        meta.data.sal.sensordepth_used=sst_depth_m;
        % sal = real(sw_salt(10*cond/sw_c3515, sst, sw_pres(sst_depth_m,latitude)));
    catch sal_err
        disp('can not calculate salinity!')
    end
end

% 2019/04/03 Ben Greenwood wind direction not transmitted by logger; 
% derive wind direction from vector components wnde and wndn
[~,wdir]=vect2md(wnde,wndn);

%% 2019/04/03 Ben Greenwood Conditionally average winds
averaged_vars='';
if (average_winds)
    [mday_avg,wnde_avg]=savg(mday,wnde,average_winds);
    [mday_avg,wndn_avg]=savg(mday,wndn,average_winds);
    [mday_avg,wspd_avg]=savg(mday,wspd,average_winds);
    [mday_avg,vane_avg]=davg(mday,vane,average_winds);
    [mday_avg,wdir_avg]=davg(mday,wdir,average_winds);
    [mday_avg,compass_avg]=davg(mday,compass,average_winds);
    averaged_vars='mday_avg wnde_avg wndn_avg wspd_avg vane_avg compass_avg wdir_avg ';
end


add_meta_data
cmd=['save '  mat_file ];
cmd=[cmd ' latitude longitude year meta ' averaged_vars];
cmd=[cmd ' mday atmp body bpr compass cond dome hrh lwr'];
cmd=[cmd ' precip sst swr tpile vane wnde wndn wsmax wsmin wspd wdir'];
if ~isempty(sal)
    cmd=[cmd ' sal'];
end
if ~isempty(wndn_corr)
    cmd=[cmd ' wndn_corr wnde_corr'];
end

if ~exist('latitude','var')
    cmd=strrep(cmd,'latitude longitude','');
end
fprintf('saving with command %s\n',cmd);
eval(cmd)



