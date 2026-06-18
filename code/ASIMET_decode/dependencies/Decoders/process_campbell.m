% Read Campbell 1 minute data (TOA5 format)
% Ben Greenwood 2019/03/20
function Campbell = read_Campbell_T0A5(infile,DATADIR,average_winds)
ctable = readtable(infile);
ctable(1,:)=[];
ctable(1,:)=[];
cstruct=table2struct(ctable,'ToScalar',true);
mday = datenum(cstruct.TIMESTAMP);
RECORD = str2double(cstruct.RECORD);
posix = str2double(cstruct.posix);
atmp = str2double(cstruct.Rotronic_temp);
hrh = str2double(cstruct.Rotronic_RH);
if isfield(cstruct,'PTB210')
    bpr = str2double(cstruct.PTB210);
elseif isfield(cstruct,'PTB110')
    bpr = str2double(cstruct.PTB110);
end
lwr = str2double(cstruct.SGR4_flux_Wm2_Avg);
lwr_volts = str2double(cstruct.SGR4_extVolts);
lwr_btemp = str2double(cstruct.SGR4_bodyTempC);
swr = str2double(cstruct.SMP21_flux_Wm2_Avg);
if isfield(cstruct,'SMP21_flux2_Wm2_Avg')
  swr2 = str2double(cstruct.SMP21_flux2_Wm2_Avg);
end
swr_volts = str2double(cstruct.SMP21_extVolts);
swr_btemp = str2double(cstruct.SMP21_bodyTempC);
wnde = str2double(cstruct.WndE_Avg);
wndn = str2double(cstruct.WndN_Avg);
wdir = str2double(cstruct.Wdir60);
wspd = str2double(cstruct.Wspd_Avg);
compass = str2double(cstruct.Comp60);
vane = str2double(cstruct.Vane60);
if isfield(cstruct,'PRC')
    precip = str2double(cstruct.PRC);
elseif isfield(cstruct,'PRC_Avg');
    precip = str2double(cstruct.PRC_Avg);
end
if isfield(cstruct,'SBE37_temp')
    sst = str2double(cstruct.SBE37_temp);
else
    sst = nan(length(mday),1);
end
if isfield(cstruct,'SBE37_cond')
    cond = str2double(cstruct.SBE37_cond);
else
    cond = nan(length(mday),1);
end
LAT = str2double(cstruct.gps_lat);
LON = str2double(cstruct.gps_lon);
if isfield(cstruct,'PS200_Avg_3_')
    load_mA = str2double(cstruct.PS200_Avg_3_);
else isfield(cstruct,'load_mA_Avg')
    load_mA = str2double(cstruct.load_mA_Avg);
end
if isfield(cstruct,'PS200_Avg_1_')
    BATT = str2double(cstruct.PS200_Avg_1_);
elseif isfield(cstruct,'VBatt_Avg')
    BATT = str2double(cstruct.VBatt_Avg);
end
temp_CR1000x = str2double(cstruct.PTemp_Avg);
if isfield(cstruct,'PS200_Avg_6_')
    temp_CH200 = str2double(cstruct.PS200_Avg_6_);
elseif isfield(cstruct,'Chg_TmpC_Avg')
    temp_CH200 = str2double(cstruct.Chg_TmpC_Avg);
end
bad=find(load_mA>10000);
load_mA(bad)=NaN;

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


meta.decoder_version='2019/05/29';
meta.decoder_author='Benjamin Greenwood, bgreenwood@whoi.edu';
meta.instruments = xml_read('campbell_instruments.xml');
meta.variables = xml_read('campbell_variables.xml');
meta.processed = datestr(now);

outdir = [ DATADIR 'processed/' strrep(infile,'.DAT','.mat') ];
save(outdir);
copyfile(infile,DATADIR);
