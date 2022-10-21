% last modified 2021/09/27; BG add SWR ASCII text output

%*********     User defined     *******************************************************************%
experiment = 'WND306';
deployment = '1';
dumpname   = 'WND306';              % directory where all output will go
tstart = datenum(2022,1,1);         % Time range of plots
tend = datenum(2022,1,1);
latitude = 41.52;		            % lat/lon only required to process Loggers
longitude = -70.67;		            % necessary to accurately calculate practical salinity
average_winds = 10;                 % number of minutes to average winds over (zero = no avg)
%%*********     Sensors defined  ******************************************************************%
% ASIMet Loggers
Loggers = {};           			% ASIMet Loggers; must have .DAT extension (do not include here)
Standalone_SD = {'WND306'};			% ASIMET standalone (SD version) Format = {'PRC273','SWR242'}
Standalone_CF = {}; 				% ASIMET standalone Compact Flash. Format= {'SWR206','BPR200'}
WXT_SD_files = {};					% WXTs (SD version); must have .DAT extension (do not include here)
Rotronic_files = {};				% Rotronic files; must have .csv extension (do not include here)
sbe39ATfiles = {};					% SBE39 files; must have .asc extension (do not include here)
Campbell_files = {};				% Campbell CR1000X 1min .DAT file; (do not include .DAT extension here)
Vaisalafile ='';					% WXT (traditional non SD version)
WXT_SN = '';						% dont put zeros
Lascar_files = {}; 					% Lascar; must have .csv extension (do not include here)
%%*************************************************************************************************%


%% Decoding and Plotting below
addpath(genpath('dependencies')); % recursively add dependencies subdirectory to path
cwd=pwd;
DATADIR = fullfile(cwd,filesep,dumpname,filesep); %fullfile(cwd,filesep);
if exist(DATADIR)~=7
    mkdir(DATADIR)
    mkdir(DATADIR,'figs');
    mkdir(DATADIR,'processed');
end

%% Process Rotronic
for i = 1:length(Rotronic_files)
    meta=[];
    infile = [ Rotronic_files{i} '.csv' ];
    decode_Rotronic;
    save([DATADIR 'processed/' Rotronic_files{i} '.mat'],'mday','hrh','atmp','meta');
    copyfile(infile,DATADIR);
end

%% Process Standalone_CF
for i = 1:length(Standalone_CF)
    meta=[];
    decode_standalone_cf(Standalone_CF{i},DATADIR);
end

%% Process Standalone_SD
for i =1:length(Standalone_SD)
    meta=[];
    decode_standalone_sd(Standalone_SD{i},DATADIR,average_winds);
end

%% Process Lascar
for i = 1:length(Lascar_files)
	meta=[];
	infile = [ Lascar_files{i} '.csv' ];
	decode_lascar;
    save([DATADIR 'processed/' Lascar_files{i} '.mat'],'Lascar','meta');
	copyfile(infile,DATADIR);
end

%% Process Campbell (2019)
for i = 1:length(Campbell_files)
    meta=[];
    infile = [ Campbell_files{i} '.DAT' ];
    process_campbell(infile,DATADIR,average_winds);
end

%% PROCESS ASIMet Logger Data
%for i = 1:length(Standalone)
%    SNs{i} = str2num(Standalone{i}(2:3));
%end
minmd=tstart;
maxmd=tend;
for i = 1:length(Loggers)
    meta=[];
    dat_file=[Loggers{i} '.DAT'];
    asc_file=[DATADIR 'processed/' Loggers{i} '.asc'];
    mat_file=[DATADIR 'processed/' Loggers{i} '.mat'];
    imet_bin_to_asc;
    imet_asc_to_mat;
    copyfile(dat_file,[DATADIR]);
 end

%% Process WXT-SD module
for i = 1:length(WXT_SD_files)
    [~,f] = fileparts(WXT_SD_files{i});
    infile = [ WXT_SD_files{i} '.DAT'];
    rawfile = [ DATADIR 'processed/' f '_raw.mat']; %'VWXT520_arr.mat';
    outfile = [ DATADIR 'processed/' f '.mat'];
    meta=[];
    meta.wxt_sn= WXT_SD_files{i};
    meta.deployment = experiment; 
    meta.start_year = deployment; 
    p.datemin = tstart; % minimum date to consider
    p.datemax = tend;  % maximum date to consider
    [stat] = get_sd_wxt520(infile, rawfile, outfile, meta, p, average_winds);
    copyfile(infile,DATADIR);
end

%% Process SBE39s
for i=1:length(sbe39ATfiles)
    meta=[];
    infile=sbe39ATfiles{i};
    get_sbe_39;
    copyfile([infile '.asc'],DATADIR);
end

%% Now call plotting routine
plot_new;
copyfile('uop_decode.m',DATADIR);
copyfile('README',DATADIR);
copyfile('dependencies',[DATADIR 'dependencies']); % recursively include code in output directory 


