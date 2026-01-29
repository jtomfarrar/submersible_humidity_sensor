% 2022/06/13 Ben Greenwood plot legend use Latex interpreter, user '.' marker for non-averaged Airmar vars
% 2022/06/01 Ben Greenwood add Airmar vars to plots
% 2019/04/03 Ben Greenwood Convert WXT winds from meteorological to oceanographic convention
% 2019/02/27 Ben Greenwood include SD module decoder data into plots
% 2018/05/24 Ben Greenwood create plot script

% load in decoded data (assume appropriately named .mat files exist
% L={};SA={};R={};W_sd={};sbe39={};C={};Lascar={};SD={};
Instruments = { Loggers,Campbell_files,Standalone_CF,Standalone_SD, ...
    Rotronic_files, WXT_SD_files,sbe39ATfiles,Lascar_files };
% data = {L, SA, R, W_sd, sbe39, C, Lascar, SD};
data = {};
for i = 1:length(Instruments)
    for j=1:length(Instruments{i})
        data{i}{j} = load([ DATADIR '/processed/' Instruments{i}{j} '.mat']);
        dateok = find(data{i}{j}.mday <= tend & data{i}{j}.mday >= tstart );
    end
end

if exist('Airmar')
	Instruments{7} = {'Airmar'};
	if (average_winds)
    	[~,Airmar.wspd_avg]=savg(Airmar.mday,Airmar.wspd,average_winds);
		indx = find(Airmar.wspd_avg==0); % remove zero records caused by Airmar duty cycling
		Airmar.wspd_avg(indx)=nan;
    	[~,Airmar.wnde_avg]=savg(Airmar.mday,Airmar.wnde,average_winds);
		indx = find(Airmar.wnde_avg==0);
		Airmar.wnde_avg(indx)=nan;
    	[~,Airmar.wndn_avg]=savg(Airmar.mday,Airmar.wndn,average_winds);
		indx = find(Airmar.wndn_avg==0);
		Airmar.wndn_avg(indx)=nan;
	    [~,Airmar.compass_avg]=davg(Airmar.mday,Airmar.compass,average_winds);
		indx = find(Airmar.compass_avg==0);
		Airmar.compass_avg(indx)=nan;
	    [Airmar.mday_avg,Airmar.wdir_avg]=davg(Airmar.mday,Airmar.wdir,average_winds);
		indx = find(Airmar.wdir_avg==0);
		Airmar.wdir_avg(indx)=nan;
	end
	data{7}{1} = Airmar;
end

% 2019/04/03 Ben Greenwood Convert WXT winds from meteorological to oceanographic convention
if length(data) >= 6
  for i = 1:length(data{6})
    data{6}{i}.wdir = mod( data{6}{i}.wdir + 180, 360 );
    data{6}{i}.wdir_avg = mod( data{6}{i}.wdir_avg + 180, 360 );
    data{6}{i}.wnde = data{6}{i}.wnde * -1;
    data{6}{i}.wnde_avg = data{6}{i}.wnde_avg * -1;
    data{6}{i}.wndn = data{6}{i}.wndn * -1;
    data{6}{i}.wndn_avg = data{6}{i}.wndn_avg * -1;
  end
end

vectors_avgd={'wnde','wndn','wspd','vane','compass'};

% Make summary plots
vars{1} = {'hrh','atmp','swr','lwr'};
units{1} = {'%','^oC','W/m^2','W/m^2'};
vars{2} = {'sst','cond','precip','bpr'};
units{2} = {'^oC','S/m','mm','mbar'};
vars{3} = {'wnde','wndn','wspd','wdir','compass'};
units{3} = {'m/s','m/s','m/s','deg','deg'};
vars{4} = {'wnde_avg','wndn_avg','wspd_avg','wdir_avg','compass_avg'};
units{4} = {'m/s','m/s','m/s','deg','deg'};

for p =1:4
    clf;
    la_list=[];
    set(gcf,'PaperPositionMode','Manual');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[8.5 11]);
    set(gcf,'PaperPosition',[ 0 0 8.5 11 ]);
    for v = 1:length(vars{p})
        ax{v}=subplot(length(vars{p}),1,v);
        legend_str=[];
        for i = 1:length(Instruments)
            for j = 1:length(Instruments{i})
              %fprintf(' -- inst: %s, marker: %s\n',Instruments{i}{j},marker);
                if isfield(data{i}{j},vars{p}{v})
                    hold on;
                    if strfind(vars{p}{v},'_avg')
                        plot(data{i}{j}.mday_avg,data{i}{j}.(vars{p}{v}));
                        legend_str{end+1}=[Instruments{i}{j} '\_avg' num2str(average_winds)];
                    else
                        if strcmp(Instruments{i}{j},'Airmar')
                            marker='.'; % since Airmar is duty cycled, use discrete marker
                        else
                            marker='-';
                        end
                        plot(data{i}{j}.mday,data{i}{j}.(vars{p}{v}),marker);
                        legend_str{end+1}=[Instruments{i}{j}];
                    end
                end
            end
        end

        if length(legend_str)
            legend(legend_str,'Interpreter','latex');
            la_list{end+1}=v; % only link axes for plots which contain data, otherwise problems
        end
        title(vars{p}{v},'Interpreter','none');
        ylabel(units{p}{v});
        set(gca,'fontsize',7);
        datetick('x',6);
    end

    suptitle({experiment,deployment});
    if length(la_list)
        linkaxes([ax{ [la_list{:}] }],'x');
        savefig([DATADIR 'figs/' experiment '_' deployment '_summary' num2str(p)]);
    end
end

% Make individual plots
vars = {'hrh','atmp','swr','lwr','sst','cond','precip','bpr','wnde','wndn','wspd','vane','compass','wdir','wdir_avg','wnde_avg','wndn_avg','wspd_avg','vane_avg','compass_avg'};
units = {'%','^oC','W/m^2','W/m^2','^oC','S/m','mm','mbar','m/s','m/s','m/s','deg','deg','deg','deg','m/s','m/s','m/s','deg','deg'};
figure(1);
for p=1:length(vars)
    clf;
    hold on;
    legend_str=[];
    for i = 1:length(Instruments)
        for j = 1:length(Instruments{i})
            if isfield(data{i}{j},vars{p})
                hold on;
                if strfind(vars{p},'_avg')
                    plot(data{i}{j}.mday_avg,data{i}{j}.(vars{p}));
                    legend_str{end+1}=[Instruments{i}{j} '\_avg' num2str(average_winds)];
                else
                    if strcmp(Instruments{i}{j},'Airmar')
                        marker='.'; % since Airmar is duty cycled, use discrete marker
                    else
                        marker='-';
                    end
                    plot(data{i}{j}.mday,data{i}{j}.(vars{p}),marker);
                    legend_str{end+1}=Instruments{i}{j};
                end
            end
        end
    end

    if length(legend_str)
        legend(legend_str,'Interpreter','latex');
        title({experiment,deployment,vars{p}},'Interpreter','none');
        ylabel(units(p));
        datetick('x',6);
        savefig([DATADIR 'figs/' experiment '_' deployment '_' vars{p}] );
    end
end
close all;


