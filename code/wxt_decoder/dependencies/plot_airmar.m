% 2022/06/01 BG integrate into uop_decode.m script; designed for WHOTS18 burnin
% 2021/05/28 BG original standalong Airmar plotting script designed for WHOTS17 burnin

% Only consider time range tstart < t < tend
good_date = find( (tstart < Airmar.mday) & (Airmar.mday < tend));
vars = {'mday','Airmar_ATMP','Airmar_BPR','Airmar_COG','Airmar_Pitch','Airmar_RH','Airmar_Roll','Airmar_SOG','Airmar_WDIR','Airmar_WSPD','Airmar_magvar','Battery','PanelTemp','Posix','RECORD','gps_height','gps_lat','gps_lon','gps_sats','gps_time'};
for i = 1:length(vars)
  Airmar.(vars{i}) = Airmar.(vars{i})(good_date);
end

% MET
figure(1);
sgtitle([ experiment ' Airmar MET'],'Interpreter','none');
subplot(5,1,1);
plot(Airmar.mday,Airmar.Airmar_ATMP,'.');
legend('ATMP');
datetick('x');
ylabel('degC');
%xlim([tstart,tend]);
subplot(5,1,2);
plot(Airmar.mday,Airmar.Airmar_BPR,'.');
legend('BPR');
datetick('x');
ylabel('mbar');
%xlim([tstart,tend]);
subplot(5,1,3);
plot(Airmar.mday,Airmar.Airmar_WDIR,'.');
legend('WDIR');
datetick('x');
ylabel('degrees from')
%xlim([tstart,tend]);
subplot(5,1,4);
plot(Airmar.mday,Airmar.Airmar_WSPD,'.');
legend('WSPD');
datetick('x');
ylabel('m/s');
%xlim([tstart,tend]);
subplot(5,1,5);
plot(Airmar.mday,Airmar.Airmar_magvar,'.');
hold on
plot(Airmar.mday,Airmar.Airmar_Hdg,'.');
legend('magvar','heading');
datetick('x');
ylabel('degrees');
%xlim([tstart,tend]);
% save .png
set(gcf,'position',[0 0 800 1200]);
saveas(gcf,[ DATADIR 'figs/' experiment '_' deployment '_Airmar_MET.png']);
fprintf(' * plot Airmar_MET.png\n');

% Position
figure(2);
sgtitle([ experiment ' Airmar Position'],'Interpreter','none');
subplot(3,1,1);
plot(Airmar.gps_lon,Airmar.gps_lat,'.');
subplot(3,1,2);
plot(Airmar.mday,Airmar.gps_height);
legend('gps height');
datetick('x');
%xlim([tstart,tend]);
subplot(3,1,3);
plot(Airmar.mday,Airmar.gps_sats,'.');
legend('gps sats');
datetick('x');
%xlim([tstart,tend]);
% save .png
set(gcf,'position',[0 0 800 1200]);
saveas(gcf,[ DATADIR 'figs/' experiment '_' deployment '_Airmar_pos.png']);
fprintf(' * plot Airmar_pos.png\n');

% Engineering
figure(3);
sgtitle([ experiment ' Airmar Engineering'],'Interpreter','none');
subplot(3,1,1);
plot(Airmar.mday,Airmar.Airmar_Pitch,'.');
hold on
plot(Airmar.mday,Airmar.Airmar_Roll,'.');
legend({'Pitch','Roll'});
datetick('x');
ylabel('degrees');
%xlim([tstart,tend]);
subplot(3,1,2);
plot(Airmar.mday,Airmar.Battery,'.');
legend('Battery');
ylabel('volts');
datetick('x');
%xlim([tstart,tend]);
subplot(3,1,3);
plot(Airmar.mday,Airmar.PanelTemp,'.');
legend('Campbell temperature');
datetick('x');
ylabel('degC');
%xlim([tstart,tend]);
% save .png
set(gcf,'position',[0 0 800 1200]);
saveas(gcf,[ DATADIR 'figs/' experiment '_' deployment '_Airmar_eng.png']);
fprintf(' * plot Airmar_eng.png\n');
