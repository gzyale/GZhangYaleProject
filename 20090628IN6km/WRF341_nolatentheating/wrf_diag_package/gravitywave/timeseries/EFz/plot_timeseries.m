addpath(genpath('/glade/u/home/gang/gz_matlab/m_map'))

clear all
clc
clf

ntime = 169;
nlevs =   19;


%% processing data
ncid=netcdf.open('EFz_Ave_VerticalProfile.nc');
dataid = netcdf.inqVarID(ncid,'efzAve');
A = netcdf.getVar(ncid,dataid);
netcdf.close(ncid)
A(A<-10000.0) = NaN;
A(A>10000.0) = NaN;

ncid=netcdf.open('/glade/p/uyle0010/20090628IN6km/WRF341/figures/gravitywave/EPrelation/16km/EFz_raw_d01_2009-06-28_00.nc');
    xlongid = netcdf.inqVarID(ncid,'XLONG');
    xlong = netcdf.getVar(ncid,xlongid);
    lon = xlong(:,1);
    xlatid = netcdf.inqVarID(ncid,'XLAT');
    xlat = netcdf.getVar(ncid,xlatid);
    lat = xlat(1,:);
netcdf.close(ncid);

figure(1)
plot(squeeze(A(15,:)),'linewidth',1.5);
grid on
set(gca,'fontsize',18,'xlim',[1 169]);
xlabel('Time','fontsize',18);
ylabel('EFz (W/m2) ', 'fontsize',18);
title('WRF NO LH: EFz at 16 km', 'fontsize',18)
line([1 169],[0 0],'color','k');
saveas(gcf,'EFz_TimeSeries_16km.png');

figure(2)
plot(squeeze(A(2,:)),'linewidth',1.5);
grid on
set(gca,'fontsize',18,'xlim',[1 169]);
xlabel('Time','fontsize',18);
ylabel('EFz (W/m2) ', 'fontsize',18);
title('WRF NO LH: EFz at 3 km', 'fontsize',18)
line([1 169],[0 0],'color','k');
saveas(gcf,'EFz_TimeSeries_03km.png');
