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

%% plot the time-mean
A_ave = squeeze(nanmean(A,2));

figure(1)
plot(A_ave',2:20,'linewidth',1.5);
grid on
set(gca,'fontsize',24);
xlabel('EFz','fontsize',18);
ylabel('Height (km)', 'fontsize',18);
title('WRF CTL: EFz ', 'fontsize',18);
line([0,0],[0,20],'color','k');
saveas(gcf,'EFz_VerticalProfile.png');

