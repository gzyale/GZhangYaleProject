addpath(genpath('/glade/u/home/gang/gz_matlab/m_map'))

clear all
clc
clf

ntime = 169;
nlevs =   19;


%% processing data
ncid=netcdf.open('EFx_Ave_VerticalProfile.nc');
dataid = netcdf.inqVarID(ncid,'efxAve');
A = netcdf.getVar(ncid,dataid);
netcdf.close(ncid)
A(A<-10000.0) = NaN;
A(A>10000.0) = NaN;

%% plot the time-mean
A_ave = squeeze(nanmean(A,2));

plot(A_ave',2:20,'linewidth',3);
line([0,0],[0,20],'color','k')
xlabel('EFx','fontsize',14);
ylabel('Height (km)', 'fontsize',14);
print('EFx_VerticalProfile','-dpng');

