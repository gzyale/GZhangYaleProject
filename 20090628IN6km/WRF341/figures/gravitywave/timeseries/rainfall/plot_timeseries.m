addpath(genpath('/glade/u/home/gang/gz_matlab/m_map'))

clear all
clc
clf

ntime = 169;
nlevs =   1;


%% processing data
fid=fopen('rainfall.bin');
A=fread(fid,inf,'single=>single'); % in mm
fclose(fid);

A(A<-10000.0) = NaN;
A(A>10000.0) = NaN;

R=zeros(1,ntime);
for k=2:ntime
   R(k)=A(k)-A(k-1); % mm/hr
end

figure(1)
plot(R,'linewidth',1.5);
grid on
set(gca,'fontsize',18,'xlim',[1 169]);
xlabel('Time','fontsize',18);
ylabel('Rainfall Rate (mm/h) ', 'fontsize',18);
title('WRF CTL: Rainfall', 'fontsize',18)
line([1 169],[0 0],'color','k');
saveas(gcf,'Rainfall_TimeSeries.png');

