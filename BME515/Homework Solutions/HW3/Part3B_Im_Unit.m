format compact
clear
close all

fntsz = 14;
lnwdth = 3;
scrnsz = get(0,'ScreenSize');

%% Data parameters
D_values = 2:2:20;         % Fiber diameters [um]


%% Prep "unit" Im(t) time course
% Load Im(t) for D=2um
D = 2;
fname = ['Im_' num2str(D) 'um.dat'];
tmp = importdata(fname);
t = tmp.data(:,1);
Im_data = tmp.data(:,2);

% Divide by 2 to get Im(t) for D=1um axon
Im_data = Im_data/2;

% Plot time course
figure
set(gca,'fontsize',fntsz)
plot(t,Im_data,'linewidth',lnwdth)
xlabel('Time (ms)')
ylabel('Im (nA)')
title('Im(t)/2 for D=2um')
xlim([6 14])
