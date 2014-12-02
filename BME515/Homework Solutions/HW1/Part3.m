% Nikki Pelot
% BME 515 - Fall 2014
% HW1
% Part 3

clear
format compact
close all

fntsz = 14;

%% *************************************************************************
% Synapse parameters
Pmax = 1;         % [unitless]
tau_s = 5;        % [ms]

dt = 0.01;           % Time step [ms]
tmax = 50;          % Simulation time [ms]
t = 0:dt:tmax;          % Time axis [ms]

%% *************************************************************************
% Alpha synapse function
Ps = Pmax.*t.*(exp(1-t./tau_s))./tau_s;
figure
plot(t,Ps)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('Ps(t) [unitless]')
title('Alpha function')

%% *************************************************************************
% Alpha synapse for multiple spike times
clear Ps
spike_times = [1 5 9 15];     % [ms]
for i = 1:numel(spike_times)
   % Ps(t) for each spike stored in columns
   Ps(:,i) = Pmax.*(t-spike_times(i)).*(exp(1-(t-spike_times(i))./tau_s))./tau_s;
   % Ps(t) is zero for t<spike_time
   Ps(1:floor(spike_times(i)/dt),i) = 0;
end
Ps_total = sum(Ps,2);

figure

subplot(2,1,1)
plot(t,Ps)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('Ps(t) [unitless]')
title('Individual time-shifted alpha functions')

subplot(2,1,2)
plot(t,Ps_total)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('Ps(t) [unitless]')
title('Total Ps(t)')