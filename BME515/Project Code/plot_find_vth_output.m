% Plot data in find_vth_output.dat
% 1-Dec-2014 (blh19)

clear; clc

a = load('find_vth_output_Istim.txt');
tvec = a(:,1) * 1e-3;
stim = a(:,2);

figure(1); clf; hold on
plot(tvec,stim,'k')
xlabel('Time (s)')
ylabel('I_{stim} (mA)')
title('Stimulation Waveform')
setfont(18)

% print -dpng find_vth_output_Istim