% plot thresholds for 1, 2, and 3 APs
% 3-dec-2014 (blh19)

clear; clc

ap1 = struct();
ap1.vm = load('vm1.txt');
ap1.istim = load('istim1.txt');

ap2 = struct();
ap2.vm = load('vm2.txt');
ap2.istim = load('istim2.txt');

ap3 = struct();
ap3.vm = load('vm3.txt');
ap3.istim = load('istim3.txt');

figure(1); clf; hold on
subplot(2,3,1); hold on
plot(ap1.istim(:,1),ap1.istim(:,2),'k')
xlabel('Time (ms)')
ylabel('I_{stim} (mA)')
title('I_{stim} = 1.83 mA')
axis([0 50 -3 1])
setfont(18)

subplot(2,3,4); hold on
plot(ap1.vm(:,1),ap1.vm(:,2),'k')
xlabel('Time (ms)')
ylabel('V_m (mV)')
setfont(18)

subplot(2,3,2); hold on
plot(ap2.istim(:,1),ap2.istim(:,2),'k')
xlabel('Time (ms)')
% ylabel('I_{stim} (mA)')
title('I_{stim} = 2.15 mA')
axis([0 50 -3 1])
setfont(18)

subplot(2,3,5); hold on
plot(ap2.vm(:,1),ap2.vm(:,2),'k')
xlabel('Time (ms)')
% ylabel('V_m (mV)')
setfont(18)

subplot(2,3,3); hold on
plot(ap3.istim(:,1),ap3.istim(:,2),'k')
xlabel('Time (ms)')
% ylabel('I_{stim} (mA)')
title('I_{stim} = 2.24 mA')
axis([0 50 -3 1])
setfont(18)

subplot(2,3,6); hold on
plot(ap3.vm(:,1),ap3.vm(:,2),'k')
xlabel('Time (ms)')
% ylabel('V_m (mV)')
setfont(18)

print -dpng fiberresponse_aps

