% plot vm and istim from neuron
% 3-dec-2014 (blh19)

clear; clc
r = @(x,c) reshape(x,c,numel(x)/c)';

istim = fileload('istim.txt');
vm = fileload('vm.txt');

istim = r(istim,2);
vm = r(vm,2);

figure(1); clf; hold on
subplot(2,1,1); hold on
plot(istim(:,1),istim(:,2),'k')
xlabel('Time (ms)')
ylabel('I_{stim} (mA)')
title('Fiber Response')
setfont(18)

subplot(2,1,2); hold on
plot(vm(:,1),vm(:,2),'k')
xlabel('Time (ms)')
ylabel('V_m (mV)')
setfont(18)

print -dpng fiberresponse