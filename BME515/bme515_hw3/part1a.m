% BME 515 HW 3 Part 1
% 30-Oct-2014 (blh19)

clear; clc
fiberDlist = 2:2:20;
workdir = pwd();

tstop = 15;
dt = 0.02;
tvec = 0:dt:tstop-dt;

nodelength = 1; % (um)
D = 20;

node1 = 6;
node2 = 30;

cd('Data')
d1 = fileload(sprintf('part1a_node%g.txt',node1)) * (nodelength*pi*D) * 1e-8; % (mA/cm^2) * (um)^2 -> (mA)
d2 = fileload(sprintf('part1a_node%g.txt',node2)) * (nodelength*pi*D) * 1e-8; % (mA/cm^2) * (um)^2 -> (mA)
cd(workdir)

figure(1); clf; hold on
subplot(2,1,1); hold on
plot(tvec,d1)
title(sprintf('Node %g',node1))
xlabel('Time (ms)')
ylabel('Current (mA)')
setfont(18)
subplot(2,1,2); hold on
plot(tvec,d2)
title(sprintf('Node %g',node2))
xlabel('Time (ms)')
ylabel('Current (mA)')
setfont(18)

% print -dpng bme515_hw3_part1a