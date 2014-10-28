% BME 515 HW 3 Part 4 - Electroneurogram (ENG)
% 30-Oct-2014 (blh19)
clear; clc; format short g

%% 4a. Histogram of Firing Rates
% uniform distribution between 10-20 Hz

% set seed
s = rng(0);

% generate uniform distribution
noaxons = 1e2;
dist = 10+10*rand(noaxons,1);

% plot histogram
figure(1); clf; hold on
hist(dist)
title('Spontaneous Basal Firing Rates')
xlabel('Firing Frequency (Hz)')
ylabel('Occurrences')
setfont(18)

print -dpng bme515_hw3_part4a