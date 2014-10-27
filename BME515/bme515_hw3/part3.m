% BME 515 HW 3 Part 3 - Compound Nerve Action Potential (CNAP)
% 30-Oct-2014 (blh19)

clear; clc

%% define parameters for a bundle of virtual axons
D = [10 6];
stdev = [2 2];
s = rng(0); % random seed generator
noaxons = 100;

% distribution of fiber diameters
dist(:,1) = stdev(1) .* randn(noaxons,1) + D(1);
dist(:,2) = stdev(2) .* randn(noaxons,1) + D(2);

% create histogram
figure(1); clf; hold on
subplot(2,1,1); hold on
hist(dist(:,1))
title(sprintf('%g axons with a fiber diameter of %gum (mean=%.1f, std=%.1f)',...
    noaxons,D(1),mean(dist(:,1)),std(dist(:,1))))
setfont(18)

subplot(2,1,2); hold on
title(sprintf('%g axons with a fiber diameter of %gum (mean=%.1f, std=%.1f)',...
    noaxons,D(2),mean(dist(:,2)),std(dist(:,2))))
hist(dist(:,2))
setfont(18)

%% create bundle of virtual axons
nonodes = 21;
axons = cell(1,2);

for a = 1:2
    for b = 1:noaxons
        thisINL = 100*dist(b,a);
        axons{a}(b,:) = -round(noaxons/2)*thisINL:thisINL:round(noaxons/2)*thisINL;
    end
end

%% action potential
k = 10; % interelectrode spacings
Re = 500; % extracellular resistivity (ohm-cm)
d = 1; % perpendicular electrode-fiber distance (mm)

