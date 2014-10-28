% BME 515 HW 3 Part 3 - Compound Nerve Action Potential (CNAP)
% 30-Oct-2014 (blh19)

clear; clc; format short g
addpath Data

%% define parameters for a bundle of virtual axons
D = [10 6];
stdev = [2 2];
s = rng(0); % random seed generator
noaxons = 100;
fiberDlist = 2:2:20;
nodelength = 1; % (um)

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

% print -dpng bme515_hw3_part3A
close(1)

%% create bundle of virtual axons
nonodes = 151;
axons = cell(1,2);

for a = 1:2
    for b = 1:noaxons
        thisINL = 100*dist(b,a);
        axons{a}(b,:) = -round(noaxons/2)*thisINL:thisINL:round(noaxons/2)*thisINL;
    end
end

%% peak positive current trendline
% calculate trendline
for k = 1:length(fiberDlist)
    D = fiberDlist(k);
    fid = sprintf('part1_%gum.txt',D);
    Im(:,k) = fileload(fid) * (nodelength*pi*D) * 1e-8; % (mA/cm^2) * (um)^2 -> (mA)
end
P = polyfit(fiberDlist,max(Im),1);

%% conduction velocity trendline
data = fileload('part2.txt');
data = reshape(data,2,numel(data)/2)';
Pcv = polyfit(data(:,1),data(:,2),1);

%% action potential
k = 10; % interelectrode spacings
Re = 500; % extracellular resistivity (ohm-cm)
d = 1; % perpendicular electrode-fiber distance (mm)

subpop = 1;
axno   = 5;
thisAxon = axons{subpop}(axno,:);
thisD = dist(axno,subpop);
[y,i]=min(abs(fiberDlist-thisD));
closestD = fiberDlist(i);
fid = sprintf('part1_%gum.txt',closestD);
Im = fileload(fid)*(nodelength*pi*closestD) * 1e-8;

% scale Im using peak positive current trendline
Im = Im/max(Im) * polyval(P,thisD);

% determine conduction velocity using trendline
thisCV = polyval(Pcv,thisD) * 1e3; % (m/s) -> (um/ms)
thisINL = 100*thisD;
Laxon = nonodes*(thisINL+nodelength); % (um)

% time parameters
tdel = 1;
pw = 1;
tstop = 25; % (ms)
dt = 0.02; % (ms)

% isolate AP part of Im
% AP = (Im > 5e-10) + (Im < -1e-10);
AP = (Im > 0.005*abs(max(Im))) + (Im < -0.001*abs(max(Im)));
AP(tdel/dt:(tdel+pw)/dt) = 0;
Im_sub = Im(AP>0)';

% iterate over time loop
tvec = 0:dt:tstop-dt;
dt = thisINL/thisCV;

% set up recording electrode
z_rec = 1000; % (um)
x_rec_ind = floor(nonodes/2);
xaxon(1) = nodelength/2;
for k=2:nonodes
    xaxon(k) = xaxon(k-1)+thisINL+nodelength;
end
xstim=xaxon(x_rec_ind);
k = 25; % inter-electrode spacing (mm)
r(:,2) = sqrt( (xaxon-xstim).^2 + z_rec^2);
r(:,1) = sqrt( (xaxon-xstim+k*1e3).^2 + z_rec^2);
r(:,3) = sqrt( (xaxon-xstim-k*1e3).^2 + z_rec^2);
noelecs = 3;

% ECAP vector
V = zeros(floor(tstop/dt),noelecs);

for a = 1:tstop/dt
    t = a*dt;
    
    % calculate the Im for all nodes at a single time point
    Im_allnodes = zeros(1,nonodes);
    if a<=nonodes
        Im_allnodes(1:a) = fliplr(Im_sub(1:a));
    elseif a>nonodes && a<numel(Im_sub)
        Im_allnodes = fliplr(Im_sub(a-nonodes+1:a));
    elseif a>numel(Im_sub)
        Im_allnodes = fliplr(Im_sub(a-nonodes:end));
        Im_allnodes = [zeros(1,nonodes-numel(Im_allnodes)) Im_allnodes];
    end
    
    if numel(Im_allnodes)~=nonodes
        error('Number of nodes incorrect.')
    end
    
    for b = 1:nonodes
        for c = 1:noelecs
            V(a,c) = V(a,c) + 4*pi*Re*r(b,c)*Im_allnodes(b)/1e4;
        end
    end
end

%% plot electrode potentials
tvec=0:dt:tstop-dt;
figure(2); clf; hold on
subplot(2,1,1); hold on
plot(tvec,V(:,1),'k')
plot(tvec,V(:,2),'b')
plot(tvec,V(:,3),'g')
legend('V1','V2','V3')
xlabel('Time (ms)'); ylabel('Potential'); setfont(18)

subplot(2,1,2); hold on
plot(tvec,V(:,2)-(V(:,1)+V(:,3))/2,'k')
xlabel('Time (ms)'); ylabel('V_m(t)'); setfont(18)