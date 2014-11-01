% BME 515 HW 3 Part 3 - Compound Nerve Action Potential (CNAP)
% 30-Oct-2014 (blh19)

clear; clc; format short g; tic
addpath Data

%% define parameters for a bundle of virtual axons
D = [6.5 0.80];
stdev = [0.47 0.09];
s = rng(0); % random seed generator
noaxons = 50;
fiberDlist = 2:2:20;
nodelength = 1; % (um)

% distribution of fiber diameters
dist(:,1) = stdev(1) .* randn(noaxons,1) + D(1);
dist(:,2) = stdev(2) .* randn(noaxons,1) + D(2);

% create histogram
figure(1); clf; hold on
subplot(2,1,1); hold on
hist(dist(:,1))
title(sprintf('%g axons with a fiber diameter of %gum (mean=%.2f, std=%.2f)',...
    noaxons,D(1),mean(dist(:,1)),std(dist(:,1))))
setfont(18)

subplot(2,1,2); hold on
title(sprintf('%g axons with a fiber diameter of %gum (mean=%.2f, std=%.2f)',...
    noaxons,D(2),mean(dist(:,2)),std(dist(:,2))))
hist(dist(:,2))
setfont(18)

% print -dpng bme515_hw3_part3A
close(1)
% return

%% create bundle of virtual axons
nonodes = 151;
axons = cell(1,2);

for a = 1:2 % subpop
    for b = 1:noaxons
        thisINL = 100*dist(b,a);
        axons{a}(b,:) = -round(noaxons/2)*thisINL:thisINL:round(noaxons/2)*thisINL;
    end
end

% reshape distance and axon node location matrices
dist = [dist(:,1); dist(:,2)];
axons = [axons{1}; axons{2}];

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

% time parameters
tdel = 1;
pw = 1;
tstop = 20; % (ms)
dt = 0.02; % (ms)
tvec = 0:dt:tstop-dt;

% set up recording electrode
z_rec = 1000; % (um)
xrec = 5000; % mm
noelecs = 3;

for axno=1%:noaxons*2 % # distributions = 2
    thisAxon = axons(axno,:);
    thisD = dist(axno);
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
    
    % isolate AP part of Im
    % AP = (Im > 5e-10) + (Im < -1e-10);
    AP = (Im > 0.005*abs(max(Im))) + (Im < -0.001*abs(max(Im)));
    AP(tdel/dt:(tdel+pw)/dt) = 0;
    Im_sub = Im(AP>0)';
    
    % set up recording electrode
    xaxon(1) = nodelength/2;
    for k=2:nonodes
        xaxon(k) = xaxon(k-1)+thisINL+nodelength;
    end
    k = 10; % inter-electrode spacing (mm)
    r(:,2) = sqrt( (xaxon-xrec).^2 + z_rec^2);
    r(:,1) = sqrt( (xaxon-xrec+k*1e3).^2 + z_rec^2);
    r(:,3) = sqrt( (xaxon-xrec-k*1e3).^2 + z_rec^2);
    
    % ECAP vector
    V = zeros(floor(tstop/dt),noelecs);
    
    % calculate the Im for all nodes at all time points
    Im_allnodes = zeros(numel(tvec),nonodes);
    
    % calculate times when AP will hit each NoR
    nodetimes = zeros(1,nonodes);
    for k=2:numel(nodetimes)
        nodetimes(k)=nodetimes(k-1)+thisINL/thisCV;
    end
    
    nodecount = 1;
    for a = 1:tstop/dt
        t = a*dt;
        
        if nodecount > nonodes
            break
        end
        
        if t > nodetimes(nodecount)
            Im_allnodes(a:a+numel(Im_sub)-1,nodecount,axno)=Im_sub';
            nodecount = nodecount+1;
        end
    end
    fprintf(sprintf('Axon %g.\n',axno))
end

fprintf('Calculating phi.\n')
for a=1:numel(tvec)
    for b=1:nonodes
        for c=1:noelecs
            phi(a,b,c,:) = 4*pi*Re*r(b,c)*Im_allnodes(a,b,:)/1e4;
        end
    end
    if mod(a,10)==0
        fprintf(sprintf('t=%gms.\n',a*dt))
    end
end

V = squeeze(sum(sum(phi,2),4));

%% plot electrode potentials
tvec=0:dt:tstop-dt;
figure(2); clf; hold on
subplot(2,1,1); hold on
plot(tvec,V(:,1),'k')
plot(tvec,V(:,2),'b')
plot(tvec,V(:,3),'g')
legend('V1','V2','V3')
title(sprintf('Tripole Recording Electrode at %gmm',xrec*1e-3))
xlabel('Time (ms)'); ylabel('Potential (V)'); setfont(18)

subplot(2,1,2); hold on
plot(tvec,V(:,2)-(V(:,1)+V(:,3))/2,'k')
xlabel('Time (ms)'); ylabel('V_m(t) (V)'); setfont(18)