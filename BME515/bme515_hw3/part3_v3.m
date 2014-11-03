% BME 515 HW 3 Part 3 - Compound Nerve Action Potential (CNAP)
% Version 2 - Speed Optimization
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
dist = bsxfun(@times,stdev,randn(noaxons,1));
dist = bsxfun(@plus,dist,D);

%% create histogram
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

% print -dpng bme515_hw3_partA
close(1);

%% peak positive current trendline
% calculate trendline
for k = 1:length(fiberDlist)
    D = fiberDlist(k);
    fid = sprintf('part1_%gum.txt',D);
    Im(:,k) = fileload(fid) * (nodelength*pi*D) * 1e-8; % (mA/cm^2) * (um)^2 -> (mA)
end
PIm = polyfit(fiberDlist,max(Im),1);

%% conduction velocity trendline
data = fileload('part2.txt');
data = reshape(data,2,numel(data)/2)';
Pcv = polyfit(data(:,1),data(:,2),1);

%% create bundle of virtual axons
nonodes = 151;

% anonymous function for absolute node locations
nodelocations = @(INL) 0:INL+nodelength:(INL+nodelength)*(nonodes-1);

% reshape distribution of fiber diameters
dist = reshape(dist,numel(dist),1);

% node locations of all axons
% output (array) is not the same size as input (vector)
% output size should be noaxons x nonodes
% locs = cell2mat(arrayfun(nodelocations,dist*100,'uni',0));

%% electrode parameters
k = 10*1e3; % inter-electrode spacing (mm -> um)
Re = 500; % extracellular resistivity (ohm-cm)
d = 1*1e3; % perpendicular electrode-fiber distance (mm->um)
xrec = 15e3; % (mm)
noelecs = 3;
nodist = 2;

%% time parameters
tdel = 1;
pw = 1;
tstop = 20; % (ms)
dt = 0.02; % (ms)
tvec = 0:dt:tstop-dt;

%% create AP Im profile
% summed voltage for each tripole contact for all time (no_tpts x noelecs)
Vsum = zeros(numel(tvec),noelecs);

% select axon
for axno=1:noaxons*2
    
    % axon-dependent parameters
    D     = dist(axno);
    INL   = D*100;
    Axon  = nodelocations(INL);
    Laxon = max(Axon); % (um)
    
    % find AP profile of Im of closest fiber diameter
    [y,i] = min(abs(fiberDlist-D));
    closestD = fiberDlist(i);
    fid = sprintf('part1_%gum.txt',closestD);
    
    % scale Im using peak positive current trendline
    Im = Im/max(Im) * polyval(PIm,D);
    
    % determine conduction velocity using trendline
    CV = polyval(Pcv,D) * 1e3; % (m/s) -> (um/ms)
    
    % isolate AP part of Im
    % AP = (Im > 5e-10) + (Im < -1e-10);
    AP = (Im > 0.005*abs(max(Im))) + (Im < -0.001*abs(max(Im)));
    AP(tdel/dt:(tdel+pw)/dt) = 0;
    Im_sub = Im(AP>0)';
    
    % calculate distance from all axon nodes to the three recording contacts
    node2v = @(xn,xv) sqrt( (xn-xv).^2 + d^2);
    xlocs = [xrec-k xrec xrec+k]; % absolute x-dist (um)
    r = bsxfun(node2v,repmat(Axon',[1 3]),xlocs); % size of 152x3 (nonodes x noelecs)
    
    % calculate times when AP will hit each NoR
    n2ntime = (INL+nodelength)/CV; % node to node time (ms)
    fprintf(sprintf('D=%1.3f um \t INL=%.4g um \t CV=%4.1f um/ms \t N2N=%g ms\n',...
        D,INL,CV,n2ntime))
    nodetimes = 0:n2ntime:(nonodes-1)*n2ntime;
    nodetimes_ind = floor(nodetimes/dt); % time index at which AP vector starts
    nodetimes_ind(1) = 1; % set AP to start at node 1 at time index 1
    
    % insert Im for all nodes at all time points for a single axon
    Im_allnodes = zeros(numel(tvec),nonodes); % (no_tpts x nonodes)
    for a=1:length(nodetimes_ind)
        Im_allnodes(nodetimes_ind(a):nodetimes_ind(a)+numel(Im_sub)-1,a) = Im_sub';
    end
    
    % calculate phi for all nodes at all time points for a single axon
    % phi should be (noelecs x nonodes x no_tpts)
    potential = @(Im,rdist) Im*Re ./ (1e4*4*pi.*rdist);
    for a=1:noelecs
        thisDist = r(:,a);
        phi{a,1} = bsxfun(potential,Im_allnodes,thisDist')';
    end
    
    % sum phi over all nodes for all time points for a single axon
    V = cell2mat(cellfun(@sum,phi,'uni',0))';
    Vsum = Vsum + V;
end

% calculate net voltage
Vnet = Vsum(:,2) - (Vsum(:,1)+Vsum(:,3))/2;

%% plot phi lines
figure(2); clf; hold on
subplot(2,1,1); hold on
plot(tvec,Vsum(:,1),'k')
plot(tvec,Vsum(:,2),'b')
plot(tvec,Vsum(:,3),'g')
legend('V1','V2','V3')
title(sprintf('Tripole Recording Electrode at %gmm',xrec*1e-3))
xlabel('Time (ms)'); ylabel('Potential (V)'); setfont(18)

subplot(2,1,2); hold on
plot(tvec,Vnet,'k')
title('Net Recorded Voltage (ECAP)')
xlabel('Time (ms)'); ylabel('V_m(t) (V)'); setfont(18)


%% timing
toc
 No newline at end of file
