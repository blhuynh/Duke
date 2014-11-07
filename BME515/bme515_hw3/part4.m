% BME 515 HW 3 Part 4 - Electroneurogram (ENG)
% 30-Oct-2014 (blh19)
clear; clc; format short g
addpath Data

%% 4a. Histogram of Firing Rates
% uniform distribution between 10-20 Hz

% set seed
s = rng(0);

% generate uniform distribution
noaxons = 1e2;
fdist = 10+(10*rand(noaxons,1));

% plot histogram
if false
    figure; clf; hold on
    hist(dist)
    title('Spontaneous Basal Firing Rates')
    xlabel('Firing Frequency (Hz)')
    ylabel('Occurrences')
    setfont(18)
end
% print -dpng bme515_hw3_part4a

%% ------------------------- 4b. Electroneurogram -------------------------

%% define parameters for a bundle of virtual axons
% D = [6.5 0.80];
% stdev = [0.47 0.09];
D = [6.5 19.8];
stdev = [0.47 std([18.7 19.3 20.1 18.4 16.2 20.6 21.0 20.8 22.5 20.6])];
s = rng(0); % random seed generator
noaxons = 50;
fiberDlist = 2:2:20;
nodelength = 1; % (um)

% distribution of fiber diameters
dist = bsxfun(@times,stdev,randn(noaxons,1));
dist = bsxfun(@plus,dist,D);

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
tstop = 500; % (ms)
dt = 0.02; % (ms)
tvec = 0:dt:tstop-dt;

%% create AP Im profile
% summed voltage for each tripole contact for all time (no_tpts x noelecs)
Vsum = zeros(numel(tvec),noelecs);

% select axon
axonpop = [1:noaxons; noaxons+1:2*noaxons];
for var1 = 1:2
    theseAxons = axonpop(var1,:);
    for axno=theseAxons
        
        % axon-dependent parameters
        D     = dist(axno);
        INL   = D*100+nodelength;
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
        
        % spontaneous basal firing rate frequency
        f = fdist(axno);
        firetimes = 0:1/f*1e3:tstop;
        firetimes_ind = floor(firetimes/dt);
        
        % calculate times when AP will hit each NoR
        n2ntime = INL/CV; % node to node time (ms)
        %         fprintf(sprintf('D=%1.3f um \t INL=%.4g um \t CV=%4.1f um/ms \t N2N=%g ms\n',...
        %             D,INL,CV,n2ntime))
        nodetimes = 0:n2ntime:(nonodes-1)*n2ntime;
        nodetimes_ind = floor(nodetimes/dt); % time index at which AP vector starts
        nodetimes_ind(1) = 1; % set AP to start at node 1 at time index 1
        
        % insert Im for all nodes at all time points for a single axon
        Im_allnodes = zeros(numel(tvec),nonodes); % (no_tpts x nonodes)
        for b=1:length(firetimes_ind)
            for a=1:length(nodetimes_ind)
                startIndex = firetimes_ind(b)+nodetimes_ind(a);
                stopIndex  = firetimes_ind(b)+nodetimes_ind(a)+numel(Im_sub)-1;
                
                % if fires close to end of tstop, only save up to tstop
                if stopIndex <= tstop/dt
                    Im_allnodes(startIndex:stopIndex,a) = Im_sub';
                elseif stopIndex > tstop/dt
                    Im_allnodes(startIndex:end,a) = Im_sub(1:end-(stopIndex-tstop/dt))';
                end
                
                %             figure(2); clf; hold on
                %             plot(tvec,Im_allnodes(1:tstop/dt,1),'k')
            end
        end
        
        % calculate phi for all nodes at all time points for a single axon
        % phi should be (noelecs x nonodes x no_tpts)
        %     potential = @(Im,rdist) 4*pi*Re.*rdist.*Im/1e4;
        potential = @(Im,rdist) Im*Re*1e4 ./ (4*pi.*rdist);
        for a=1:noelecs
            thisDist = r(:,a);
            phi{a,1} = bsxfun(potential,Im_allnodes,thisDist')';
        end
        
        % sum phi over all nodes for all time points for a single axon
        V = cell2mat(cellfun(@sum,phi,'uni',0))';
        Vsum = Vsum + V(1:tstop/dt,:);
    end
    % calculate net voltage
    Vnet(var1,:) = Vsum(:,2) - (Vsum(:,1)+Vsum(:,3))/2;
end

%% debugging
figure(1); clf; hold on
subplot(2,1,1); hold on
plot(tvec,Vnet(1,:),'k')
title('Electroneurogram (Population 1)')
xlabel('Time (ms)'); ylabel({'Net Recorded','Voltage (V)'}); setfont(18)

subplot(2,1,2); hold on
plot(tvec,Vnet(2,:),'k')

title('Electroneurogram (Population 2)')
xlabel('Time (ms)'); ylabel({'Net Recorded','Voltage (V)'}); setfont(18)
% print('-dpng','bme515_hw3_part4_2')