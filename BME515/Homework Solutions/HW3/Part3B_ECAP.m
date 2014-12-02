tic
clear all
%close all
format compact

fntsz = 14;
lnwdth = 2;
scrnsz = get(0,'ScreenSize');
%% Parameters
% Prepare fiber diameters in nerve bundle (um)
load Dvalues.mat
%Dvalues = Dvalues(1:100);
%Dvalues = Dvalues(101:200);

% Electrical parameters
rhoe = 500;    % resistivity of external conducting medium (ohm cm)

% Load Im(t) for D=2um (nA)
D = 2;
fname = ['Im_' num2str(D) 'um.dat'];
tmp = importdata(fname);
t = tmp.data(:,1);
Im_data = tmp.data(:,2);
% Divide by 2 to get Im(t) for D=1um axon
Im_data = Im_data/2;
% Only use Im_data from t=15 to 22ms
Im_data = Im_data(t>6 & t<14);

% Downsample to dt=20us instead of 1us
%Im_data = Im_data(1:20:end);

% Temporal parameters (ms)
tstart = 0;
tend = 75;
dt = 0.02; % same dt as for Im_data
t = tstart:dt:tend;

% Axon length (um)
axon_length = 200000;

% Spacing between each pair of electrodes; "k" in the homework instructions (um)
electrode_spacing = 100;
% Distance from start of axon to each electrode in tripolar recording (um)
x_elec1 = 10000; 
x_elec2 = x_elec1+electrode_spacing;         % placed over center node; defines half of axon length
x_elec3 = x_elec2+electrode_spacing*2;
% Perpendicular distance from electrode to axon (um)
z_elec = 1000;


%% Simulation
% Prepare arrays for recording; each value is already summed across nodes
volrec1 = zeros(numel(Dvalues),numel(t));
volrec2 = zeros(numel(Dvalues),numel(t));
volrec3 = zeros(numel(Dvalues),numel(t));

for i = 1:numel(Dvalues)
   D = Dvalues(i);
   INL = 100*D;                           % internodal distance (um)
   l = 1;                                 % length of node (um)
   nnodes = floor(axon_length/(INL+l));   % number of nodes in axon
   
   % x coordinates of nodes
   x_nodes = 0:(INL+l):axon_length-(INL+l);

   % Distance from each node to each recording electrode
   r_elec1 = sqrt((x_elec1-x_nodes).^2 + (z_elec)^2);
   r_elec2 = sqrt((x_elec2-x_nodes).^2 + (z_elec)^2);
   r_elec3 = sqrt((x_elec3-x_nodes).^2 + (z_elec)^2);
   
   % Prepare array for keeping track of transmembrane currents
   imem = zeros(numel(t),nnodes);

   % Conduction velocity is 800um/ms per um of fiber diameter
   CV = 800*D;
   
   % Delay of start times for current waveform between adjacent nodes (ms)
   imem_dt = (INL+l)/CV;

   % ******************************************
   % *** METHOD 1: Convolution (slower)
%    % Compute current time courses for each node (time shifted)
%    for ii = 1:nnodes
%       % Index of start time of current time course
%       % Start time for node N = ii*imem_dt
%       % Convert to index by dividing by dt
%       imem_start = zeros(numel(t),1);
%       imem_start(floor(ii*imem_dt/dt)) = 1;
%       
%       % Convolve Im(t) with delta at start time
%       imem_tmp = conv(imem_start,Im_data.*D);
%       
%       % Save into Im array (# time points x # nodes)
%       imem(:,ii) = imem_tmp(1:numel(t),1);
%    end
%    
   % ******************************************
   
   % ******************************************
   % *** METHOD 2: Indexing (faster)
   % Compute current time courses for each node (time shifted) and the
   % distance between each node and each recording electrode
   for ii = 1:nnodes
      % Index of start time of current time course
      imem_start = floor(ii*imem_dt/dt);
      % Only include current time courses that fall within simulation time
      if (imem_start + numel(Im_data) - 1)*dt < t(end)
         imem_end = imem_start + numel(Im_data) - 1;
      else
         imem_end = numel(t);
      end
      % Populate current time course for the node; use entire time course
      % (Current Values), unless we're running over the simulation time
      imem((imem_start:imem_end-1),ii) = Im_data(1:numel(imem_start:imem_end-1)).*D;
   end
   % ******************************************
   
   % Compute the evoked potential with a tripolar electrode for each time
   % point (summing over all nodes) (nV)
   volrec1_tmp = 10^4 .* rhoe .* imem ./ (4 * pi .* (ones(numel(t),1)*r_elec1));
   volrec2_tmp = 10^4 .* rhoe .* imem ./ (4 * pi .* (ones(numel(t),1)*r_elec2));
   volrec3_tmp = 10^4 .* rhoe .* imem ./ (4 * pi .* (ones(numel(t),1)*r_elec3));
   volrec1(i,:) = sum(volrec1_tmp,2);
   volrec2(i,:) = sum(volrec2_tmp,2);
   volrec3(i,:) = sum(volrec3_tmp,2);
end

volrec = volrec2 - (volrec1 + volrec3)./2;

%% Plot results
figure
set(gca,'fontsize',fntsz)
plot(t,sum(volrec(:,:),1),'linewidth',lnwdth)
title('Recorded ECAP')
xlabel('Time (ms)')
ylabel('V (nV)')

toc