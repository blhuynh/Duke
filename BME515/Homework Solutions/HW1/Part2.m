% Author: N.Pelot
% BME 515 - Fall 2014
% HW1
% Part 2
% Intrinsically bursting neuron

format compact
clear
close all

fntsz = 14;

%% Model parameters
% Misc parameters
Nneurons = 1;

% a, b, c, d are [unitless]
a = 0.02;         % Time scale of the recovery variable u; smaller a = slower recovery
b = 0.2;          % Defines resting potential & sensitivity of u to subthreshold v
c = -55;          % Defines after-spike reset value of v
d = 4;            % Defines after-spike reset value of u

vrest = -70;    % [mV]

% Temporal parameters
tmax = 250;       % total simulation time [ms]
dt = 0.1;         % time step [ms]
t = dt:dt:tmax;   % time axis [ms]

% Synapse parameters
% 0.05 and 2.5 or 3
Pmax = 1;         % [unitless]
%Pmax_total = 100;%0.5;      % [unitless]
tau_s = 2;        % [ms]

% Applied current parameters
Istim = 20;        % stimulation current amplitude [units?]
del = 10;        % stimulation current start time [ms]
dur = 100;             % Duration of applied current Iin [ms]
I = zeros(numel(t),Nneurons);  % Applied current as a function of time [uA/cm^2]

%% Time marching
% Initial conditions
v = zeros(numel(t),Nneurons);
v(1,:) = vrest;
u = squeeze(b.*v(1,:));         % Only plotting v(t) so don't need u(t) time course

% Time loop
for t_ind = 1:numel(t)-1
   for neuron_ind = 1:Nneurons
      if ((t(t_ind) >= del(neuron_ind)) && (t(t_ind) <= del(neuron_ind) + dur(neuron_ind)))
         I(t_ind,neuron_ind) = Istim(neuron_ind);
      end

      % Euler's forward method
      v(t_ind+1,neuron_ind) = v(t_ind,neuron_ind) + dt*(0.04*v(t_ind,neuron_ind)^2 +...
         5*v(t_ind,neuron_ind) + 140 - u(neuron_ind) + I(t_ind,neuron_ind));

      if (v(t_ind+1,neuron_ind) >= 30) % Auxiliary after-spike resetting
         v(t_ind+1,neuron_ind) = c;
         u(neuron_ind) = u(neuron_ind) + d;
      else
         % Euler's forward method
         u(neuron_ind) = u(neuron_ind) + dt*(a*(b*v(t_ind+1,neuron_ind)-u(neuron_ind)));
      end
   end
end

%% Plotting
figure
subplot(2,1,1)
plot(t,v)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('V_m (mV)')
title('Intrinsically bursting neuron')

subplot(2,1,2)
plot(t,I)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('I_s_t_i_m')
title('Stimulus current')
ylim([-2*abs(Istim) 2*abs(Istim)])