% Nikki Pelot
% BME 515 - Fall 2014
% HW1
% Part 4

format compact
clear
close all

fntsz = 14;
scrnsz = get(0,'screensize'); % left bottom width height

%% Model parameters
% Misc parameters
Nneurons = 2;

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
Istim = [30 30];        % stimulation current amplitude [units?]
del = [0 10];        % stimulation current start time [ms]
dur = [tmax tmax];             % Duration of applied current Iin [ms]
I = zeros(numel(t),Nneurons);  % Applied current as a function of time [uA/cm^2]


% AP parameters
Nspikes = zeros(Nneurons,1);            % Count the number of spikes (number of times Vthreshold crossed with a rising edge)
max_Nspikes = 100;            % Just for the sake of pre-allocating arrays
spike_times = zeros(max_Nspikes,Nneurons);   % Keep track of the spike times for each neuron
Ps = zeros(numel(t),max_Nspikes,Nneurons);

%% Time marching
% Initial conditions
v = zeros(numel(t),Nneurons);
v(1,:) = vrest;
u = squeeze(b.*v(1,:));         % Only plotting v(t) so don't need u(t) time course

% Time loop
for t_ind = 1:numel(t)-1
   % Ps for neuron 1 uses spike times from neuron 2
   for spike_ind = 1:Nspikes(2)
      Ps(t_ind,spike_ind,1) = (Pmax.*(t(t_ind)-spike_times(spike_ind,2)).*...
         (exp(1-(t(t_ind)-spike_times(spike_ind,2))./tau_s))./tau_s);
      % Ps(t) is zero for t<spike_time
      Ps(1:floor(spike_times(spike_ind,2)/dt),spike_ind,1) = 0;
   end
   
   % Ps for neuron 2 uses spike times from neuron 1
   for spike_ind = 1:Nspikes(1)
      Ps(t_ind,spike_ind,2) = (Pmax.*(t(t_ind)-spike_times(spike_ind,1)).*...
         (exp(1-(t(t_ind)-spike_times(spike_ind,1))./tau_s))./tau_s);
      % Ps(t) is zero for t<spike_time
      Ps(1:floor(spike_times(spike_ind,1)/dt),spike_ind,2) = 0;
   end
   
     Ps_total = squeeze(sum(Ps,2));
%        tmp = find(Ps_total>Pmax_total);
%    Ps_total(tmp) = Pmax_total;
   
   for neuron_ind = 1:Nneurons
      if ((t(t_ind) >= del(neuron_ind)) && (t(t_ind) <= del(neuron_ind) + dur(neuron_ind)))
         I(t_ind,neuron_ind) = Istim(neuron_ind);
      end

      % Euler's forward method
      v(t_ind+1,neuron_ind) = v(t_ind,neuron_ind) + dt*(0.04*v(t_ind,neuron_ind)^2 +...
         5*v(t_ind,neuron_ind) + 140 - u(neuron_ind) + I(t_ind,neuron_ind)) - ...
         Ps_total(t_ind,neuron_ind);

      if (v(t_ind+1,neuron_ind) >= 30) % Auxiliary after-spike resetting
         v(t_ind+1,neuron_ind) = c;
         u(neuron_ind) = u(neuron_ind) + d;
         
         % Record the spike
         Nspikes(neuron_ind) = Nspikes(neuron_ind) + 1; %  Count the number of spikes
         spike_times(Nspikes(neuron_ind),neuron_ind) = t(t_ind);  % Keep track of spike times
         
      else
         % Euler's forward method
         u(neuron_ind) = u(neuron_ind) + dt*(a*(b*v(t_ind+1,neuron_ind)-u(neuron_ind)));
      end
   end
end

%% Plotting
figure
set(gcf,'OuterPosition',[scrnsz(1)/3 scrnsz(2) scrnsz(3)/3 scrnsz(4)])

subplot(3,1,1)
plot(t,v)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('V_m (mV)')
title('Two oscillating intrinsically bursting neurons')

subplot(3,1,2)
plot(t,I,'linewidth',3)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('I_s_t_i_m')
title('Stimulus current')
ylim([-2*abs(Istim(1)) 2*abs(Istim(1))])
legend('Neuron 1','Neuron 2','location','southeast')

subplot(3,1,3)
plot(t,Ps_total)
set(gca,'fontsize',fntsz)
xlabel('Time (ms)')
ylabel('Ps(t) [unitless]')
title('Total Ps(t)')