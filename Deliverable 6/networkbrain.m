% General Network of Neurons
% Written by Bryan Howell (9/23/12)
% Current capabilities
% 1. Can specify the electrical properties of each neuron
% 2. Can specify both the rise and decay time constants of every synapse
% 3. Can specify the conductance of every synapse
% 4. Can stimulate all three neurons
% 5. Includes a relative refractory period option

clear all;
clc;
addpath('/Users/blhuynh/Documents/F2014/BME503/Project 1/Network Model')

global nrn_param;
global syn_param;
global stim_param;
nrn_param  = struct();
syn_param  = struct();
stim_param = struct();

% number of neurons
nonrn = 2;
nrn_param.nn = nonrn;

% ---------------------- Stimulus Parameters ------------------------------

% times in ms
to = 0;    % initial time (ms)
tf = 50;  % final time (ms)
dt = 0.05; % initial time step (ms)
t  = to:dt:tf;
no_tpts = length(t);

stim_param.dur   = zeros(nonrn,1); % duration of stimulus pulse (ms)
stim_param.delay = zeros(nonrn,1);  % delay before stimulating (ms)
stim_param.mag   = zeros(nonrn,1); % magnitude of stimulus pusle (nA)

% stimulate neuron 1
% stim_param.dur   = [50;0;0;0]; %This needs to be same size as 
% stim_param.delay = [5;0;0;0];
% stim_param.mag   = [5;0;0;0];

stim_param.dur(1,1)=50;
stim_param.delay(1,1)=5;
stim_param.mag(1,1)=5; %function input from bugnetworkmodel (Sensor L)
stim_param.mag(2,1)=5; %function input from bugnetworkmodel (Sensor R)

% ------------------------ Neuron Parameters ------------------------------

% model selection: 0 = linear; 1 = quadratic; 2 = exponential
nrn_param.model = 0;

if(nrn_param.model)
    % max allowed potential for quad. and exp. models (mV)
    nrn_param.vth_ss = -20*ones(nonrn,1);
else
    % max allowed potential for linear model (mV)
    nrn_param.vth_ss = -50*ones(nonrn,1);
end

% time constants and electrical properties
nrn_param.tau_a  = 20*ones(nonrn,1);  % adaptation time constant (ms);
nrn_param.del_ga = zeros(nonrn,1);    % conductivity increase per spike (uS)
nrn_param.del_ga(:) = 0;              % no adaptation
nrn_param.v_rest = -65*ones(nonrn,1); % resting potential (mV)
nrn_param.tau_m  = 10*ones(nonrn,1);  % membrane time constant (ms)
nrn_param.r_m    = 10*ones(nonrn,1);  % membrane resistance (Mohm)  


% action potential parameters
v_ap   = 60; % max membrane potential of AP (mV)
dur_rp = 1;  % duration of relative refractory period (ms)

% ------------------------ Network Parameters -----------------------------

% tau rise and decay
syn_param.tau_r = 0.05*ones(nonrn,nonrn);
syn_param.tau_d = 1*ones(nonrn,nonrn);

% conductivity and reversal potentials:(i,j) = neuron i projecting to j
syn_param.gp    = zeros(nonrn,nonrn);  % conductivity matrix
syn_param.erev  = zeros(nonrn,nonrn);  % default excitatory synapses
syn_param.apdel = zeros(nonrn,nonrn); % ap conduction delay
syn_param.apdel(1,2) = 0.5;
% syn_param.apdel(2,1) = 0.3;

syn_param.gp(1,2) = 0.5;
%syn_param.gp(2,3) = 0.5;
%syn_param.gp(3,1)= 0.5; 
syn_param.erev(1,2) = 0;
%syn_param.erev(2,3) = 0;
%syn_param.erev(3,1) = -75;

syn_param.gscale = scale_gpeak(syn_param);

% --------------------- Model Input and Output ----------------------------

nrn = struct();
A_size = [nonrn^2,no_tpts];
nrn.v     = zeros(nonrn,no_tpts); % membrane potential in mV
nrn.g_a   = zeros(nonrn,no_tpts); % dynamic threshold
nrn.upre  = zeros(nonrn,no_tpts); % presynaptic spike times
nrn.upreprop = zeros(nonrn^2,no_tpts); % postsynaptic spike times
nrn.gsyn  = zeros(nonrn^2,no_tpts); % recording values in G
nrn.zsyn  = zeros(nonrn^2,no_tpts); % recording values in Z

nrn.v(:,1) = nrn_param.v_rest(:); % initial condition (resting potential)

fire_ap = ( nrn.v(:,1) > nrn_param.vth_ss(:) );
time_ap = zeros(nonrn,1)-2*dur_rp;
v_int   = zeros(nonrn,1);

% --------------------- Time Marching -------------------------------------

o = ones(1,nonrn);

tic;
for k = 2:no_tpts
    
    t_from_ap = t(k) - time_ap;
    ck_rr     = (t_from_ap >=0) & (t_from_ap < dur_rp);    
   
    v_int = nrn.v(:,k-1);
    if( any(fire_ap) )
        v_int(fire_ap) = nrn_param.v_rest(fire_ap);
    end
    
    % calculate derivative
    [dv,dga,dg,dz] = dydt_network(t(k-1),v_int,nrn.g_a(:,k-1),nrn.upreprop(:,k-1),...
                                  nrn.gsyn(:,k-1),nrn.zsyn(:,k-1));
    
    % take time step
    nrn.v(:,k)    = v_int + dv*dt;
    nrn.g_a(:,k)  = nrn.g_a(:,k-1) + dga*dt;
    nrn.gsyn(:,k) = nrn.gsyn(:,k-1) + dg*dt;
    nrn.zsyn(:,k) = nrn.zsyn(:,k-1) + dz*dt;
    nrn.v(ck_rr,k) = nrn_param.v_rest(ck_rr); % reset v to rest if in refractory period
    
    % see if ap was fired
    fire_ap = ( nrn.v(:,k) > nrn_param.vth_ss(:) );
    if(any(fire_ap))
        
        time_ap(fire_ap)    = t(k); % record ap timestamp
        nrn.v(fire_ap,k)    = v_ap; % set ap v value
        nrn.g_a(fire_ap,k)  = nrn.g_a(fire_ap,k) + nrn_param.del_ga(fire_ap); % adapt
        nrn.upre(fire_ap,k) = 1/dt; % record spike event (i.e. impulse)

        Upre = nrn.upre(:,k)*o;
        syn_fire = (Upre(:) > 0);
        
        s_pre  = find(syn_fire) + nonrn^2*(k-1);          %This section is for creating a delay
        s_move = round( syn_param.apdel(syn_fire)/dt );   % set up for no delay with current
        s_preprop = shift_spike(A_size,s_pre,s_move);        % settings
        nrn.upreprop(s_preprop) = 1/dt;
        
    end

    
end
toc;

% plotting
figure;
subplot(nonrn,1,1);
title('Network Response','FontSize',36,'FontWeight','b');
hold on;
for ii = 1:nonrn
    
    subplot(nonrn,1,ii);
    plot(t,nrn.v(ii,:),'k','LineWidth',2);hold on;
    ylabel(['V',num2str(ii),' (mV)'],'FontSize',30);
    
    if( stim_param.mag(ii) && stim_param.dur(ii) )
        
        i_t = stim_param.delay(ii) + [0,stim_param.dur(ii)];
        plot(i_t,[-90,-90],'b','LineWidth',5);
        L3 = ['I_{mag} = ',num2str(stim_param.mag(ii)),' nA'];
        legend('V_m(t)',L3);
        
    end
    
    set(gca,'FontSize',22);
    ylim([-100,100]);
    xlim([to,tf]);
    setfont
end

figure(2)
plot(t,10.*nrn.gsyn(3,:))
hold on
plot(t, nrn.v(2,:))
plot(t, nrn.v(1,:),'r')
plot(t, nrn.upreprop(3,:),'g')
legend('values in G (Neuron 3)','Neuron 1','Neuron 2','Postsynaptic Spike Times (Neuron 3)',...
    'Location','SW')
set(gca,'FontSize',22);
xlabel('time (ms)','FontSize',30);
setfont

% figure;
% subplot(nonrn,1,1);
% title('Network Response','FontSize',36,'FontWeight','b');
% hold on;
% for ii = 1:nonrn
%     
%     subplot(nonrn,1,ii);
%     plot(t,nrn.g_a(ii,:),'k','LineWidth',3);hold on;
%     ylabel(['ga_',num2str(ii),' (\muS)'],'FontSize',30);
%     
%     set(gca,'FontSize',22);
%     ylim([0,8]);
%     xlim([to,tf]);
%     
% end
% xlabel('time (ms)','FontSize',30);
% 
% return
% 
% figure;
% plot(t,g','LineWidth',3);
% Title('Network Conductance','FontSize',30,'FontWeight','b');
% ylabel('Synaptic Conductance (\muS)','FontSize',30);
% xlabel('Time (ms)','FontSize',30);
% set(gca,'FontSize',28);
% L1 = 'g11 = 0 \muS';
% L2 = 'g21 = 0.5 \muS';
% L3 = 'g31 = 0.2 \muS';
% L4 = 'g12 = 1 \muS';
% L5 = 'g22 = 0 \muS';
% L6 = 'g32 = 0.1 \muS';
% L7 = 'g13 = 0 \muS';
% L8 = 'g23 = 0.05 \muS';
% L9 = 'g33 = 0 \muS';
% legend(L1,L2,L3,L4,L5,L6,L7,L8,L9,'Location','NE');
