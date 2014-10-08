function [motorL,motorR,v] = brain_avoid(sensorL,sensorR,DT,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two neuron controller for BraitenBug
%
% Avoids the light source
%
% Neuron 1 receives input from the LEFT sensor and drives the LEFT motor
% Neuron 2 receives input from the RIGHT sensor and drives the RIGHT motor
%
% Neurons are modeled as leaky integrate-and-fire units
%
% Mark Nelson, Feb 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent v % membrane voltage

% fixed parameters
nNeurons = 2;
GAIN = 0.075;
VTHR = 5;
VSPK = 20;
TAU = 0.3; 


% get sensor input 
% (sensor input range: 0.0 - 1.0)
input = [sensorL sensorR];

% add some steady-state input current, otherwise 
% the bug won't move when it's far from the light, 
% resulting in a very boring simulation! 

% loop over neurons
for n = 1:nNeurons
    
    % update membrane voltage
    if(v(n)==VSPK)
        v(n) = 0;
    else
        dvdt = (v(n) + GAIN*input(n))/TAU;
        v(n) = v(n) + dvdt*DT;
    end
    
    % check spike threshold
    if(v(n) >= VTHR)
        v(n)= VSPK;
    end
end

% drive the motors with spike output
motorL = (v(1)==VSPK);
motorR = (v(2)==VSPK);