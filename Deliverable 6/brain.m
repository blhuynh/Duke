function [motorL, motorR] = brain(sensorL, sensorR)
%Inputs
% sensorL, sensorR: magnitude of the signal (distance from food)
%
%Outputs
% motorL, motorR: magnitude of output motor

exponent=1;
scalar=500;
% compute a relationship of sensor input to motor output
motorL=(sensorL*scalar)^exponent;      
motorR=(sensorR*scalar)^exponent;
% fprintf([num2str([motorL motorR]),'\n'])