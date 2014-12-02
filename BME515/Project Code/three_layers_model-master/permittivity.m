function [epsilon, e, s] = permittivity(tissue,freq)
%epsilon: complex relative permittivity
%e: relative permittivity
%s: conductivity

e_0 = 8.85e-12; %F/m
x = tissue;

freq(freq==0) = 1; %approx d.c. as a 1Hz signal to avoid Inf
freq = abs(freq);
epsilon = x.e_inf + x.sigma./(1i*2*pi*freq*e_0);
for j = 1:4
    s = num2str(j);
    a = x.(['alpha_' s]);
    delta =  x.(['de_' s]);
    tau =  x.(['tau_' s]);
    epsilon = epsilon + delta./(1 + (1i*2*pi*freq*tau).^(1-a));
end

% e = real(epsilon);
e = zeros(size(epsilon));
% s = -imag(epsilon)*2*pi.*freq*e_0;
s = x.sigma*ones(size(epsilon));