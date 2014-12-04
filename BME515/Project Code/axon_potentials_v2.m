%% Initialize 
clear

a = 3.175/2;   % electrode radius in cm 
% r = -a:0.01:a; % coordinate along electrode in cm
r = load('nodelocations.txt') / 10; % (mm -> cm)
z = 2;         % depth coordinate in cm
I = 25;        % current in mA
sigma = 1/12;  % tissue conductivity in S/m
 % convert to S/cm
 sigma = sigma*(100); % tissue conductivity in S/cm

%% Current density
J = I./(pi.*a^2); % mA/cm^2

%% Voltage at electrode
V0 = zeros(1,length(r));
for i = 1:length(r)
    V0(i) = (pi*J.*((a^2-r(i).^2).^0.5))./(2*sigma);   % units of mV 
end

% plot to check
figure(1); clf
plot(r, V0)
xlabel('Coordinate r along electrode (cm)')
ylabel('Potential at coordinate r mV)')
title('Potential V_{0} at electrode')

%% Voltage within medium
V = zeros(1,length(r));
for i = 1:length(r)
    V(i) = (2*V0(i)./pi).*asin((2*a)./(sqrt((r(i)-a).^2+z.^2)+sqrt((r(i)+a).^2+z.^2)));
end

% Plot V(r,z) to check
figure(2); clf
plot(r,V)
xlabel('Coordinate r along electrode (cm)')
ylabel('Potential at coordinate r (mV)')
title('Potential within a medium at depth z = ?? cm')