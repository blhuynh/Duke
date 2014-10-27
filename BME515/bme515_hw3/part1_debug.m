% BME 515 HW 3 Part 1
% 30-Oct-2014 (blh19)

clear; clc
nocols = 10; 
fid1 = 'part1_debug.txt';
data = fileload(fid1);
data = reshape(data,nocols,numel(data)/nocols)';

ina_node50 = data(:,1);     % last node
ina_node25 = data(:,2);     % center node
ik_node50 = data(:,3);      % last node
ik_node25 = data(:,4);      % center node
il_node50 = data(:,5);      % last node
il_node25 = data(:,6);      % center node
vm_node50 = data(:,7);      % last node (mV)
vm_node25 = data(:,8);      % center node (mV)
i_membrane(:,1) = data(:,9);
i_membrane(:,2) = data(:,10); 

tstop = 25;
dt = 0.02; % (ms)
tvec = 0:dt:tstop-dt;

cm = 2; % (uF/cm^2)

% calculate capacitive current
ic_node50 = diff(vm_node50)*cm/dt*1e-3; % (uF/cm^2 * mV/ms) -> (mA/cm^2)
ic_node25 = diff(vm_node25)*cm/dt*1e-3; 

% total transmembrane current
itot_node25 = ina_node25+ik_node25+il_node25;
itot_node25(1:length(ic_node25)) = itot_node25(1:length(ic_node25))+ic_node25;
itot_node50 = ina_node50+ik_node50+il_node50;
itot_node50(1:length(ic_node50)) = itot_node50(1:length(ic_node50))+ic_node50;

%% Transmembrane Current
% 1) capacitive current
% 2) transient inward current
% 3) delayed outward current
% ina inward, ik outward

figure(1); clf; hold on
subplot(2,1,1); hold on
plot(tvec,ina_node25,'k')
plot(tvec,ik_node25,'b')
plot(tvec,il_node25,'g')
plot(tvec(1:length(ic_node25)),ic_node25,'c')
plot(tvec,itot_node25,'r')
legend('I_{Na}','I_K','I_{leak}','I_c','I_m','Location','NE')
xlabel('Time (ms)')
ylabel('Membrane Current (mA)') % does Nikki want mA or mA/cm^2 ??
title('Center Node')
setfont(24)

subplot(2,1,2); hold on
plot(tvec,ina_node50,'k')
plot(tvec,ik_node50,'b')
plot(tvec,il_node50,'g')
plot(tvec(1:length(ic_node50)),ic_node50,'c')
plot(tvec,itot_node50,'r')
legend('I_{Na}','I_K','I_{leak}','I_c','I_m','Location','NE')
xlabel('Time (ms)')
ylabel('Membrane Current (mA)') % does Nikki want mA or mA/cm^2 ??
title('Last Node')
setfont(24)

%% Membrane Voltage
figure(2); clf; hold on
plot(tvec,vm_node25,'k')
plot(tvec,vm_node50,'b')
legend('Node 25','Node 50')
xlabel('Time (ms)')
ylabel('V_m(t)')
title('Membrane Voltage')
setfont(24)

%% i_membrane
figure(3); clf; hold on
plot(tvec,i_membrane(:,1),'k')
plot(tvec,i_membrane(:,2),'b')
legend('Close to Electrode','Far from Electrode')
title('i_{membrane}')
xlabel('Time (ms)')
ylabel('I_m(t)')
axis tight
setfont(24)