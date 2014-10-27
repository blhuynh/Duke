% BME 515 HW 3 Part 1
% 30-Oct-2014 (blh19)

clear; clc
fiberDlist = 2:2:20;
workdir = pwd();

nodelength = 1; % (um)

cd Data
for k = 1:length(fiberDlist)
    D = fiberDlist(k);
    fid = sprintf('part1_%gum.txt',D);
    Im(:,k) = fileload(fid) * (nodelength*pi*D) * 1e-8; % (mA/cm^2) * (um)^2 -> (mA)
end
cd(workdir)

tstop = 25;
dt = 0.02; % (ms)
tvec = 0:dt:tstop-dt;

%% Transmembrane Current
% 1) capacitive current
% 2) transient inward current
% 3) delayed outward current
% ina inward, ik outward

for k=1:length(fiberDlist)
    lgd{k} = sprintf('%g um',fiberDlist(k));
end
plottypes = {'k','y','m','c','r','g','b','k--','k:','k-.'};

figure(1); clf; hold on
for k=1:length(fiberDlist)
    plot(tvec,Im(:,k),plottypes{k})
end

xlabel('Time (ms)')
ylabel('I_m (mA)')
title('Transmembrane Current')
legend(lgd,'Location','NW')
axis tight
setfont(24)
print -dpng bme515_hw3_part1c

%% Peak Positive Current
% calculate trendline
P = polyfit(fiberDlist,max(Im),1);
Y = polyval(P,fiberDlist);
[r2,rmse] = rsquare(max(Im),Y);

trendline = sprintf('y=%.2gx+%.2g',P(1),P(2));
stats = sprintf('R2=%.2g RMSE=%.2g',r2,rmse);
fprintf([trendline,'\n',stats,'\n'])

figure(2); clf; hold on
plot(fiberDlist,max(Im),'k.')
plot(fiberDlist,Y,'k--')
legend('Data','Fit','Location','NW')
xlabel('Fiber Diameter (um)')
ylabel('I_m (mA)')
title({'Peak Positive Current';trendline;stats})
setfont(24)

print -dpng bme515_hw3_part1d