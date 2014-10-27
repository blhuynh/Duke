% BME 515 HW 3 Part 2 - Conduction Speed
% 30-Oct-2014 (blh19)

clear; clc

%% load data
data = fileload('part2.txt');
data = reshape(data,2,numel(data)/2)';

%% trendline
P = polyfit(data(:,1),data(:,2),1);
X = linspace(min(data(:,1)),max(data(:,1)));
Y = polyval(P,X);

%% R2
[r2,rmse] = rsquare(data(:,2),polyval(P,data(:,1)));

trendline = sprintf('y=%.2gx+%.2g',P(1),P(2));
stats = sprintf('R2=%.2g RMSE=%.2g',r2,rmse);
fprintf([trendline,'\n',stats,'\n'])
%% plot data
% N.B. This is not representative of a mammalian peripheral nerve axon.

figure(1); clf; hold on
plot(data(:,1),data(:,2),'k.')
plot(X,Y,'k')
legend('Data','Fit','Location','NW')
title({'Conduction Velocity for Varying Fiber Diameters',trendline,stats})
xlabel('Fiber Diameter (\mum)')
ylabel('Conduction Velocity (m/s)')
setfont(24)

print -dpng bme515_hw3_part2