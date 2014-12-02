format compact
clear
close all

fntsz = 14;
lnwdth = 3;
scrnsz = get(0,'ScreenSize');

%% CV data from NEURON
D_values = 2:2:20;         % Fiber diameters [um]
CV = [2.644 5.287 8.023 10.698 13.372 16.047 18.721 21.647 24.07 26.744];

%% Fit data to a line
p = polyfit(D_values,CV,1);

%% Plot
figure
hold on
set(gca,'fontsize',fntsz)
plot(D_values,CV,'ko','markersize',10,'linewidth',lnwdth)
plot(D_values, p(1).*D_values + p(2),'k--','linewidth',lnwdth)
xlabel('D (um)')
ylabel('CV (m/s)')
h = legend('Data from NEURON',...
   ['Trendline: CV[m/s]=' num2str(p(1)) '*D[um] + ' num2str(p(2))],...
   'location','northwest');
set(h,'fontsize',12);
