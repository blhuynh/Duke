format compact
clear
close all

fntsz = 14;
lnwdth = 3;
scrnsz = get(0,'ScreenSize');

%% Data parameters
D_values = 2:2:20;         % Fiber diameters [um]

%% Read in the Im(t) data
for D_ind = 1:numel(D_values)
   fname = ['Im_' num2str(D_values(D_ind)) 'um.dat'];
   tmp = importdata(fname);
   if (D_ind==1)
      t = tmp.data(:,1);
   end
   Im_data(D_ind,:) = tmp.data(:,2);
end

% Downsample to dt=20us instead of 1us
%t = t(1:20:end);
%Im_data = Im_data(:,1:20:end);

%% Obtain peak current values
Im_peaks = max(Im_data,[],2);
% Fit data to a line
p = polyfit(D_values',Im_peaks,1);

%% Plot
% Plot Im(t) for all D values
figure
set(gca,'fontsize',fntsz)
plot(t,Im_data,'linewidth',lnwdth)
xlabel('Time (ms)')
ylabel('Im (nA)')
legend('D=2um','D=4um','D=6um','D=8um','D=10um','D=12um','D=14um',...
   'D=16um','D=18um','D=20um')
xlim([6 15])

% Plot peak Im vs fiber diameter
figure
hold on
set(gca,'fontsize',fntsz)
plot(D_values,Im_peaks,'ko','markersize',10,'linewidth',lnwdth)
plot(D_values, p(1).*D_values + p(2),'k--','linewidth',lnwdth)
xlabel('D (um)')
ylabel('max(Im) (nA)')
h = legend('Data from NEURON',...
   ['Trendline: max(Im)[nA]=' num2str(p(1)) '*D[um] + ' num2str(p(2))],...
   'location','northwest');
set(h,'fontsize',12);