% Create axon potentials across two-layer thing
% 1-Dec-2014 (blh19)

% clear all;
clc;

%% ------------------------- DC fiber body --------------------------------

D = 1;
arclength=150*1e3; % mm -> um

numnode=1;
nummysa=2;
numflut=2;
numstin=3;

nodelength=1;
mysalength=3;
flutlength=0.0003*D^6 - 0.0163*D^5 + 0.3425*D^4 - 3.4752*D^3 + 17.1212*D^2 - 31.4582*D + 27.4112;
inl=0.1207*D^4 - 5.0396*D^3 + 67.8440*D^2 - 219.7125*D + 375.3634;
stinlength=(inl-nodelength-nummysa*mysalength-numflut*flutlength)/numstin;

numele_inl=8;
nonodes=floor(arclength/inl)+1;
numtotele=(nonodes-1)*numele_inl+1;

%% ------------------------- coordinates ----------------------------------

p1=nodelength/2+mysalength/2;
p2=mysalength/2+flutlength/2;
p3=flutlength/2+stinlength/2;
oneseg=[p1,p2,p3,stinlength*ones(1,numstin-1),p3,p2,p1];
allseg=[0,repmat(oneseg,[1,nonodes-1])];
x_allseg=1e-3*cumsum(allseg); % um -> mm

if(mod(numtotele,2)~=0)
    midnode=(numtotele+1)/2;
    x_allseg=x_allseg-x_allseg(midnode);
else
    midnode=(numtotele+1)/2;
    x_allseg=x_allseg-x_allseg(end)/2;
end

xyz_axon=zeros(numtotele,3);
xyz_axon(:,1)=0; % x location
xyz_axon(:,2)=0; % y location
xyz_axon(:,3)=x_allseg; % z location
xyz_axon=xyz_axon';

% plot(xyz_axon(3,:),zeros(size(xyz_axon(3,:))),'k.');
% plot3(xyz_axon(2,:),xyz_axon(2,:),xyz_axon(3,:),'k.');

%% -------------------- Calculate Potentials ------------------------------
% conversion ratios
mmPerInch = unitsratio('mm','in');

skinthickness = 2.25; % (mm)
fatthickness = 18; % (mm)
sigma1=1/1000 * 1e-3; % skin conductivity (S/m -> S/mm);
sigma2=1/27.33 * 1e-3; % fat conductivity
z = skinthickness+fatthickness;
sigmalumped = 1/(skinthickness/(z*sigma1) + fatthickness/(z*sigma2));
iunit=25e-3; % (mA -> A)

a = 1.25/2*mmPerInch;   % electrode radius in mm 
I = 1;                  % current in mA

%% Wiley and Webster
% Current density
J = I./(pi.*a^2); % mA/mm^2

r = abs(xyz_axon(3,:));

% calculate potential at electrode-medium boundary
V0 = zeros(size(r));
V0(r<a) = pi * J * sqrt(a^2-r(r<a).^2) / (2*sigmalumped) * 1e-3; % mV -> V

% calculate potential at fiber
V = 2/pi * max(V0) .* asin(2*a ./ ( sqrt((r-a).^2+z^2) + sqrt((r+a).^2+z^2) ) );

% figure(1); clf; hold on
% ax(1) = subplot(2,1,1);
% plot(xyz_axon(3,:),V0,'k')
% title('Electrode Surface Potentials (z=0mm)')
% xlabel('Distance (r)'); ylabel('Voltage (V)')
% gzoom; setfont(18)
% 
% ax(2) = subplot(2,1,2);
% plot(xyz_axon(3,:),V,'k')
% title(sprintf('Axon Potentials (z=%gmm)',z))
% xlabel('Distance (r)'); ylabel('Voltage (V)')
% gzoom; setfont(18)
% 
% linkaxes(ax,'xy')
% print -dpng potentials

figure(3); clf; hold on
plot(xyz_axon(3,:),V,'k')
title(sprintf('Axon Potentials (z=%gmm)',z))
xlabel('Distance (r)'); ylabel('Voltage (V)')
gzoom; setfont(18)
print -dpng potentials


dlmwrite('phiaxon.txt',V,'delimiter',' ','precision',8);

%% activating function
vnor = V(1:8:end);
dnor = xyz_axon(3,1:8:end);

d2 = diff(diff(vnor));

figure(2); clf; hold on
plot(dnor(2:end-1),d2,'k')
title('Activating Function')
xlabel('Position (mm)'); ylabel('\delta^2 V')
axis tight; gzoom
setfont(18)
print -dpng activatingfunction