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
% intradural = 17.7, extradural = 22.7

mono=0;
if(mono==1)
    cfac=0;
    ptag='mp';
else
    cfac=1;
    ptag='bp';
end

% conversion ratios
mmPerInch = unitsratio('mm','in');

xe1=0;
ye1=1;
ze1=0;
xe2=0;
ye2=ye1;
ze2=0.5*mmPerInch; % (inch->mm)

sigma1=16e-3 * 1e-3; % skin conductivity (S/m -> S/mm);
sigma2=2733e-2 * 1e-3; % fat conductivity
iunit=1e-3; % (mA -> A)

re1=[xe1;ye1;ze1];
re2=[xe2;ye2;ze2];

% DC body
r1_axon=sqrt( sum((xyz_axon-re1*ones(1,numtotele)).^2,1) );
r2_axon=sqrt( sum((xyz_axon-re2*ones(1,numtotele)).^2,1) );
phi_axon=(iunit/4/pi/(sigma1+sigma2))*(1./r1_axon-cfac*1./r2_axon);

plot(xyz_axon(3,:),phi_axon,'k');

dlmwrite(['phiaxon_',ptag,'_',num2str(D),'um_ye',num2str(ye1),'mm.txt'],phi_axon,'delimiter',' ','precision',8);
