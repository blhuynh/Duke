function [V, varargout] = threelayers(V_0,d,h,x,z,y,t,source_param)
%%% Function V = threelayers(V_0,d,h,x,z,y,t)
%%% Calculates the potentials in a three-layer planar volume conductor
%%% model as in L Mesin & R Merletti, "Distribution od Electrical Stimulation
%%% Current in a Planar Miltilayer Anisotropic Tissue", IEEE Trans in
%%% Biomed Eng, vol. 55, no 2., pp.660-670, 2008.
%%% Inputs:
%%%     V_0: voltage or current input in contact with skin layer. 
%%%         It can be a matrix with the 3-D Fourier spatial & time 
%%%         frequency transform of any geometry (e.g. Bessel for an 
%%%         ellipse), or a scalar for point source stimulation at 
%%%         (x,y,z)=(0,0,0). Null frequencies (all dimensions)
%%%         need to be centered (i.e. use fftshift if generated
%%%         with fft Matlab function)
%%%     d: skin thickness (m)
%%%     h: fat thickness (m)
%%%     x,z: spatial sampling in (m) of the y=0 plane. Need to be vectors!
%%%         For ex.: x = (-128:4:124)*1e-3); as used in Mesin & Merletti
%%%     y: depth (m) 
%%%         Interfaces: 
%%%             vacuum-skin: y=0
%%%             skin-fat: y=d
%%%             fat-muscle: y=d+h
%%%     t: time vector in ms, or frequency in Hz for time-harmonic case
%%% Outputs:
%%%     V: voltage (V) in space and time (or phasor if time-harmonic)
%%%     varargout: Current density: Jy, and Jx, Jz (optional)
%%%
%%% Leo Medina, Jan 2013

source_type = source_param.source_type;
Nsamples = length(t);
ifft_mode = 'symmetric';
if Nsamples == 0
    error('Empty time vector');
elseif Nsamples==1
    f = t; %Time-harmonic, argument taken as freq
    if size(V_0,3)>1
        error('Time-harmonic: V_0 shouldnt depend on time');
    end
    ifft_mode = 'nonsymmetric'; %V is a phasor, i.e. may need angle
else
    dt = t(2)-t(1);
    fs = 1/dt*1e3; %Sampling freq in Hz
    NFFT = 2^nextpow2(Nsamples);
    f = fs/2*linspace(0,1,NFFT/2+1);
    f = [-fliplr(f(2:end-1)) f];
end

H = d+h; %coordinate system
w = 2*pi*f;
e_0 = 8.85e-12; %F/m
nfreq = length(f);

%%Tissue parameters
tissue_parameters;
[~, epsilonS, sigmaS] = permittivity(skin,f);
alphaS = sigmaS + 1i*w.*epsilonS*e_0;
[~, epsilonF, sigmaF] = permittivity(fat,f);
alphaF = sigmaF + 1i*w.*epsilonF*e_0;
[~, epsilonMT, sigmaMT] = permittivity(muscle,f); %transversal
alphaM = sigmaMT + 1i*w.*epsilonMT*e_0;
[~, epsilonML, sigmaML] = permittivity(muscle_along,f); %longitudinal
alphaML = sigmaML + 1i*w.*epsilonML*e_0;

%%Spatial frequency domain
[kx, kz, ky] = threelayers_space_sampling(x,z);
nkx = size(kx,2);
nkz = size(kz,1);

%add frequency dimension
ky = repmat(ky,[1 1 nfreq]);
kx = repmat(kx,[1 1 nfreq]);
kz = repmat(kz,[1 1 nfreq]);
alphaS = repmat(reshape(alphaS,[1 1 nfreq]),[nkz nkx 1]);
alphaF = repmat(reshape(alphaF,[1 1 nfreq]),[nkz nkx 1]);
alphaM = repmat(reshape(alphaM,[1 1 nfreq]),[nkz nkx 1]);
alphaML = repmat(reshape(alphaML,[1 1 nfreq]),[nkz nkx 1]);

if strcmpi(source_type,'voltage_eti') %note this is only for disk (mono- and bi-polar)
    %%Eletrode-Skin Interface
    Yeti = electrode_skin_interface(abs(f)); %needs to be normalized per area, see model   
    Yeti = repmat(reshape(Yeti,[1 1 nfreq]),[nkz nkx 1]);
end

%compute kya (note frequency dependency
kya = sqrt(kx.^2+(alphaML./alphaM).*kz.^2);

%Clear unnecessary variables
nout = max(nargout,1) - 1;
if nout<2
    clear kx alphaML;
end

%%Computing coefficients
switch source_type
    case 'voltage'
        Fx = (alphaF.*ky.*tanh(ky*H)+alphaM.*kya)./(alphaM.*kya.*tanh(ky*H)+alphaF.*ky);        
        AfplusBf = alphaS.*V_0./(alphaS+(alphaS-alphaF).*sinh(ky*d).^2+0.5*(alphaF-alphaS).*Fx.*sinh(2*ky*d));        
        Bm = AfplusBf.*(cosh(ky*H)-Fx.*sinh(ky*H))./exp(-kya*H);
        Af = (alphaF.*ky-alphaM.*kya).*exp(-(kya+ky)*H).*Bm./(2*alphaF.*ky);
        Bf = AfplusBf - Af;
        As = .5*(V_0.*(1-coth(ky*d))+AfplusBf.*(coth(ky*d)-Fx));
        Bs = V_0 - As;
        
    case 'voltage_eti'        
        Fx = (alphaF.*ky.*tanh(ky*H)+alphaM.*kya)./(alphaM.*kya.*tanh(ky*H)+alphaF.*ky);
        Fx(isinf(Fx)) = exp(100);
        AfplusBf = -alphaS.*Yeti.*V_0./( sinh(ky*d).^2.*(alphaF.*Yeti+alphaS.^2.*ky.*Fx) - cosh(ky*d).^2.*(alphaS.*alphaF.*ky.*Fx+alphaS.*Yeti) + sinh(ky*d).*cosh(ky*d).*(alphaF.*alphaS.*ky + (alphaS-alphaF).*Fx.*Yeti - alphaS.^2.*ky) );
        AfplusBf(isinf(AfplusBf)) = exp(100);
        Bm = AfplusBf.*(cosh(ky*H)-Fx.*sinh(ky*H))./exp(-kya*H);
        Bm(isinf(Bm)) = exp(100);
        Af = (alphaF.*ky-alphaM.*kya).*exp(-kya*H).*Bm./(2*alphaF.*ky.*exp(ky*H));
        Af(isinf(Af)) = exp(100);
        Bf = AfplusBf - Af;
        As = .5*(V_0+(Yeti+alphaS.*ky).*(AfplusBf.*(cosh(ky*d)-Fx.*sinh(ky*d))-V_0.*cosh(ky*d))./(Yeti.*sinh(ky*d)+alphaS.*ky.*cosh(ky*d)));        
        As(isinf(As)) = exp(100);
        Bs = (V_0.*Yeti + As.*(alphaS.*ky - Yeti))./(Yeti + alphaS.*ky);        
        Bs(isinf(Bs)) = exp(100);
        
    case 'current'
        Fx = (alphaF.*ky.*tanh(ky*H)+alphaM.*kya)./(alphaM.*kya.*tanh(ky*H)+alphaF.*ky);
        AfplusBf = (V_0./ky)./(alphaF.*Fx+(alphaS-alphaF).*sinh(2*ky*d)/2+(alphaF-alphaS).*Fx.*sinh(ky*d).^2);
        Bm = AfplusBf.*(cosh(ky*H)-Fx.*sinh(ky*H))./exp(-kya*H);
        Af = (alphaF.*ky-alphaM.*kya).*exp(-(kya+ky)*H).*Bm./(2*alphaF.*ky);
        Bf = AfplusBf - Af;
        As = .5*((V_0./(alphaS.*ky)).*(tanh(ky*d)-1)+AfplusBf.*(1-Fx.*tanh(ky*d)));
        Bs = V_0./(ky.*alphaS) + As;
        
    otherwise
        error([source_type ': Unknown source type']);
end

%%Computing potential by inverting 2D Fourier transfor at y=y(i) plane
V = zeros(length(y),nkz,nkx,nfreq);
if nout>0
    Jy = zeros(length(y),nkz,nkx,nfreq);
    if nout==3
        Jx = zeros(length(y),nkz,nkx,nfreq);
        Jz = zeros(length(y),nkz,nkx,nfreq);
    end
    if nout~=1 && nout~=3
        error('Incorrect number of outputs: it can be 1, 2 or 4');
    end
end
for i=1:length(y)   
    if y(i)>=0 && y(i)<d %skin
        aux = As.*exp(ky*y(i))+Bs.*exp(-ky*y(i));        
        if exist('Jy','var')
            auxJy = -alphaS.*ky.*(As.*exp(ky*y(i))-Bs.*exp(-ky*y(i)));
        end
        if exist('Jx','var') || exist('Jz','var')
            alpha_aux = alphaS;
            alpha_aux_z = alphaS;
        end
    elseif y(i)>=d && y(i)<H %fat
        aux = Af.*exp(ky*y(i))+Bf.*exp(-ky*y(i));
        if exist('Jy','var')
            auxJy = -alphaF.*ky.*(Af.*exp(ky*y(i))-Bf.*exp(-ky*y(i)));
        end
        if exist('Jx','var') || exist('Jz','var')
            alpha_aux = alphaF;
            alpha_aux_z = alphaF;
        end
    elseif y(i)>=H %muscle
        aux = Bm.*exp(-kya*y(i));
        if exist('Jy','var')
            auxJy = alphaM.*kya.*Bm.*exp(-kya*y(i));
        end
        if exist('Jx','var') || exist('Jz','var')
            alpha_aux = alphaM;
            alpha_aux_z = alphaML;
        end
    else
        warning(['Not evaluated at y = ' num2str(y(i))]);
    end
    if (strcmpi(source_type,'voltage') || strcmpi(source_type,'voltage_eti'))
        aux(ky==0) = V_0(ky==0);
    end
    ind0 = find(isnan(aux) | isinf(aux));
    for j=1:length(ind0)
        ind = [-4:-1] + ind0(j); %Note Matlab goes by columns
        if any(ind<=0)
            aux(ind0(j)) = 0;
        elseif ~issorted(kz(ind))
            aux(ind0(j)) = 0;
        else
            aux(ind0(j)) = interp1(abs(kz(ind)),abs(aux(ind)),kz(ind0(j)),'cubic','extrap');
        end
    end    
    if exist('Jy','var')        
        ind0 = find(isnan(auxJy) | isinf(auxJy));
        for j=1:length(ind0)    
            ind = [-4:-1] + ind0(j); %Note Matlab goes by columns
            if any(ind<=0)
                auxJy(ind0(j)) = 0;
            elseif ~issorted(kz(ind))
                auxJy(ind0(j)) = 0;
            else
                auxJy(ind0(j)) = interp1(abs(kz(ind)),abs(auxJy(ind)),abs(kz(ind0(j))),'cubic');
            end
        end
    end
 
    V(i,:,:,:) = ifftn(ifftshift(aux),ifft_mode);
    if exist('Jy','var')
        Jy(i,:,:,:) = ifftn(ifftshift(auxJy),ifft_mode);
    end
    if exist('Jx','var') || exist('Jz','var')
        Jx(i,:,:,:) = ifftn(ifftshift(1i.*kx.*alpha_aux.*aux),ifft_mode);
        Jz(i,:,:,:) = ifftn(ifftshift(1i.*kz.*alpha_aux_z.*aux),ifft_mode);
    end
end

V = V(:,:,:,1:Nsamples);
if exist('Jy','var')
    Jy = Jy(:,:,:,1:Nsamples);
    varargout(1) = {Jy};
end
if exist('Jx','var') || exist('Jz','var')
    Jx = Jx(:,:,:,1:Nsamples);
    Jz = Jz(:,:,:,1:Nsamples);
    varargout(2) = {Jx};
    varargout(3) = {Jz};
end