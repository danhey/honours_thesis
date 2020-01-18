%--------------------------------------------------------------------%
%------------ Transmission spectra of Etalon ------------------------%
%------------ with PML Boundary conditions---------------------------%
%------------ all units in eV,fs, nm --------------------------------%
%--------------------------------------------------------------------%

clc
close all
clear all

%--------------------------------------------------------------------%
h=0.6582119514 ; % Planck's constant
dt=0.001; % time step in femto seconds
tmax=30; % maximum time limit
t=-25:dt:tmax; % time array
c0=299.792458; % speed of light in nm/fs
S=0.5; % courant factor should be less than 1
dz=c0*dt/S; % step size determined by time step
mu0=2.013354451e-4; % permeability of free space in (V fs^2/e nm)
ep0=55.26349597e-3; % permittivity of free space in (e / V nm)

M=500; % no. of spatial grid points
Z=(0:M-1).*dz; % space axis
ep(1:M)=ep0; % permitivity array
mu(1:M)=mu0; % pemeability array

%-----------------------------------------
%-----------insert device

Mstart=200;
Mend=300;
RP=32; % point where reflection is measured
TP=468;  % point where transmission is measured
ep(Mstart:Mend)=13.0*ep0; % insert etalon

%---- PML absorbing boundary condition----

sigma(1:M)=0; % initialize conductivity array
d=30; % width of PML layer
m=3; % polynomial order for grading sigma array (pp 292, Taflove)
neta=sqrt(mu0/ep0); 
R=1e-8; % required reflectivity
sigma_max=-(m+1)*log(R)/(2*neta*d*dz);
Pright=((1:d+1)./d).^m*sigma_max;
sigma(M-d:M)=Pright; % lossy conductivity profile
sigma(1:d+1)=fliplr(Pright);
sigma_star(1:M)=sigma.*mu0./ep0; % Eq 7.8 Taflove, pp 275

%------------- PML constants ----------------------------------------%

A=((mu-0.5*dt*sigma_star)./(mu+0.5*dt*sigma_star)); 
B=(dt/dz)./(mu+0.5*dt*sigma_star);                          
C=((ep-0.5*dt*sigma)./(ep+0.5*dt*sigma)); 
D=(dt/dz)./(ep+0.5*dt*sigma);                     

%!-----------------Source parameters---------------------------------%

wlev=1.5d0 ;   % freq of light in eV
wl=(wlev/h) ;  % freq of light in (1/fs)
delta=5.0d0 ; % width of electric field
Zs=70; % position index of source
Esrc=exp(-((t).^2)/(delta^2)).*cos(wl.*t); % Electric field source
Am=sqrt(ep(Zs)/mu(Zs));
s=((S*dz)/c0)-dt/2;
Hsrc=Am*exp(-((t+s).^2)/(delta^2)).*cos(wl.*(t+s));

%-------------------------------------
% Initialize Fourier Transform, Ref. Dr.Raymond Rumph, UTEP
fmax=10;
Nfreq=500;
FREQ=linspace(0,fmax/h,Nfreq);
K=exp(-1i*dt*FREQ);
REF=zeros(1,Nfreq);
TRN=zeros(1,Nfreq);
SRC=zeros(1,Nfreq);

%---------------------------------------------------------------------%
%------initialize fields in space array

Hy(1:M)=0.0; 
Ex(1:M)=0.0;

%---- begin time loop----
fh = figure(1);
set(fh, 'Color', 'white'); 

for n=1:length(t)
        Hy(1:M-1)=A(1:M-1).*Hy(1:M-1)-B(1:M-1).*(Ex(2:M)-Ex(1:M-1));
        Hy(Zs-1)=Hy(Zs-1)-B(Zs-1)*Esrc(n);
        Ex(2:M-1)=C(2:M-1).*Ex(2:M-1)-D(2:M-1).*(Hy(2:M-1)-Hy(1:M-2));
        Ex(Zs)=Ex(Zs)-D(Zs)*Hsrc(n);
        
        if mod(n,50)==0
            plot(Hy)
            plot(Ex)
            drawnow();
        end
    for nf = 1 : Nfreq 
     REF(nf) = REF(nf) + (K(nf)^n)*Ex(RP); 
     TRN(nf) = TRN(nf) + (K(nf)^n)*Ex(TP); 
     SRC(nf) = SRC(nf) + (K(nf)^n)*Esrc(n); 
    end 
end

%-----------------------------------------------------------------------%
REF = REF*dt; 
TRN = TRN*dt; 
SRC = SRC*dt;

REF = abs(REF./SRC).^2; 
TRN = abs(TRN./SRC).^2;
CON = REF+TRN;

%-----------------------------------------------------------------------%

figure(1)
plot(FREQ*h,REF)
hold on
plot(FREQ*h,TRN,'r-')
plot(FREQ*h,TRN+REF,'k-')
xlabel('E(eV)')
ylabel('Amplitude')
title('Impulse response of Etalon')
legend on
legend('Rx','Tx','Tx+Rx')
hold off
fh = figure(1);
set(fh, 'color', 'white');

%-------------------------------------------------------------------------