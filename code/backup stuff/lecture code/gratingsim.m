% Lecture22_sawtooth.m
% INITIALIZE MATLAB
close all;
clc;
clear all;
% UNITS
meters = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches = 2.54 * centimeters;
feet = 12 * inches;
seconds = 1;
hertz = 1/seconds;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;
% CONSTANTS
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;
% SOURCE PARAMETERS
NFREQ = 500;
fmax = 15 * gigahertz;
FREQ = linspace(5,15,NFREQ) * gigahertz;
f0 = 10 * gigahertz;
lam0 = c0/f0;
% GRATING PARAMETERS
L = 1.5 * centimeters;
d = 1.0 * centimeters;
er1 = 1.0;
er2 = 9.0;
% GRID PARAMETERS
nmax = sqrt(max([er1 er2]));
NRES = 10;
NPML = [0 0 20 20];
BUF = 0.5*lam0 * [1 1];

% COMPUTE INITIAL GRID RESOLUTION
lam0_min = c0/max(FREQ);
dx = lam0_min/nmax/NRES;
dy = lam0_min/nmax/NRES;

% SNAP GRID TO CRITIAL DIMENSIONS
Nx = 2*ceil(L/dx/2) + 1;
dx = L/Nx;
Ny = ceil(d/dy);
dy = d/Ny;

% COMPUTE GRID SIZE
Sx = L;
Sy = sum(BUF) + d;
Ny = ceil(Sy/dy) + NPML(3) + NPML(4) + 5;
Sy = Ny*dy;

% INITIALIZE MATERIALS TO FREE SPACE
URxx = ones(Nx,Ny);
URyy = ones(Nx,Ny);
ERzz = er1 * ones(Nx,Ny); 

% COMPUTE POSITION INDICES
ny1 = NPML(3) + 3 + round(BUF(1)/dy);
ny2 = ny1 + round(d/dy) - 1;

% ADD GRATING
for ny = ny1 : ny2
f = (ny - ny1 + 1)/(ny2 - ny1 + 2);
nx = round(f*Nx);
nx2 = Nx;
nx1 = nx2 - nx + 1;
ERzz(nx1:nx2,ny) = er2;
end

ERzz(:,ny2+1:Ny) = er2;
% COMPUTE STABLE TIME STEP
dmin = min([dx dy]);
dt = dmin/(2*c0);
% SOURCE POSITION
ny_src = NPML(3) + 2;

% COMPUTE TIME PARAMETERS
tau = 0.5/fmax;
t0 = 6*tau;
A = sqrt(ERzz(1,ny_src)/URyy(1,ny_src));
delt = 0.5*dy/c0 + dt/2;

% TIME STEPS
proptime = nmax*Sy/c0;
simtime = 2*t0 + 10*proptime;
STEPS = ceil(simtime/dt);

% COMPUTE GAUSSIAN SOURCES
ta = [0:STEPS-1]*dt;
Ez_src = exp(-((ta-t0)/tau).^2);
Hx_src = A*exp(-((ta-t0+delt)/tau).^2);

% COMPUTE FOURIER TRANSFORM KERNELS
K = exp(-1i*2*pi*dt*FREQ);
K0 = exp(-1i*2*pi*dt*f0);

% INITIALIZE STEADY-STATE FIELDS
Eref = zeros(Nx,NFREQ);
Etrn = zeros(Nx,NFREQ);
SRC = zeros(1,NFREQ);
Eref0 = zeros(Nx,1);
Etrn0 = zeros(Nx,1);
ssSRC = 0;

% POSITION OF RECORD PLANES
ny_ref = NPML(3) + 1;
ny_trn = Ny - NPML(4);

% COMPUTE REFRACTIVE INDICES IN RECORD PLANES
nref = sqrt(ERzz(1,ny_ref)*URxx(1,ny_ref));
ntrn = sqrt(ERzz(1,ny_trn)*URxx(1,ny_trn));

%% Compute PML
% 2x grid formulation
Nx2 = 2*Nx;
Ny2 = 2*Ny;

% COMPUTE sigx
sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end

% COMPUTE sigy
sigy = zeros(Nx2,Ny2);
for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

%% Compute update coefficients
% These are outside the main loop so that they're not calculated at every
% time-step

% Hx coefficients
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./URxx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./URxx ./ mHx0;

% Hy coefficients
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./URyy ./ mHy0;

% Dz coefficients
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

% Ez coefficient
mEz1 = 1./ERzz;

% INITIALIZE FIELDS
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
Dz = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);
% INITIALIZE CURL
CEx = zeros(Nx,Ny);
CEy = zeros(Nx,Ny);
CHz = zeros(Nx,Ny);
% INITIALIZE INTEGRATION TERMS
ICEx = zeros(Nx,Ny);
ICEy = zeros(Nx,Ny);
IDz = zeros(Nx,Ny);

for T = 1:STEPS
    % Compute CEx
for ny = 1 : Ny-1
for nx = 1 : Nx
CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy;
end
end
for nx = 1 : Nx
CEx(nx,Ny) = (Ez(nx,1) - Ez(nx,Ny))/dy;
end

% Compute CEy
for nx = 1 : Nx-1
for ny = 1 : Ny
CEy(nx,ny) = - (Ez(nx+1,ny) - Ez(nx,ny))/dx;
end
end
for ny = 1 : Ny
CEy(Nx,ny) = - (Ez(1,ny) - Ez(Nx,ny))/dx;
end

% TF/SF
CEx(:,ny_src-1) = CEx(:,ny_src-1) - Ez_src(T)/dy;

% Update Integration Terms
ICEx = ICEx + CEx;
ICEy = ICEy + CEy;

% Update Hx and Hy
Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;

% Compute CHz
CHz(1,1) = (Hy(1,1) - Hy(Nx,1))/dx ...
- (Hx(1,1) - Hx(1,Ny))/dy;
for nx = 2 : Nx
CHz(nx,1) = (Hy(nx,1) - Hy(nx-1,1))/dx ...
- (Hx(nx,1) - Hx(nx,Ny))/dy;
end
for ny = 2 : Ny
CHz(1,ny) = (Hy(1,ny) - Hy(Nx,ny))/dx ...
- (Hx(1,ny) - Hx(1,ny-1))/dy;
for nx = 2 : Nx
CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1,ny))/dx ...
- (Hx(nx,ny) - Hx(nx,ny-1))/dy;
end
end

% TF/SF
CHz(:,ny_src) = CHz(:,ny_src) + Hx_src(T)/dy;

% Update Integration Term
IDz = IDz + Dz;
% Update Dz
Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
% Update Ez
Ez = mEz1.*Dz;

% Update Fourier Transforms
for nfreq = 1 : NFREQ
    Eref(:,nfreq) = Eref(:,nfreq) + (K(nfreq)^T)*Ez(:,ny_ref)*dt;
    Etrn(:,nfreq) = Etrn(:,nfreq) + (K(nfreq)^T)*Ez(:,ny_trn)*dt;
    SRC(nfreq) = SRC(nfreq) + (K(nfreq)^T)*Ez_src(T)*dt;
end

if mod(T,50)==0
imagesc(Ez);
drawnow();
end
end

% % LOOP OVER FREQUENCY
% for nfreq = 1 : NFREQ
% % Compute Wave Vector Components
% lam0 = c0/FREQ(nfreq); %free space wavelength
% k0 = 2*pi/lam0; %free space wave number
% kyinc = k0*nref; %incident wave vector
% m = [-floor(Nx/2):floor(Nx/2)]'; %spatial harmonic orders
% kx = - 2*pi*m/Sx; %wave vector expansion
% kyR = sqrt((k0*nref)^2 - kx.^2); %ky in reflection region
% kyT = sqrt((k0*ntrn)^2 - kx.^2); %ky in transmission region
% % Compute Reflectance
% ref = Eref(:,nfreq)/SRC(nfreq); %normalize to source
% ref = fftshift(fft(ref))/Nx; %compute spatial harmonics
% ref = real(kyR/kyinc) .* abs(ref).^2; %compute diffraction eff.
% REF(nfreq) = sum(ref); %compute reflectance
% % Compute Transmittance
% trn = Etrn(:,nfreq)/SRC(nfreq); %normalize to source
% trn = fftshift(fft(trn))/Nx; %compute spatial harmonics
% trn = real(kyT/kyinc) .* abs(trn).^2; %compute diffraction eff.
% TRN(nfreq) = sum(trn); %compute transmittance
% end

for nfreq = 1:NFREQ
    ref2 = abs((Eref(:,nfreq)/SRC(nfreq))).^2;
    trn2 = abs((Etrn(:,nfreq)/SRC(nfreq))).^2;
    REF2(nfreq) = sum(ref2);
    TRN2(nfreq) = sum(trn2);
end

figure(2); hold on
plot(FREQ,REF2)
plot(FREQ,TRN2)
plot(FREQ,REF2+TRN2);
hold off

% Eref0 = Eref0 + (K0^T)*Ez(:,ny_ref)*dt;
% Etrn0 = Etrn0 + (K0^T)*Ez(:,ny_trn)*dt;
% ssSRC = ssSRC + (K0^T)*Ez_src(T)*dt;

% COMPUTE ENERGY CONSERVATION
CON0 = Eref0 + Etrn0;
CON = REF + TRN;
% REPORT RESULTS AT DESIGN FREQUENCY
disp(['Reflectance = ' num2str(100*Eref0,'%4.1f') '%']);
disp(['Transmittance = ' num2str(100*Etrn0,'%4.1f') '%']);
disp(['Conservation = ' num2str(100*CON0,'%4.1f') '%']);

% INITIALIZE FIGURE WINDOW
close all;
fig = figure('Color','w');
% PLOT LINEAR REFLECTANCE, TRANSMITTANCE AND ENERGY CONSERVATION
h = plot(FREQ/gigahertz,100*REF,'-r','LineWidth',2);
hold on;
plot(FREQ/gigahertz,100*TRN,'-b','LineWidth',2);
plot(FREQ/gigahertz,100*CON,':k','LineWidth',2);
hold off;
axis([FREQ(1)/gigahertz FREQ(NFREQ)/gigahertz 0 105 ]);
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
h = legend('Reflectance','Transmittance','Conservation');
set(h,'Location','NorthEastOutside');
xlabel('Frequency (GHz)');
ylabel('%','Rotation',0,'HorizontalAlignment','right');