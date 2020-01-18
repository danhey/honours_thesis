% Lecture22_bend.m
% INITIALIZE MATLAB
clc;
clear all;
% UNITS
meters      = 1;
micrometers = 1e-6;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
% CONSTANTS
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;
% SOURCE
lam0 = 1.55 * micrometers;
f0   = c0/lam0;
k0   = 2*pi/lam0;
% RIB WAVEGUIDE PARAMETERS
nair  = 1.0;
nsub  = 1.52;
ncore = 1.9;
w = 0.5 * lam0;         
h = 0.5 * lam0;         
rbend = 9 * micrometers;
% GRID PARAMETERS
nmax = max([nair nsub ncore]);
NRES = 10;
BUF  = 5*lam0 * [1 2 1 1];
NPML = [40 40 40 40];

% COMPUTE SLAB GRID
Sz = 7*h;               
%size of grid
dz = lam0/NRES/nmax;    
%initial resolution
Nz = ceil(h/dz);        
%snap grid to h
dz = h/Nz;              
%revised resolution
Nz = ceil(Sz/dz);       
%total number of cells
% INITIALIZE MATERIALS TO FREE SPACE
ERzz = ones(1,Nz);
URxx = ones(1,Nz);
URyy = ones(1,Nz);
% COMPUTE POSITION INDICES
nz1 = round(3*h/dz);
nz2 = nz1 + round(h/dz) - 1;
% CLADDING
neff2 = 0.3*nair + 0.7*nsub;
% CORE SLAB
ERzz = nair^2 * ones(1,Nz);
ERzz(nz1:nz2)  = ncore^2;
ERzz(nz2+1:Nz) = nsub^2;
[Ez,Hx,neff1]  = ezmode(URxx,URyy,ERzz,k0*dz);

% DEFAULT RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;
% SNAP GRID TO CRITICAL DIMENSIONS
Nx = ceil(w/dx);        dx = w/Nx;
Ny = ceil(w/dy);        dy = w/Ny;
% COMPUTE GRID SIZE
a  = rbend + w/2;
Sx = BUF(1) + a + BUF(2);
Nx = ceil(Sx/dx) + NPML(1) + NPML(2);
Sx = Nx*dx;
Sy = BUF(3) + a + BUF(4);
Ny = ceil(Sy/dy) + NPML(3) + NPML(4);
Sy = Ny*dy;
% COMPUTE GRID AXES
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;

% INITIALIZE MATERIALS TO FREE SPACE
URxx = ones(Nx,Ny);
URyy = ones(Nx,Ny);
ERzz = ones(Nx,Ny);

% COMPUTE START AND STOP INDICES OF BEND WINDOW
nx  = round(a/dx);
nx1 = NPML(1) + round(BUF(1)/dx);
nx2 = nx1 + nx - 1;
ny  = round(a/dy);
ny1 = NPML(3) + round(BUF(3)/dx);
ny2 = ny1 + ny - 1;

% COMPUTE CENTER OF CURVATURE
x0 = nx2*dx;
y0 = ny1*dy;

% CONSTRUCT BEND
[Y,X] = meshgrid(ya,xa);
R1    = ((X - x0).^2 + (Y - y0).^2) <= (rbend + w/2)^2;
R2    = ((X - x0).^2 + (Y - y0).^2) >= (rbend - w/2)^2;
R     = R1 .* R2;

% CLIP BEND WINDOW
R(1:nx1-1,:)  = 0;
R(nx2+1:Nx,:) = 0;
R(:,1:ny1-1)  = 0;
R(:,ny2+1:Ny) = 0;

% ADD INPUT WAVEGUIDE
nx  = round(w/dx);
nxa = nx1;
nxb = nxa + nx - 1;
R(nxa:nxb,1:ny1-1) = 1;

% ADD OUTPUT WAVEGUIDE
ny  = round(w/dy);
nya = ny2 - ny + 1;
nyb = ny2;
R(nx2+1:Nx,nya:nyb) = 1;

% BUILD MATERIALS
ERzz = ones(Nx,Ny);

% COMPUTE STABLE TIME STEP
dmin = min([dx dy]);
dt   = dmin/(2*c0);

% SNAP TIME STEP SO WAVE PERIOD IS AN INTEGER NUMBER OF STEPS
period = 1/f0;
Nt     = ceil(period/dt);
dt     = period/Nt;

% COMPUTE SOURCE POSITION
nx1_src = NPML(1) + 1;
nx2_src = Nx - NPML(2);
ny_src  = NPML(3) + 2;

% EXTRACT MATERIALS ACROSS INPUT SLAB WAVEGUIDE
urxx = URxx(nx1_src:nx2_src,ny_src);
uryy = URyy(nx1_src:nx2_src,ny_src);
erzz = ERzz(nx1_src:nx2_src,ny_src);

% ANALYZE WAVEGUIDE
[Ez_src,Hx_src,neff,EZR,mref] = ezmode(urxx,uryy,erzz,k0*dx);
emax = max(abs(Ez_src));

delt = 0.5*neff*dy/c0 + dt/2;

d     = sqrt(Sx^2 + Sy^2);
tprop = nmax*d/c0;
tsim  = 2*tprop;
STEPS = ceil(tsim/dt);

tau = 3/f0;
t0  = 3*tau;
ta  = [0:STEPS-1]*dt;
ramp = exp(-((ta - t0)./tau).^2);
ind  = find(ta>=t0);
ramp(ind) = 1;

% KERNEL FOR f0
K = exp(-1i*2*pi*dt*f0);

% INITIALIZE STEADY-STATE FIELDS

% POSITION OF RECORD PLANES
nyref = NPML(3) + 1;
nxtrn = Nx - NPML(2);

% EXTRACT MATERIALS
ny1 = NPML(3) + 1;
ny2 = Ny - NPML(4);
nx  = Nx - NPML(2);
urxx = URxx(nx,ny1:ny2);
uryy = URyy(nx,ny1:ny2);
erzz = ERzz(nx,ny1:ny2);

% ANALYZE OUTPUT WAVEGUIDE
[Ez,Hx,n,EZT,mtrn] = ezmode(urxx,uryy,erzz,k0*dy);

% NUMBER OF POINTS ON 2X GRID
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

% COMPTUE HX UPDATE COEFFICIENTS
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0  = (1/dt) + sigHy/(2*e0);
mHx1  = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2  = - c0./URxx./mHx0;
mHx3  = - (c0*dt/e0) * sigHx./URxx ./ mHx0;

% COMPUTE HY UPDATE COEFFICIENTS
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0  = (1/dt) + sigHx/(2*e0);
mHy1  = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2  = - c0./URyy./mHy0;
mHy3  = - (c0*dt/e0) * sigHy./URyy ./ mHy0;

% COMPUTE DZ UPDATE COEFFICIENTS
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0  = (1/dt) + (sigDx + sigDy)/(2*e0) ...
+ sigDx.*sigDy*(dt/4/e0^2);
mDz1  = (1/dt) - (sigDx + sigDy)/(2*e0) ...
- sigDx.*sigDy*(dt/4/e0^2);
mDz1  = mDz1 ./ mDz0;
mDz2  = c0./mDz0;
mDz4  = - (dt/e0^2)*sigDx.*sigDy./mDz0;

% COMPUTE EZ UPDATE COEFFICIENT
mEz1  = 1./ERzz;

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
IDz  = zeros(Nx,Ny);
aTest = zeros(1,STEPS);

% POSITION OF RECORD PLANES
nyref = NPML(3) + 1;
nxtrn = Nx - NPML(2);

NFREQ = 1000;
% INITIALIZE STEADY-STATE FIELDS
REF = zeros(Nx,NFREQ);
TRN = zeros(Nx,NFREQ);
% KERNEL FOR f0
FREQ  = linspace(0,f0,NFREQ);  
K     = exp(-i*2*pi*dt.*FREQ);
%K = exp(-1i*2*pi*dt*f0);
SRC = zeros(1,NFREQ);
t0 = 25.0;
tau = 8.0;

for T=1:STEPS
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

% TF/SF Correction
ezsrc = exp(-0.5*((t0-(T-1))/tau).^2);
CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - ezsrc/dy;

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

% TF/SF Correction
hxsrc =  real(Hx_src*exp(-1i*2*pi*f0*(T*dt + delt)));
CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

% Update Integration Term
IDz = IDz + Dz;

Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;

Ez = mEz1.*Dz;


% FINISH TRANSFORMS

% Update Graphical Status
if ~mod(T,50)    
    % draw field
    imagesc(real(Ez));

    % force MATLAB to draw graphics
    drawnow;
end

        for i = 1:NFREQ
        REF(:,i) = REF(:,i) + (K(i).^(T)).*Ez(:,nyref);
        TRN(:,i) = TRN(:,i) + (K(i).^(T)).*Ez(:,359);
        SRC(i) = SRC(i) + (K(i).^T).*ezsrc;
        end
end


ref2 = abs(REF./SRC).^2;
trn2 = abs(TRN./SRC).^2;

for nfreq = 1:NFREQ
   % ref2 = abs((REF(:,nfreq)/SRC(nfreq)));
    %trn2 = abs((TRN(:,nfreq)/SRC(nfreq)));
    REF2(nfreq) = sum(ref2(:,nfreq));
    TRN2(nfreq) = sum(trn2(:,nfreq));
end

figure(2); hold on
plot(FREQ,REF2)
plot(FREQ,TRN2)
plot(FREQ,REF2+TRN2);
legend('REF','TRN','CON');
hold off