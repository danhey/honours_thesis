%% Initialise MATLAB
clc;
clear all;
clear variables;

% Mode of input wave (1 is fundamental)
mode = 1;

% Define constants
L0 = 1e-6;                  %Normalisation constant (micrometers)
e0 = 8.85418782e-12 * L0;   %Permittivity (F/L0)
u0 = 1.25663706e-6 * L0;    %Permeability ()
c0 = sqrt(1/(e0*u0));       %Speed of light (L0/s)

%% Source position (L0 units)

src_x = (1:2);
src_y = (1:2);

%% Define grid

xrange = [0 10];             % x boundaries in L0
yrange = [0 10];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
dx = 0.01;                  % Resolution along x
dy = 0.01;                  % Resolution along y

% Calculate matrix size based on given resolution
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];

%% Define problem

% SOURCE
wavelength = (1/0.129);     % Initial wavelength
wavelength_target = (1/0.198);

Omega = 2*pi*c0* (1/wavelength_target - 1/wavelength); 

frequency = c0/wavelength;
wavenumber = 2*pi/wavelength;


%% Compute grid

% % DEFAULT RESOLUTION
nmax = 3.5;
NRES = 10;
NPML = [20 20 20 20];

%% Build device
                      
% INITIALIZE MATERIALS TO FREE SPACE
mu_space_xx = ones(N);
mu_space_yy = ones(N);
eps_space = ones(N);

%% Set up the permittivity.
eps_wg = 12.25; 
eps_clad = 1;

%% Waveguide dimensions
wg_upper = 3.05;
wg_lower = 1.95;

within_wg = @(x, y) x > wg_lower & x < wg_upper;
central_guide = @(x, y) y>24.2 & y<35.8 & x >1.5 & x<3.5;

eps_space = assign_val(eps_space, xrange, yrange, within_wg, eps_wg);
eps_space = assign_val(eps_space, xrange, yrange, central_guide, eps_wg);

%% Compute source
% Courant factor time step
dmin = min([dx dy]);
dt = dmin/(2*c0);
Nt = (1/frequency)/dt;

% COMPUTE SOURCE POSITION
nx1_src = NPML(1) + 1;
nx2_src = N(1) - NPML(2);
ny_src = NPML(3) + 2;

% EXTRACT MATERIALS ACROSS INPUT SLAB WAVEGUIDE
urxx = mu_space_xx(nx1_src:nx2_src,ny_src);
uryy = mu_space_yy(nx1_src:nx2_src,ny_src);
erzz = eps_space(nx1_src:nx2_src,ny_src);

% ANALYZE WAVEGUIDE
[Ez_src,Hx_src,neff,~,~] = ezmode(urxx,uryy,erzz,wavenumber*dx, mode);
emax = max(abs(Ez_src));

% COMPUTE NUMBER OF TIME STEPS
d = sqrt(N(1)^2+N(2)^2);
tprop = nmax*d/c0;
tsim = 2*tprop;
STEPS = ceil(tsim/dt);

% COMPUTE RAMP FUNCTION
delt = 0.5*neff*dy/c0 + dt/2;


%% Compute PML
% 2x grid formulation
Nx2 = 2*N(1);
Ny2 = 2*N(2);

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
mHx2 = - c0./mu_space_xx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./mu_space_xx ./ mHx0;

% Hy coefficients
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./mu_space_yy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./mu_space_yy ./ mHy0;

% Dz coefficients
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

% Ez coefficient
mEz1 = 1./eps_space;

%% Pre-allocate variable memory
% Electromagnetic fields (TE)
Hx = zeros(N); Hy = zeros(N);
Dz = zeros(N); Ez = zeros(N);

% Curls
CEx = zeros(N); CEy = zeros(N); CHz = zeros(N);

% Integration
ICEx = zeros(N); ICEy = zeros(N); IDz = zeros(N);

%% Pre-allocate modulation space handles
% Defines a function to create arbitrary shapes on epsilon space
modLeft = @(x, y) y> 5.2 & y<24.2 & x > 2.5 & x < wg_upper;     %rectangle
modRight = @(x, y) y> 35.8 & y<54.8 & x > 2.5 & x < wg_upper;   %rectangle

%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);

%% FOURIER stuff?

%% Calculation loop
for T=1:12000
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    % Ez source  
    ezsrc =  real(Ez_src*exp(-f*T*dt));
    CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - ezsrc/dy;

    % Integration terms
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    
    % Update Hx and Hy
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;

    % Curl Hz
    CHz(1,1) = (Hy(1,1) - Hy(N(1),1))/dx - (Hx(1,1) - Hx(1,N(2)))/dy;
    CHz(2:Nx,1) = (Hy(2:Nx,1) - Hy(1:Nx-1,1))./dx  - (Hx(2:Nx,1) - Hx(2:Nx,Ny))./dy;    
    CHz(1,2:Ny) = (Hy(1,2:Ny) - Hy(Nx,2:Ny))./dx - (Hx(1,2:Ny) - Hx(1,1:Ny-1))./dy;
    CHz(2:Nx,2:Ny) = (Hy(2:Nx,2:Ny) - Hy(1:Nx-1,2:Ny))./dx - (Hx(2:Nx,2:Ny) - Hx(2:Nx,1:(Ny-1)))./dy;

    % Hx source
    hxsrc = real(Hx_src*exp(-f*(T*dt)));
    CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;
    
    
    %% Draw field
    % Update field graph every 20 steps. Comment this section out when
    % benchmarking.
    if (mod(T,20) == 0)
        visreal(Ez,xrange,yrange);
        drawnow();
    end
    
    %% Find steady state time
    aTest(T,:) = Ez(50,50);
end

L = length(time);
tr = linspace(min(time),max(time),L);
vr = resample(aTest,tr);
Ts = tr(2)-tr(1);
Fs = 1/Ts;                                          % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
FTvr = fft(vr)/L;                                   % Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
Iv = 1:length(Fv);

figure(2);
plot(Fv, abs(FTvr(Iv))*2)