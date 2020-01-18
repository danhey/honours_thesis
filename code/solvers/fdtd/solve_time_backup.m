%% Initialise
clc; clear all;
mode = 1;
constants;

%% Define grid
xrange = [0 5];                 % x boundaries in L0
yrange = [0 60];                % y boundaries in L0
Npml = [20 20];                 % [Nx_pml, Ny_pml]
NPML = [20 20 20 20];
dx = 0.04;                      % Resolution along x
dy = 0.04;                      % Resolution along y
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];

%% Define source
wavelength = (1/0.199);         % Initial wavelength
wavelength_target = (1/0.129);  % Target wavelength
Omega = 2*pi*c0* (1/(1/0.199) - 1/(1/0.129)); 
frequency = c0/wavelength;
frequency2 = c0/wavelength_target;
omeg = 2*pi/wavelength;

nx1_src = ceil(2.5/dy - 82/2);
nx2_src = ceil(2.5/dy + 82/2);
ny_src = ceil(59/dy);
%% Pre-allocate material space
mu_space_xx = ones(N);      % mu_xx
mu_space_yy = ones(N);      % mu_yy
eps_space = ones(N);        % eps_zz
eps_wg = 12.25;             % eps_wg

%% Waveguide dimensions: defines initial permittivity
wg_upper = 3.05;            % W.g. of width 1.1a
wg_lower = 1.95;
within_wg = @(x, y) x > wg_lower & x < wg_upper;
central_guide = @(x, y) y>24.2 & y<35.8 & x >1.5 & x<3.5;
eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, central_guide, eps_wg);

%% Modulation space: defines region of modulation
modLeft = @(x, y) y> 5.2 & y<24.2 & x > 2.5 & x < wg_upper;
modRight = @(x, y) y> 35.8 & y<54.8 & x > 2.5 & x < wg_upper; 
modulation_array = {modLeft,modRight};
phase_array = [0,pi/2];

%% Compute source
dmin = min([dx dy]);
dt = dmin/(2*c0);

% Calculate modes using dispersion solver
[~, ~, ~, Ez_src, ~, ~] = beta_sym(L0, wavelength, (1.1/2), 1, eps_wg, 1, dx, 41); 

%% Compute PML: (code adapted from R. Rumpf's lectures on FDTD)
Nx2 = 2*N(1); Ny2 = 2*N(2);
sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
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

%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);
mod_array_length = length(modulation_array);

Fs = 1/dt;
Nyquist = Fs/2; %Nyquist frequency
NFREQ = 5000;
FREQ  = linspace(0,4e14,NFREQ);  
%0.5e15
K     = exp(-i*2*pi*dt.*FREQ);
REF = zeros(Nx,NFREQ);
DET = zeros(Nx,NFREQ);

STEPS = 50000;
%% Calculation loop
for T=1:10000
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    % Ez source  
    ezsrc =  real(Ez_src*exp(-f*(T-1)*dt));
    CEx(nx1_src:nx2_src,ny_src) = CEx(nx1_src:nx2_src,ny_src) - ezsrc/dy;

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
    %hxsrc = real(Hx_src*exp(-f*((T-1)*dt)));
    %CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;
    
    %% Modulate permittivity
    for i = 1:mod_array_length
        eps_space = assign_space(eps_space, xrange, yrange, modulation_array{i}, (eps_wg + 0.1*eps_wg*cos(Omega * T * dt+phase_array(i)))); 
    end
    mEz1 = 1./eps_space;
    
    %% Draw field: draw Ez every 100 steps
    if (mod(T,100) == 0)
        visreal(Ez,xrange,yrange);
        drawnow();
    end

    for i=1:NFREQ
      % REF(:,i) = REF(:,i) + (K(i)^(T)).*Ez(20:120,554);
       REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,58/dy)*dt;
       DET(:,i) = DET(:,i) + (K(i)^T)*Ez(:,ny_src-1)*dt;
    end
end

for i = 1:NFREQ
    REF2(i) = sum(abs(REF(:,i)).^2);
    DET2(i) = sum(abs(DET(:,i)).^2);
end

figure(2); hold on
FREQ = 1./(c0./FREQ);
plot(FREQ,REF2./max(REF2));
plot(FREQ,DET2./max(DET2));
hold off