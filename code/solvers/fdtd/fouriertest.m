%% Initialise MATLAB
clc; clear all;

% Mode of input wave (1 is fundamental)
mode = 1;

% Import constants
constants

%% Source position (L0 units)


%% Define grid

xrange = [0 5];             % x boundaries in L0
yrange = [0 5];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
NPML = [40 40 40 40];
dx = 0.02;
dy = 0.02;
src_x = ceil(2.5/dx);
src_y = ceil(2.5/dy);
% Calculate matrix size based on given resolution
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];

%% Define problem

% Source characteristics
wavelength = (1);     % Initial wavelength
wavelength_target = (1/0.1985);

% Modulation frequency
Omega = 2*pi*c0* (1/wavelength_target - 1/wavelength); 

frequency = c0/wavelength;
frequency2 = c0/wavelength_target;
omeg = 2*pi/wavelength;


%% Build device
                      
% Pre-allocate material space
mu_space_xx = ones(N);
mu_space_yy = ones(N);
eps_space = ones(N);

%% Define material permittivity
eps_wg = 12.25;
eps_clad = 1;

%% Waveguide dimensions
wg_upper = 3.05;
wg_lower = 1.95;

within_wg = @(x, y) x > wg_lower & x < wg_upper;

%eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);


%% Modulation space
mod_x = [2, 21];
%modLeft = @(x, y)x > 2.5 & x < wg_upper & y > mod_x(1) & y < mod_x(2);
%modRight = @(x, y) x > 1.95 & x < 2.5 & y > mod_x(1) & y < mod_x(2);

modulation_array = {};
phase_array = [0,pi/2];

%% Compute source
% Courant factor time step
dmin = min([dx dy]);
dt = dmin/(2*c0);

% Source position (legacy)
nx1_src = NPML(1) + 1;
nx2_src = N(1) - NPML(2);
ny_src = NPML(3) + 2;

% Extract 1D slab for analysis
urxx = mu_space_xx(nx1_src:nx2_src,ny_src);
uryy = mu_space_yy(nx1_src:nx2_src,ny_src);
erzz = eps_space(nx1_src:nx2_src,ny_src);

% Calculate modes
[Ez_src,Hx_src,neff,~,~] = ezmode(urxx,uryy,erzz,omeg*dx, mode);
emax = max(abs(Ez_src));

% Time steps (legacy)
d = sqrt(N(1)^2+N(2)^2);
tprop = 12.25*d/c0;
tsim = 2*tprop;
STEPS = ceil(tsim/dt);


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
t0 = 25.0;
tau = 8.0;


%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);
mod_array_length = length(modulation_array);
N_lambda=c0/(frequency);
% FOURIER
Fs = 1/dt;
Nyquist = Fs/2; %Nyquist frequency
NFREQ = 5000;
FREQ  = linspace(0,3.5e15,NFREQ);  
K     = exp(-i*2*pi*dt.*FREQ);
REF = zeros(Nx,NFREQ);
%% Calculation loop
for T=1:2000
    %source(T) = exp(-((t0-(T-1))/tau).^2);
    source = sin(((2*pi*(c0/(N_lambda))*(T)*dt)));
    Ez(src_x,src_y) = source;
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    % Ez source  
%     ezsrc =  real(Ez_src*exp(-f*(T-1)*dt));
%     CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - ezsrc/dy;

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
%     hxsrc = real(Hx_src*exp(-f*((T-1)*dt)));
%     CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;
    
    
    
    %% Modulate permittivity
%     for i = 1:mod_array_length
%         eps_space = assign_space(eps_space, xrange, yrange, modulation_array{i}, (eps_wg + 0.1*eps_wg*cos(Omega * T * dt+phase_array(i)))); 
%     end
%     mEz1 = 1./eps_space;    % Update Ez coefficient since epsilon changes
    
    %% Draw field
    % Update field graph every 20 steps. Ironically, graphing the field is
    % the slowest part of the code.
    if (mod(T,20) == 0)
        visreal(Ez,xrange,yrange);
        %contour(Ez);
       % plot(REF);
       caxis([-0.5 0.5]);
        drawnow();
    end
    
    for i=1:NFREQ
      % REF(:,i) = REF(:,i) + (K(i)^(T)).*Ez(20:120,554);
       REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,209)*dt;
    end
   % aTest(T) = Ez(62,554);
   % time(T) = T*dt;
end

%REF = REF .* (dt);

for i = 1:NFREQ
    REF2(i) = sum(abs(REF(:,i)).^2);
end

figure(2);
plot(FREQ,REF2);
% L = length(time);
% tr = linspace(min(time),max(time),L);
% vr = resample(aTest,tr);
% Ts = tr(2)-tr(1);
% Fs = 1/Ts;                                          % Sampling Frequency
% Fn = Fs/2;                                          % Nyquist Frequency
% FTvr = fft(vr)/L;                                   % Fourier Transform
% Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
% Iv = 1:length(Fv);
% 
% figure(2);
% plot(Fv, abs(FTvr(Iv))*2)
