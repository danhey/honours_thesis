%% Initialise MATLAB
clc; clear all;

% Mode of input wave (1 is fundamental)

% Import constants
L0 = 1e-6;
e0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
u0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(e0*u0);  % speed of light in vacuum in L0/sec


%% Define grid

xrange = [-2.5 2.5];             % x boundaries in L0
yrange = [0 10];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
NPML = [20 20 20 20];
dx = 0.025;                  % Resolution along x
dy = 0.025;                  % Resolution along y

% Calculate matrix size based on given resolution
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];

%% Define problem
%% Source position (legacy)
nx1_src = 100-40;
nx2_src = 100+40;
ny_src = 1/dy;

% Source characteristics
wavelength = 1.55;     % Initial wavelength
%wavelength_target = (1/0.1985);
%wavelength = (1/0.129031);
wavelength_target = 1.5;

% Modulation frequency
Omega = 2*pi*c0* (1/wavelength_target - 1/wavelength); 
omegaf = Omega/(2*pi);
frequency = c0/wavelength;
frequency2 = c0/wavelength_target;
frequency3 = c0/wavelength + 2*omegaf;
frequency4 = frequency - omegaf;
frequency5 = frequency - 2*omegaf;
omeg = 2*pi/wavelength;


%% Build device
                      
% Pre-allocate material space
mu_space_xx = ones(N);
mu_space_yy = ones(N);

%% Define material permittivity
eps_wg = 4; 
eps_clad = 1; 

eps_space = eps_clad*ones(N);



%% Waveguide dimensions
wg_upper = 0.5001; 
wg_lower = -wg_upper; 

within_wg = @(x, y) x > wg_lower & x < wg_upper; 

eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);


%% Modulation space
delta_max = 0.1; % Relative permittivity modulation strength
delta = 0.1*e0; 

mod_x = [1.5, 9.001]; % Length of the modulation

modLeft = @(x, y) x > 0 & x < wg_upper & y > mod_x(1) & y < mod_x(2); 
modRight = @(x, y) x > wg_lower & x < 0 & y > mod_x(1) & y < mod_x(2); 

modulation_array = {modLeft,modRight};
phase_array = [0,pi/2];

%% Compute source
% Courant factor time step
dmin = min([dx dy]);
dt = dmin/(2*c0);

% Calculate modes
[~, ~, ~, Ez_src, ~, ~] = beta_sym(L0, wavelength, (1/2), 1, eps_wg, 0, dx, 40); 
emax = max(abs(Ez_src));



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



%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);
mod_array_length = length(modulation_array);

% Fs = 1/dt;
% Nyquist = Fs/2; %Nyquist frequency
% NFREQ = 5000;
% FREQ  = linspace(0,0.5e15,NFREQ);  
% K     = exp(-i*2*pi*dt.*FREQ);
% REF = zeros(Nx,NFREQ);

K_1 = exp(-i*2*pi*dt.*frequency);
K_2 = exp(-i*2*pi*dt.*frequency2);
K_3 = exp(-i*2*pi*dt.*frequency3);
K_4 = exp(-i*2*pi*dt.*frequency4);
K_5 = exp(-i*2*pi*dt.*frequency5);

mode_1 = zeros(Nx,Ny);
mode_2 = zeros(Nx,Ny);
mode_3 = zeros(Nx,Ny);
mode_4 = zeros(Nx,Ny);
mode_5 = zeros(Nx,Ny);

STEPS = 50000;
%% Calculation loop
for T=1:STEPS
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
    %CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) + ezsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;
    
    %% Modulate permittivity
    for i = 1:mod_array_length
        eps_space = assign_space(eps_space, xrange, yrange, modulation_array{i}, (eps_wg + 0.1*e0*cos(Omega * T * dt+phase_array(i)))); 
    end
    mEz1 = 1./eps_space;    % Update Ez coefficient since epsilon changes
    
    %% Draw field
    % Update field graph every 20 steps. Ironically, graphing the field is
    % the slowest part of the code.
    if (mod(T,100) == 0)
        visreal(Ez,xrange,yrange);
        title(T);
        %contour(Ez);
        drawnow();
    end
    
%     for i=1:NFREQ
%       % REF(:,i) = REF(:,i) + (K(i)^(T)).*Ez(20:120,554);
%        REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,58/dx)*dt;
%     end
    for i = 1:N(2)
        mode_1(:,i) = mode_1(:,i) + (K_1^(T)).*Ez(:,i)*dt;
        mode_2(:,i) = mode_2(:,i) + (K_2^(T)).*Ez(:,i)*dt;
        mode_3(:,i) = mode_3(:,i) + (K_3^(T)).*Ez(:,i)*dt;
        mode_4(:,i) = mode_4(:,i) + (K_4^(T)).*Ez(:,i)*dt;
        mode_5(:,i) = mode_5(:,i) + (K_5^(T)).*Ez(:,i)*dt;
    end
    
end

for i = 1:N(2)
    mode_1_test(i) = sum(abs(mode_1(:,i)).^2);
    mode_2_test(i) = sum(abs(mode_2(:,i)).^2);
end

plot(1:N(2),mode_1_test);
hold on
plot(1:N(2),mode_2_test);
hold off
% for i = 1:NFREQ
%     REF2(i) = sum(abs(REF(:,i)).^2);
% end

