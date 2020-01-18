%% Initialise MATLAB
 clear all;

% Mode of input wave (1 is fundamental)
mode = 1;

% Import constants
constants

%% Source position (L0 units)
% Not yet working..
src_x = (1:2);
src_y = (1:2);

%% Define grid

xrange = [-2 2];             % x boundaries in L0
yrange = [0 14.04];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
NPML = [20 20 20 20];
dx = 0.02;                  % Resolution along x
dy = 0.02;                  % Resolution along y

% Calculate matrix size based on given resolution
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];


%nx1_src = NPML(1) + 1;
%nx2_src = N(1) - NPML(2);
nx2_src = (0.10+2)/dy;
nx1_src = (-0.10+2)/dy;
ny_src = 22;
%% Define problem

% Source characteristics
wavelength = (1/0.6468);     % Initial wavelength
wavelength_target = (1/0.8879);
k1 = 1.836;
k2 = 1.367;
% Modulation frequency
Omega = 2*pi*c0* (1/0.8879 - 1/0.6468); 
Qtarg = 2*pi*(k1-k2);
frequency = c0/wavelength;
frequency2 = c0/wavelength_target;
omeg = 2*pi/wavelength;
% Initial wavelength
%% Build device
                      
% Pre-allocate material space
mu_space_xx = ones(N);
mu_space_yy = ones(N);
eps_space = ones(N);

%% Define material permittivity
eps_wg = 12.25;
eps_clad = 1;

%% Waveguide dimensions
wg_upper = 0.11;
wg_lower = -0.11;
arm_sep = 0.58;
mod_start = 4.51;
mod_end = 9.53;
within_wg1 = @(x, y) x > wg_lower & x < wg_upper & y < 2;
within_wg2 = @(x, y) x > wg_lower & x < wg_upper & y >12.02;
within_wg3 = @(x, y) x > 0.45 & x < 0.67 & y > mod_start & y < mod_end;
within_wg4 = @(x, y) x < -0.45 & x > -0.67 & y < mod_end & y>mod_start;

arm1 = @(x,y) x < 0.224*y - 0.338 & y <4.51 & y > 2 & x > 0.224*y - 0.558;
arm2 = @(x,y) x > -0.224*y + 0.338 & y <4.51 & y > 2 & x < -0.224*y + 0.558;
arm3 = @(x,y) x < -0.224*y + 2.80248 & y <12.02 & y > 9.52 & x > -0.224*y +2.58248;
arm4 = @(x,y) x > 0.224*y - 2.80248 & y <12.02 & y > 9.52 & x < 0.224*y -2.58248;
eps_space = assign_space(eps_space, xrange, yrange, within_wg1, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, within_wg2, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, within_wg3, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, within_wg4, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, arm1, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, arm2, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, arm3, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, arm4, eps_wg);
%% Modulation space

mod1 = @(x, y) x > 0.29 & x < 0.4 & y > 7 & y < 17.04;
mod2 = @(x, y) x > 0.4 & x < 0.51 & y > 7 & y < 17.04;
mod3 = @(x, y) x < -0.29 & x > -0.4 & y > 7 & y < 17.04;
mod4 = @(x, y) x < -0.4 & x > -0.51 & y > 7 & y < 17.04;
modulation_array = {mod1,mod2,mod3,mod4};
phase_array = [0,pi,0,pi];

%% Compute source
% Courant factor time step
dmin = min([dx dy]);
dt = dmin/(2*c0);

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



%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);
mod_array_length = length(modulation_array);


% K_1 = exp(-i*2*pi*dt.*frequency);
% K_2 = exp(-i*2*pi*dt.*frequency2);
% modlength = ceil((mod_end-mod_start)/dx);
% mode_1 = zeros(Nx,modlength);
% mode_2 = zeros(Nx,modlength);
% FOURIER
NFREQ = 500;
FREQ  = linspace(2.5e14,2.7e14,NFREQ);  
%0.5e15
K     = exp(-i*2*pi*dt.*FREQ);
REF = zeros(length(ceil((diff(xrange)/dx)/2 - 0.11/dy):ceil((diff(xrange)/dx)/2 + 0.11/dy)),NFREQ);
DET = zeros(Nx,NFREQ);
SRC = zeros(1,NFREQ);
STEPS = 10000;
tau = 1/(2*frequency);
t0 = 6*tau;
t = (0:STEPS-1)*dt;
Esrc = exp(-((t-t0)/tau).^2); 
%% Calculation loop
for T=1:STEPS
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    % Ez source  
    %ezsrc =  real(Ez_src*exp(-f*(T-1)*dt));
    ezsrc = Esrc(T);
    CEx(nx1_src:nx2_src,22) = CEx(nx1_src:nx2_src,22) - ezsrc/dy;

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
    for i = (mod_start/dy:mod_end/dy)
       eps_space(ceil((0.45+0.9)/dx:(0.56+0.9)/dx),ceil(i)) = (eps_wg + 1*cos(Omega * T * dt + i*Qtarg*dx));
       eps_space(ceil((0.56+0.9)/dx:(0.67+0.9)/dx),ceil(i)) = (eps_wg + 1*cos(Omega * T * dt + i*Qtarg*dx + pi));
       eps_space(ceil((0.45+2*0.9+0.22)/dx:(0.56+2*0.9+0.22)/dx),ceil(i)) = (eps_wg + 1*cos(Omega * T * dt + i*Qtarg*dx + pi));
        eps_space(ceil((0.56+2*0.9+0.22)/dx:(0.67+2*0.9+0.22)/dx),ceil(i)) = (eps_wg + 1*cos(Omega * T * dt + i*Qtarg*dx));
    end
    mEz1 = 1./eps_space;    % Update Ez coefficient since epsilon changes
    
    %% Draw field
    % Update field graph every 20 steps. Ironically, graphing the field is
    % the slowest part of the code.
    if (mod(T,20) == 0)
        visreal(Ez,xrange,yrange);
        title(T);
        caxis([-0.1 0.1]);
        %contour(Ez);
        drawnow();
    end
    
%     %Fourier
%     for i = 1:((mod_end-mod_start)/dy)
%         mode_1(:,i) = mode_1(:,i) + (K_1^(T)).*Ez(:,i+299)*dt;
%         mode_2(:,i) = mode_2(:,i) + (K_2^(T)).*Ez(:,i+299)*dt;
%     end

    %Total fourier amplitudes of all freqs
%     for i=1:NFREQ
%       % REF(:,i) = REF(:,i) + (K(i)^(T)).*Ez(20:120,554);
%        REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,9/dy)*dt;
%        DET(:,i) = DET(:,i) + (K(i)^T)*Ez(:,ny_src+1)*dt;
%     end
    for i=1:NFREQ
       REF(:,i) = REF(:,i) + (K(i)^T)*Ez(ceil((diff(xrange)/dx)/2 - 0.11/dy):ceil((diff(xrange)/dx)/2 + 0.11/dy),N(2)-22)*dt;
       SRC(i) = SRC(i) + (K(i)^T)*ezsrc*dt;
       %DET(:,i) = DET(:,i) + (K(i)^T)*Ez(:,ny_src)*dt;
    end
end

% for i = 1:modlength
% mode_1_test(i) = sum(abs(mode_1(:,i)).^2);
% mode_2_test(i) = sum(abs(mode_2(:,i)).^2);
% end
% figure(2);
% hold on
% plot(1:modlength,mode_1_test./max(mode_1_test));
% plot(1:modlength,mode_2_test./max(mode_2_test));
% hold off

for i = 1:NFREQ
    REF2(i) = sum((REF(:,i)));
    SRC2(i) = sum((SRC(i)));
    %DET2(i) = sum(abs(DET(:,i)).^2);
end

figure(2); hold on
FREQ = 1./(c0./FREQ);
plot(FREQ,abs(REF2./SRC2).^2);
%plot(FREQ,DET2./max(DET2));
hold off

