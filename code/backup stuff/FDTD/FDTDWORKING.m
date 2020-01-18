%% Initialise MATLAB
clc;

% Mode of input wave (1 is fundamental)
mode = 1;

% Define constants
L0 = 1e-6;                  %Normalisation constant
e0 = 8.85418782e-12 * L0;   %Permittivity (F/L0)
u0 = 1.25663706e-6 * L0;    %Permeability ()
c0 = sqrt(1/(e0*u0));       %Speed of light (L0/s)


%% Define grid
xrange = [0 5];             % x boundaries in L0
yrange = [0 60];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
dx = 0.05;
dy = 0.05;

N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];

%% Define problem

% SOURCE
wavelength = (1/0.129);
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
URxx = ones(N);
URyy = ones(N);
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


% COMPUTE SOURCE POSITION
nx1_src = NPML(1) + 1;
nx2_src = N(1) - NPML(2);
ny_src = NPML(3) + 2;

% EXTRACT MATERIALS ACROSS INPUT SLAB WAVEGUIDE
urxx = URxx(nx1_src:nx2_src,ny_src);
uryy = URyy(nx1_src:nx2_src,ny_src);
erzz = eps_space(nx1_src:nx2_src,ny_src);

% ANALYZE WAVEGUIDE
[Ez_src,Hx_src,neff,EZR,mref] = ezmode(urxx,uryy,erzz,wavenumber*dx, mode);
emax = max(abs(Ez_src));

% COMPUTE NUMBER OF TIME STEPS
d = sqrt(N(1)^2+N(2)^2);
tprop = nmax*d/c0;
tsim = 2*tprop;
STEPS = ceil(tsim/dt);

% COMPUTE RAMP FUNCTION
tau = 3/frequency;
t0 = 3*tau;
ta = [0:STEPS-1]*dt;
ramp = exp(-((ta - t0)./tau).^2);
ind = find(ta>=t0);
ramp(ind) = 1;
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

% Hx coefficients
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./URxx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./URxx ./ mHx0;

% COMPUTE HY UPDATE COEFFICIENTS
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./URyy ./ mHy0;

% COMPUTE DZ UPDATE COEFFICIENTS
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

% COMPUTE EZ UPDATE COEFFICIENT
mEz1 = 1./eps_space;

%% Pre-allocate variable memory
% Fields
Hx = zeros(N);
Hy = zeros(N);
Dz = zeros(N);
Ez = zeros(N);
% Curls
CEx = zeros(N);
CEy = zeros(N);
CHz = zeros(N);
% Integration
ICEx = zeros(N);
ICEy = zeros(N);
IDz = zeros(N);

%% Pre-allocate modulation space handlers
modLeft = @(x, y) y> 5.2 & y<24.2 & x > 2.5 & x < wg_upper;
modRight = @(x, y) y> 35.8 & y<54.8 & x > 2.5 & x < wg_upper;
figure();

f = 1i*2*pi*frequency;

Nx = N(1);
Ny = N(2);

%% Main loop
for T=1:STEPS
    %% Update curls
    % NB should vectorise this at some point
    %CEx(1:Nx,1:Ny-1) = (Ez(1:N(1),1:N(2)-1) - Ez(1:N(1),1:N(2)-1))./dy;
    for ny = 1 : N(2)-1
        for nx = 1 : N(1)
            CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy;
        end
    end
    for nx = 1 : N(1)
        CEx(nx,N(2)) = (Ez(nx,1) - Ez(nx,N(2)))/dy;
    end
    %CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;

    % Compute CEy
    for nx = 1 : N(1)-1
        for ny = 1 : N(2)
            CEy(nx,ny) = - (Ez(nx+1,ny) - Ez(nx,ny))/dx;
        end
    end
    for ny = 1 : N(2)
        CEy(N(1),ny) = - (Ez(1,ny) - Ez(N(1),ny))/dx;
    end

    %% Update Ez src
    % TF/SF Correction
    %ezsrc = ramp(T) * real(Ez_src*exp(-1i*2*pi*frequency*T*dt));
    ezsrc =  real(Ez_src*exp(-f*T*dt));
    CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - ezsrc/dy;

    %% Update Integration Terms
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    
    % Update Hx and Hy
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;

    % Compute CHz
    CHz(1,1) = (Hy(1,1) - Hy(N(1),1))/dx - (Hx(1,1) - Hx(1,N(2)))/dy;
    for nx = 2 : N(1)
        CHz(nx,1) = (Hy(nx,1) - Hy(nx-1,1))/dx ...
        - (Hx(nx,1) - Hx(nx,N(2)))/dy;
    end
    
    for ny = 2 : N(2)
        CHz(1,ny) = (Hy(1,ny) - Hy(N(1),ny))/dx ...
        - (Hx(1,ny) - Hx(1,ny-1))/dy;
        for nx = 2 : N(1)
            CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1,ny))/dx ...
            - (Hx(nx,ny) - Hx(nx,ny-1))/dy;
        end
    end

    % TF/SF Correction
    % hxsrc =  real(Hx_src*exp(-1i*2*pi*frequency*(T*dt + delt)));
    hxsrc = real(Hx_src*exp(-f*(T*dt + delt)));
    CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;

    %% Draw field
    %Draw field every 20*dt seconds
    if (mod(T,10) == 0)
%         visreal(Ez,xrange,yrange);
%         drawnow();
        imagesc(Ez);
        axis image;
        colormap('b2r');
        colorbar;
        %caxis([-emax, emax]);
        title([num2str(T*dt) 's of ' num2str(STEPS*dt) ]);
        drawnow;
    end
    % end
    
    %% Update modulation
    eps_space = assign_val(eps_space, xrange, yrange, modLeft, (12.25 + 0.1*12.25*cos(Omega * T * dt))); 
    eps_space = assign_val(eps_space, xrange, yrange, modRight, (12.25 + 0.1*12.25*cos(Omega * T * dt + pi/2))); 
    mEz1 = 1./eps_space;

end