function [Ez, Hx, Hy] = solve_time(L0, wavelength, Omega, xrange, yrange, src ,,
constants;	% Import constants

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

%% Compute PML (adapted from code provided by R. Rumpf)
Nx2 = 2*N(1);
Ny2 = 2*N(2);
% sigma x
sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
% sigma y
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

% Ez coefficient (Dz and Ez are calculated separately since when we modify
% eps_space we must also re-calculate its update coefficient. Arranging the
% code likes this allows for a relatively inexpensive piecewise division
% instead of re-calculating all of Dz coefficients)
mEz1 = 1./eps_space;

%% Pre-allocate variable memory
% Electromagnetic fields (TE)
Hx = zeros(N); Hy = zeros(N);
Dz = zeros(N); Ez = zeros(N);

% Curls
CEx = zeros(N); CEy = zeros(N); CHz = zeros(N);

% Integration
ICEx = zeros(N); ICEy = zeros(N); IDz = zeros(N);

%% Pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);
mod_array_length = length(modulation_array);

% Setup Fourier kernel
Fs = 1/dt;
Nyquist = Fs/2;
NFREQ = 5000;
FREQ  = linspace(0,0.5e15,NFREQ);  
K     = exp(-i*2*pi*dt.*FREQ);
REF = zeros(Nx,NFREQ);

%% Calculation loop. Field updates are fully vectorised for speed
for T=1:STEPS
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    % Ez source  
    ezsrc =  real(Ez_src*exp(-f*(T-1)*dt));
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
    hxsrc = real(Hx_src*exp(-f*((T-1)*dt)));
    CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

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
    
    % Update Ez coefficient (as discussed above)
    mEz1 = 1./eps_space;    
    
    % Update field graph every 20 steps. Drawing graphics is expensive, so
    % for benchmarking this should be commented out.
    if (mod(T,50) == 0)
        visreal(Ez,xrange,yrange);
        drawnow();
    end
    
    % Calculate Fourier transform
    for i=1:NFREQ
       REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,58/dx)*dt;
    end
end

% Complete Fourier transform and obtain power spectrum
REF2 = zeros(NFREQ,1);
for i = 1:NFREQ
    REF2(i) = sum(abs(REF(:,i)).^2);
end