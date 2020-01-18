clear all;
clc;

L0 = 1e-6;  % length unit: um
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

xrange = [0 50];  % x boundaries in L0
yrange = [0 8];  % y boundaries in L0
N = [300 300];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

input_wavelength = 1.55023; % Input wavelength
target_wavelength = 1.50; % Target transition wavelength

hx = diff(xrange) / (N(1)); 
hy = diff(yrange) / (N(2)); 

eps_structure = ones(N);

%% Compute grid

% COMPUTE GRID AXES
xa = [0:N(1)-1]*hx;
ya = [0:N(2)-1]*hy;

%% define waveguide
wg_upper = 2.34; 
wg_lower = 2.06; 
eps_wg = 4;
%function handler for waveguide
within_wg = @(x, y) y > wg_lower & y < wg_upper & x<8; 
eps_structure = assign_val(eps_structure, xrange, yrange, within_wg, eps_wg);

%Display structure
% Display the structure
figure; visabs(eps_structure, xrange, yrange); 
xlabel('x (\mum)'); ylabel('y (\mum)'); 

%% Field arrays
% INITIALIZE FIELDS
Hx = zeros(N);
Hy = zeros(N);
Dz = zeros(N);
Ez = zeros(N);
% INITIALIZE CURL
CEx = zeros(N);
CEy = zeros(N);
CHz = zeros(N);
% INITIALIZE INTEGRATION TERMS
ICEx = zeros(N);
ICEy = zeros(N);
IDz = zeros(N);

Nx = N(1);
Ny = N(2);

%% Main loop
for T=1:simulation_steps
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
    ezsrc = ramp(T) * real(Ez_src*exp(-1i*2*pi*frequency*T*dt));
    CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - ezsrc/dy;

    % Update Integration Terms
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;

    % Update Hx and Hy
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;

    % Compute CHz
    CHz(1,1) = (Hy(1,1) - Hy(Nx,1))/dx - (Hx(1,1) - Hx(1,Ny))/dy;
    for nx = 2 : Nx
        CHz(nx,1) = (Hy(nx,1) - Hy(nx-1,1))/dx  - (Hx(nx,1) - Hx(nx,Ny))/dy;
    end
    for ny = 2 : Ny
        CHz(1,ny) = (Hy(1,ny) - Hy(Nx,ny))/dx   - (Hx(1,ny) - Hx(1,ny-1))/dy;
            for nx = 2 : Nx
                CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1,ny))/dx  - (Hx(nx,ny) - Hx(nx,ny-1))/dy;
            end
    end

    % TF/SF Correction
    hxsrc = ramp(T) * real(Hx_src*exp(-1i*2*pi*frequency*(T*dt + delt)));
    CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - hxsrc/dy;

    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;

    % Draw field
    imagesc(Ez);
    axis equal tight;
    axis image;
%     caxis([-emax emax]);
    %colormap('b2r');
    colorbar;
    title([num2str(T) ' of ' num2str(simulation_steps) ]);
    drawnow;
    % end

    %% Modulate the structure for next loop
    ERzz(93:132,150:900) = n_waveguide^2 + 1.225*cos(omega*T*dt);
    ERzz(132:171,150:900) = n_waveguide^2 + 1.225*cos(omega*T*dt+pi/2);


% COMPUTE EZ UPDATE COEFFICIENT
mEz1 = 1./ERzz;

end