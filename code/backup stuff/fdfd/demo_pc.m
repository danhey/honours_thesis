% demo_pc
%
% This MATLAB code demonstrates the finite-difference frequency-domain
% method by modeling transmission and reflectance from a photonic crystal.
%
% Raymond C. Rumpf, Ph.D.
% Associate Professor of Electrical and Computer Engineering
% University of Texas at El Paso
% EL Paso, TX 79968
%
% Example code written for short course on
% "Introduction to Optical Simulation Using the Finite-Difference
% Frequency-Domain Method."

% INITIALIZE MATLAB
close all; clc;
clear all;

% UNITS
micrometers = 1;
nanometers  = micrometers / 1000;
radians     = 1;
degrees     = pi/180;

% INITIALIZE FIGURE WINDOW
ss = get(0,'ScreenSize');
figure('Position',[1 0.056*ss(4) ss(3) 0.87*ss(4)]);

%% DEFINE SIMULATION PARAMETERS

% SOURCE
LAMBDA = [500:10:1500] * nanometers;    %wavelength range to simulate
theta  = 45 * degrees;                  %angle of incidence
pol    = 'H';                           %polarization

% PHOTONIC CRYSTAL PARAMETERS
a   = 0.125 * micrometers;              %lattice constant
r   = 0.35*a;                           %hole radius
nL  = 1.0;                              %low refractive index
nH  = 3.5;                              %high refractive index
NP1 = 2;                                %lattice periods before PC
NP2 = 10;                               %lattice periods of PC
NP3 = 2;                                %lattice periods after PC
NPx = 10;                               %number of lattice perids to draw

% GRID
nmax    = max([nL nH]);                 %maximum refractive index
NRES    = 20;                           %grid resolution parameter
NPML    = 20;                           %size of perfectly matched layer
buf_ylo = 1 * micrometers;              %space between PC and PML
buf_yhi = 1 * micrometers;              %space between PC and PML

%% COMPUTE OPTIMIZED GRID

% INITIAL GRID RESOLUTION
dx = min([LAMBDA]) /nmax/NRES;
dy = dx;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);        dx = a/nx;    %snap x-grid to lattice constant
ny = ceil(a/dy);        dy = a/ny;    %snap y-grid to lattice constant

% GRID SIZE
NP = NP1 + NP2 + NP3;
Sx = a;                                 %physical size along x
Sy = buf_ylo + NP*a + buf_yhi;          %physical size along y
Nx = round(Sx/dx);                      %grid size along x
Ny = round(Sy/dy) + 2*NPML;             %grid size along y

% ENSURE Nx IS ODD
Nx = 2*round(Nx/2) + 1;
dx = Sx/Nx;

% 2X GRID PARAMETERS
Nx2 = 2*Nx;                             %grid size is twice
Ny2 = 2*Ny;     

dx2 = dx/2;                             %grid spacing is halved
dy2 = dy/2;

% GRID AXES
xa  = [0:Nx-1]*dx;      xa  = xa - mean(xa);
ya  = [0:Ny-1]*dy;
xa2 = [0:Nx2-1]*dx2;    xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2;

% REPORT GRID SIZE
disp(['Grid size is:  Nx = ' num2str(Nx)]);
disp(['               Ny = ' num2str(Ny)]);

% CLEAR TEMPORARY VARIABLES
clear dx1 dy1 nx ny;

%% BUILD PHOTONIC CRYSTAL ON 2X GRID

% CONSTRUCT UNIT CELL
    % Initialize
    Ny2_uc = round(a/dy2);
    UC     = zeros(Nx2,Ny2_uc);
    % Grid
    ya2_uc = [0:Ny2_uc-1]*dy2;
    ya2_uc = ya2_uc - mean(ya2_uc);
    [Y,X]  = meshgrid(ya2_uc,xa2);
    % Add Hole
    UC = (X.^2 + Y.^2) <= r^2;
    % Convert to Dielectric Constant
    UC = nH*(UC==0) + nL*(UC==1);
    
% BUILD PHOTONIC CRYSTAL
N2X = nL*ones(Nx2,Ny2);
ny1 = 2*NPML + round(buf_ylo/dy2);          %index of top of lattice
for np = NP1+1 : NP1+NP2
    nya = ny1 + (np-1)*Ny2_uc;
    nyb = nya + Ny2_uc - 1;
    N2X(:,nya:nyb) = UC;
end
N2X(:,nyb+1:Ny2) = nL;

% SHOW DEVICE
xa3 = [0:Nx2*NPx-1]*dx2;
P = zeros(Nx2*NPx,Ny2);
for npx = 1 : NPx
    nx1 = (npx-1)*Nx2 + 1;
    nx2 = nx1 + Nx2 - 1;
    P(nx1:nx2,:) = N2X;
end
subplot(141);
imagesc(xa3/micrometers,ya2/micrometers,P');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('DEVICE');
colorbar;
axis equal tight;
drawnow;

% COMPUTE DIELECTRIC AND MAGNETIC FUNCTIONS
ER2 = N2X.^2;                           %er = n^2
UR2 = ones(Nx2,Ny2);                    %no magnetic response

% CLEAR TEMPORARY VARIABLES
clear Ny2_uc UC ya2_uc X Y R N2X ny1 nya nyb;

%% SIMULATE PHOTONIC CRYSTAL

% INITIALIZE RECORD VARIABLES
NLAM = length(LAMBDA);                  %determine how many simulations
REF  = zeros(NLAM,1);                   %initialize reflection record
TRN  = zeros(NLAM,1);                   %initialize transmission record

% CALCULATE AXIS VECTOR
xa3 = [0:Nx*NPx-1]*dx;

% ITERATE OVER WAVELENGTH
dT = 0;
t0 = clock;
for nlam = 1 : NLAM
    
    % Get Next Wavelength
    lam0 = LAMBDA(nlam);
    
    % Compute Source
    k0   = 2*pi/lam0;
    kinc = k0*nL*[sin(theta);cos(theta)];
    
    % Call FDFD2D
    RES2 = [dx2 dy2];
    [R,T,m,F] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,pol);
    
    % Record Transmission and Reflection Response (%)
    REF(nlam) = 100*sum(R(:));          %add all spatial harmonics
    TRN(nlam) = 100*sum(T(:));          %   to compute total power
    
    % Update Graphics
    if etime(clock,t0)>dT
        % compute field over NPx unit cells
        k0    = 2*pi/lam0;
        kxinc = k0*nL*sin(theta);
        P     = zeros(Nx*NPx,Ny);
        for npx = 1 : NPx
            nx1 = (npx-1)*Nx + 1;
            nx2 = nx1 + Nx - 1;
            P(nx1:nx2,:) = F * exp(-i*kxinc*Nx*dx*(npx-1));
        end
        % show field
        subplot(142);
        imagesc(xa3/micrometers,ya/micrometers,real(P)');
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        title('FIELD');
        colorbar;
        axis equal tight;
        % show spectra
        subplot(1,4,3:4);
        plot(LAMBDA(1:nlam),REF(1:nlam),'-r'); hold on;
        plot(LAMBDA(1:nlam),TRN(1:nlam),'-b');
        plot(LAMBDA(1:nlam),REF(1:nlam)+TRN(1:nlam),':k'); hold off;
        axis([min(LAMBDA) max(LAMBDA) 0 105]);
        title('SPECTRAL RESPONSE');
        xlabel('Wavelength (\mum)');
        ylabel('%   ','Rotation',0);
        % draw now!
        drawnow;
        % update time
        t0 = clock;
    end
end





