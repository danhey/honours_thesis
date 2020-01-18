% demo_gmr
%
% This MATLAB code demonstrates the finite-difference frequency-domain
% method by modeling transmission and reflectance from a guided-mode
% resonance filter.
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
LAMBDA = [400:1:700] * nanometers;      %wavelength range to simulate
theta  = 0 * degrees;                   %angle of incidence
pol    = 'E';                           %polarization or mode

% GMR
T  = 134 * nanometers;                  %thickness of GMR grating
L  = 314 * nanometers;                  %period of GMR grating
n1 = 1.00;                              %refractive index of superstrate
n2 = 1.52;                              %refractive index of substrate
nL = 2.00;                              %low refractive index of GMR
nH = 2.10;                              %high refractive index of GMR
f  = 0.5;                               %duty cycle of GMR grating

% GRID
nmax    = max([n1 n2 nL nH]);           %maximum refractive index
NRES    = 20;                           %grid resolution parameter
NPML    = 20;                           %size of perfectly matched layer
buf_ylo = 500 * nanometers;             %space between GMR and PML
buf_yhi = 500 * nanometers;             %space between GMR and PML

%% COMPUTE OPTIMIZED GRID

% INITIAL GRID RESOLUTION
dx1 = min(LAMBDA) /nmax/NRES;           %must resolve smallest wavelength
dy1 = dx1;

dx2 = f*L/NRES;                         %must resolve finest dimension
dy2 = T/NRES;                           

dx = min(dx1,dx2);                      %choose smallest numbers
dy = min(dy1,dy2);

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);        dx = L/nx;      %snap x-grid to GMR period
ny = ceil(T/dy);        dy = T/ny;      %snap y-grid to GMR thickness

% GRID SIZE
Sx = L;                                 %physical size along x
Sy = buf_ylo + T + buf_yhi;             %physical size along y
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

%% BUILD GMR DEVICE ON 2X GRID

% COMPUTE STRUCTURE INDICES
nx1 = round((Nx2-f*L/dx2)/2);           %position of left of high index
nx2 = round((Nx2+f*L/dx2)/2);           %position of right of high index

ny1 = 2*NPML + round(buf_ylo/dy2);      %position of top of GMR
ny1 = 2*round(ny1/2) + 1;               %make position index odd
ny2 = ny1 + round(T/dy2) - 1;           %compute bottom position

% BUILD GMR
N2X = n1 * ones(Nx2,Ny2);               %initialize entire grid with n1
N2X(:,ny1:ny2)       = nL;              %fill GMR with nL
N2X(nx1:nx2,ny1:ny2) = nH;              %add nH to GMR
N2X(:,ny2+1:Ny2)     = n2;              %add n2 to superstrate region

% SHOW DEVICE
subplot(141);
imagesc(xa/micrometers,ya/micrometers,N2X');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('DEVICE');
colorbar;
drawnow;

% COMPUTE DIELECTRIC AND MAGNETIC FUNCTIONS
ER2 = N2X.^2;                           %er = n^2
UR2 = ones(Nx2,Ny2);                    %no magnetic response

% CLEAR TEMPORARY VARIABLES
clear nx1 nx2 ny1 ny2 N2X;

%% SIMULATE GUIDED-MODE RESONANCE FILTER

% INITIALIZE RECORD VARIABLES
NLAM = length(LAMBDA);                  %determine how many simulations
REF  = zeros(NLAM,1);                   %initialize reflection record
TRN  = zeros(NLAM,1);                   %initialize transmission record

% ITERATE OVER WAVELENGTH
for nlam = 1 : NLAM
    
    % Get Next Wavelength
    lam0 = LAMBDA(nlam);
    
    % Compute Source Vector
    k0   = 2*pi/lam0;
    kinc = k0*n1*[sin(theta);cos(theta)];
    
    % Call FDFD2D
    RES2 = [dx2 dy2];
    [R,T,m,F] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,pol);
    
    % Record Transmission and Reflection Response (%)
    REF(nlam) = 100*sum(R(:));          %add all spatial harmonics
    TRN(nlam) = 100*sum(T(:));          %   to compute total power
    
    % Show Field
    subplot(142);
    imagesc(xa/micrometers,ya/micrometers,real(F)');
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title('FIELD');
    colorbar;
    
    % Show Spectra
    subplot(1,4,3:4);
    plot(LAMBDA(1:nlam)/nanometers,REF(1:nlam),'-r'); hold on;
    plot(LAMBDA(1:nlam)/nanometers,TRN(1:nlam),'-b');
    plot(LAMBDA(1:nlam)/nanometers,REF(1:nlam)+TRN(1:nlam),':k'); hold off;
    axis([min(LAMBDA)/nanometers max(LAMBDA)/nanometers 0 105]);
    xlabel('Wavelength (nm)');
    ylabel('%   ','Rotation',0);
    title('SPECTRAL RESPONSE');
    
    drawnow;                            %update graphics now!
end





