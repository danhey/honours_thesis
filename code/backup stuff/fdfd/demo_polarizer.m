% demo_polarizer
%
% This MATLAB code demonstrates the finite-difference frequency-domain
% method by modeling transmission and reflectance from a broadband 
% wire grid polarizer.
%
% See Phot. Technol. Lett., Vol. 22, No. 9, pp. 697-699, 2008
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
LAMBDA = [0.3:0.05:10] * micrometers;     %wavelength range to simulate
theta  = 0 * degrees;                   %angle of incidence

% BROADBAND WIRE GRID POLARIZER
L  = 80 * nanometers;                   %polarizer grating period
h  = 60 * nanometers;                   %polarizer grating height
t1 = 40 * nanometers;                   %top metal layer thickness
t2 = 30 * nanometers;                   %bottom metal layer thickness
f  = 0.5;                               %grating duty cycle

n1 = 1.00;                              %refractive index of superstrate
n2 = 1.35;                              %refractive index of substrate
nm = 2 - 1i*20;                          %refractive index of metal

% GRID
nmax    = max(real([n1 n2 nm]));        %maximum refractive index
NRES    = 20;                           %grid resolution parameter
NPML    = 20;                           %size of perfectly matched layer
buf_ylo = 0.1 * micrometers;            %space between polarizer and PML
buf_yhi = 0.1 * micrometers;            %space between polarizer and PML

%% COMPUTE OPTIMIZED GRID

% INITIAL GRID RESOLUTION
dx = min(LAMBDA) /nmax/NRES;            %must resolve smallest wavelength
dy = dx;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);        dx = L/nx;      %snap x-grid to GMR period
ny = ceil(h/dy);        dy = h/ny;      %snap y-grid to GMR thickness

% GRID SIZE
Sx = L;                                 %physical size along x
Sy = buf_ylo + h + t1 + buf_yhi;        %physical size along y
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

%% BUILD DEVICE ON 2X GRID

% COMPUTE STRUCTURE INDICES
nx1 = round((Nx2-f*L/dx2)/2);           %position of left of high index
nx2 = round((Nx2+f*L/dx2)/2);           %position of right of high index

ny1  = 2*NPML + round(buf_ylo/dy2);     %top of top metal layer
ny1  = 2*round(ny1/2) + 1;              %make position index odd
ny2  = ny1 + round(t1/dy2) - 1;         %bottom of top metal layer
ny3  = ny2 + 1;                         %top of grating
ny4  = ny3 + round(h/dy2) - 1;          %bottom of grating
ny3b = ny4 - round(t2/dy2) + 1;         %top of bottom metal layer

% BUILD POLARIZER
N2X = n1 * ones(Nx2,Ny2);               %initialize entire grid with n1
N2X(nx1:nx2,ny1:ny2)    = nm;           %top metal layer
N2X(nx1:nx2,ny3:ny4)    = n2;           %grating tooth
N2X(1:nx1-1,ny3b:ny4)   = nm;           %left bottom metal layer
N2X(nx2+1:Nx2,ny3b:ny4) = nm;           %left bottom metal layer
N2X(:,ny4+1:Ny2)        = n2;           %substrate

% SHOW DEVICE
subplot(2,3,1);
imagesc(xa/micrometers,ya/micrometers,real(N2X'));
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

%% SIMULATE DEVICE

% INITIALIZE RECORD VARIABLES
NLAM = length(LAMBDA);                  %determine how many simulations
REFE  = zeros(NLAM,1);                  %initialize reflection record
TRNE  = zeros(NLAM,1);                  %initialize transmission record
REFH  = zeros(NLAM,1);                  %initialize reflection record
TRNH  = zeros(NLAM,1);                  %initialize transmission record

% ITERATE OVER WAVELENGTH
for nlam = 1 : NLAM
    
    % Get Next Wavelength
    lam0 = LAMBDA(nlam);
    
    % Compute Source
    k0   = 2*pi/lam0;
    kinc = k0*n1*[sin(theta);cos(theta)];
    
    % Call FDFD2D
    RES2 = [dx2 dy2];
    [RE,TE,m,FE] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,'E');
    [RH,TH,m,FH] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,'H');
    
    % Record Transmission and Reflection Response (%)
    REFE(nlam) = 100*sum(RE(:));          %add all spatial harmonics
    TRNE(nlam) = 100*sum(TE(:));          %   to compute total power
    REFH(nlam) = 100*sum(RH(:));          %add all spatial harmonics
    TRNH(nlam) = 100*sum(TH(:));          %   to compute total power
    
    % Show Field
    subplot(232);
    imagesc(xa/micrometers,ya/micrometers,real(FE)');
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title('E_z FIELD');
    colorbar;
    
    subplot(233);
    imagesc(xa/micrometers,ya/micrometers,real(FH)');
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    title('H_z FIELD');
    colorbar;
    
    % Show Spectra
    subplot(2,3,4:6);
    plot(LAMBDA(1:nlam),REFE(1:nlam),'-r'); hold on;
    plot(LAMBDA(1:nlam),TRNE(1:nlam),'-b');
    plot(LAMBDA(1:nlam),REFH(1:nlam),':r');
    plot(LAMBDA(1:nlam),TRNH(1:nlam),':b'); hold off;
    axis([min(LAMBDA) max(LAMBDA) 0 105]);
    legend('E_{ref}','E_{trn}','H_{ref}','H_{trn}','Location','NorthEastOutside');
    xlabel('Wavelength (\mum)');
    ylabel('%   ','Rotation',0);
    title('SPECTRAL RESPONSE');
    
    drawnow;                            %update graphics now!
end





