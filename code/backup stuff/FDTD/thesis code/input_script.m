%% Input script

%Normalised dimension
a = 1e-6; 

x_range = [0,10];
y_range = [0,10];
resolution = 0.05;

pml_layer = [20,20];

%Simulation settings
sidebands = 2;
input_mode = 0;

%CONSTANTS
EPS_0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
MU_0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
C_0 = 1/sqrt(EPS_0*MU_0);  % speed of light in vacuum in L0/sec


omega0 = 2*pi*C_0 / wvlen0; 

Omega = 2*pi*C_0* (1/wvlen1 - 1/wvlen0); 

wvlen_m1 = 2*pi*C_0/(omega0 - Omega); 

%Solver
[Ez, Hx, Hy, omega] = solve_mf_fdfd(L0, input_wavelength, mod_frequency, sidebands, x_range, y_range, eps_space, mod_region, mod_phase, Jz0, pml_layer);