
clear all; close all; clc;

Nsb = 2;  % Total number of frequencies is 2*Nsb+1
mode = 0; 
i = 1;
for h = 0.02:0.01:1

%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: nm

xrange = [0 10];  % x boundaries in L0
yrange = [-2.5 2.5];  % y boundaries in L0
% N = [151 51];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

wvlen0 = 1.55; % Input wavelength
wvlen1 = 1.50; % Target transition wavelength


hx = h;
hy = h;
N = [ceil(diff(xrange)/hx) ceil(diff(yrange)/hy)];

src_max1 = 1i;
src_max2 = 1; 

%%
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

omega0 = 2*pi*c0 / wvlen0; 

Omega = 2*pi*c0* (1/wvlen1 - 1/wvlen0); 

wvlen_m1 = 2*pi*c0/(omega0 - Omega); 

%% Set up the permittivity.
eps_wg = 4; 
eps_clad = 1; 

eps_space = eps_clad*ones(N);



%% Waveguide dimensions
wg_upper = 0.5001; 
wg_lower = -wg_upper; 

within_wg = @(x, y) y > wg_lower & y < wg_upper; 

eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);

%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

delta_max = 0.1; % Relative permittivity modulation strength
delta = 0.1*eps0; 

mod_x = [1.5, 9.001]; % Length of the modulation

within_mod1 = @(x, y) y > 0 & y < wg_upper & x > mod_x(1) & x < mod_x(2); 
within_mod2 = @(x, y) y > wg_lower & y < 0 & x > mod_x(1) & x < mod_x(2); 

mod_reg = assign_space(mod_reg, xrange, yrange, within_mod1, delta); 
mod_reg = assign_space(mod_reg, xrange, yrange, within_mod2, -delta); 


%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 1; 


% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, wg_upper, eps_clad, eps_wg, mode, hy, pts); 

% Location of the source
src_x = 1; 
src_y = -0.0001; 

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = 6.1164; 
Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0/sqrt(norm_P); 

% Display the structure

%% Solve for field distributions
% solver

[A,b] = condition(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_space, mod_reg, mod_phi, Jz0, Npml);
tic
x =A\b;
time(i) = toc;
res(i) = h;
i = i+1
end
