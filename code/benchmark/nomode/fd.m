clear all; close all; clc;
i=1;
for n = 1:3
    i = 1;
for h = 0.02:0.01:1
Nsb = 0;  % Total number of frequencies is 2*Nsb+1
src_x = 1; 
%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: um

xrange = [0 5];  % x boundaries in L0
yrange = [0 5];  % y boundaries in L0
Npml = [1 1];  % [Nx_pml, Ny_pml]

hx = h;
hy = h;

N = [ceil(diff(xrange)/hx) ceil(diff(yrange)/hy)];

% 
% hx = diff(xrange) / (N(1)); 
% hy = diff(yrange) / (N(2)); 

src_max1 = 1i;
src_max2 = 1; 

%%
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec
wvlen0 = 1/0.129; % Input wavelength
wvlen1 = 1/0.1985; % Target transition wavelength
omega0 = 2*pi*c0 / wvlen0; 

Omega = 2*pi*c0* (1/wvlen1 - 1/wvlen0); 
%Omega = 1.3092e+14;

wvlen_m1 = 2*pi*c0/(omega0 - Omega); 

%% Set up the permittivity.
eps_wg = 12.25; 
eps_clad = 1;
eps_space = ones(N);



%% Waveguide dimensions
wg_upper = 3.05;
wg_lower = 1.95;

within_wg = @(x, y) y > wg_lower & y < wg_upper;
central_guide = @(x, y) x>24.2 & x<35.8 & y >1.5 & y<3.5;

eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);
eps_space = assign_space(eps_space, xrange, yrange, central_guide, eps_wg);

%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

% delta_max = 1; % Relative permittivity modulation strength
% delta = .1*(eps_wg*eps0); 
% 
% modLeft = @(x, y) x> 5.2 & x<24.2 & y > 2.5 & y < wg_upper;
% modRight = @(x, y) x> 35.8 & x<54.8 & y > 2.5 & y < wg_upper;
% 
% % 
% mod_reg = assign_space(mod_reg, xrange, yrange, modLeft, delta); 
% mod_reg = assign_space(mod_reg, xrange, yrange, modRight, delta); 
% mod_phi = assign_space(mod_phi, xrange, yrange, modRight, pi/2); 

%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 1; 

mode = 0; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, wg_upper-wg_lower, 1, eps_wg, mode, hy, pts); 

% Location of the source

src_y = 2.5000; 

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = max(jz0); 
Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0./sqrt(6.114); 


%% Solve for field distributions

% solver
[A,b] = bench(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_space, mod_reg, mod_phi, Jz0, Npml);
tic
x = A\b;
time(i,n) = toc;
res(i,n) = h;
i = i+1;
end
end