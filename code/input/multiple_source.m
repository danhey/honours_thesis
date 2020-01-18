clear all; close all; clc;

Nsb = 0;  % Total number of frequencies is 2*Nsb+1

%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: nm

xrange = [0 12];  % x boundaries in L0
yrange = [0 12];  % y boundaries in L0
% N = [151 51];  % [Nx Ny]
N = [300 300];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

wvlen0 = 1.55023; % Input wavelength
wvlen1 = 1.50; % Target transition wavelength


hx = diff(xrange) / (N(1)); 
hy = diff(yrange) / (N(2)); 

src_max1 = 1i;
src_max2 = 1; 

%%
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

omega0 = 2*pi*c0 / wvlen0; 

Omega = 10e9;

wvlen_m1 = 2*pi*c0/(omega0 - Omega); 

%% Set up the permittivity.
eps_wg = 4; 
eps_clad = 1; 
eps_ring = complex(4,-5e-4);

eps_space = eps_clad*ones(N);



%% Waveguide dimensions
wg_upper = 2.34; 
wg_lower = 2.06; 

within_wg = @(x, y) y > wg_lower & y < wg_upper; 

ring_radius_inner = 2.9;
ring_radius_outer = 3.18;
spacing = 0.48;

within_ring = @(x, y) (x-6).^2+(y-6).^2<ring_radius_outer^2 & (x-6).^2+(y-6).^2>ring_radius_inner^2;

%eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);
%eps_space = assign_space(eps_space, xrange, yrange, within_ring, eps_ring);


%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

delta_max = 0.1; % Relative permittivity modulation strength
delta = delta_max * eps0; 

mod_x = [1.5, 9.001]; % Length of the modulation

within_ring_mod1 = @(x, y) (x-6).^2+(y-6).^2<3.04^2 & (x-6).^2+(y-6).^2>ring_radius_inner^2;
within_ring_mod2 = @(x, y) (x-6).^2+(y-6).^2<3.18^2 & (x-6).^2+(y-6).^2>3.04^2;

%mod_reg = assign_space(mod_reg, xrange, yrange, within_ring_mod1, delta); 
%mod_reg = assign_space(mod_reg, xrange, yrange, within_ring_mod2, -delta); 


%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 6; 

mode = 0; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, wg_upper, eps_clad, eps_wg, mode, hy, pts); 

% Location of the source
[Y,X] = meshgrid(yrange,xrange);
R     = sqrt(X.^2 + Y.^2);
fsrc  = exp(1i*(1/wvlen0)*1*R)./sqrt(R);
% CONSTRUCT Q
Q = R <(0.25*wvlen0);
Jz0 = diag(sparse(Q(:)));
%% Solve for field distributions

% solver
[Ez, Hx, Hy, omega] = solve_frequency(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_space, mod_reg, mod_phi, Jz0, Npml);


%% Extract modal power amplitude from simulation
a_sim = zeros(2*Nsb+1, N(1)); 
beta = zeros(2*Nsb+1, 1); 
ez = zeros(2*pts+1, 2*Nsb+1); 
amp = zeros(2*Nsb+1, 1); 
x = linspace(xrange(1), xrange(2), N(1)); 
kappa = zeros(2*Nsb, 1); 


figure; 
for i = 1 : 2*Nsb+1
    subplot(2*Nsb+1, 1, i); 
    visreal(Ez{i}, xrange, yrange); 
    
    Ez_i = Ez{i}; 
    wvlen_i = 2*pi*c0 / omega(i); 
    
    [beta(i), ~, ~, ez(:, i), amp(i)] = beta_sym(L0, wvlen_i, wg_upper, eps_clad, eps_wg, mod(i-1-Nsb, 2), hy, pts);
    
    for j = 1:N(1)
        a_sim(i, j) = max(max(abs(Ez_i(j, :)))) / amp(i); 
    end
    
end

