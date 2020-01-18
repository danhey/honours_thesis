clear all; close all; clc;

Nsb = 3;  % Total number of frequencies is 2*Nsb+1
ring_radius_inner = 2.9;
ring_radius_outer = 3.18;
spacing = 0.48;
number_rings = 1;
ring_middle_y = 5;

%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: nm

xrange = [0 12];  % x boundaries in L0
yrange = [0 10];  % y boundaries in L0
% N = [151 51];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

wvlen0 = 1.55; % Input wavelength
wvlen1 = 1.50; % Target transition wavelength

hx = 0.05;
hy = 0.05;

N = [ceil(diff(xrange)/hx) ceil(diff(yrange)/hy)];

src_max1 = 1i;
src_max2 = 1; 

%%
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

omega0 = 2*pi*c0 / wvlen0; 
n_0 = sqrt(4);
Omega = c0/(n_0 * ring_radius_outer);

wvlen_m1 = 2*pi*c0/(omega0 - Omega); 

%% Set up the permittivity.
eps_wg = n_0^2; 
eps_clad = 1; 
eps_ring = n_0^2;

eps_space = eps_clad*ones(N);



%% Waveguide dimensions
wg_left = 2.34; 
wg_right = 2.06; 

within_wg = @(x, y) x > wg_right & x < wg_left; 

for i = 0:number_rings-1
    within_ring = @(x,y) (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2<ring_radius_outer^2 & (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2>ring_radius_inner^2;
    eps_space = assign_space(eps_space, xrange, yrange, within_ring, eps_ring);
end

eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);

%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

delta_max = 0.1; % Relative permittivity modulation strength
delta = delta_max * eps0 * eps_ring; 


for i = 0:number_rings-1
      within_ring = @(x,y) (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2<ring_radius_outer^2 & (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2>ring_radius_inner^2;
      within_ring_phase =  @(x,y) (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2<ring_radius_outer^2 & (x-(6+2*i*ring_radius_outer+i*spacing)).^2+(y-ring_middle_y).^2>ring_radius_inner^2;
      
      mod_reg = assign_space(mod_reg, xrange, yrange, within_ring, delta); 
      mod_phi = assign_space(mod_phi, xrange, yrange, within_ring_phase, (i+1)*(pi/2)); 
end


%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 2; 

mode = 0; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, 0.28, eps_clad, eps_wg, mode, hy, pts); 

% Location of the source
src_x = 2.2; 
src_y = 7; 

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = 6.1164; 
%Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0/sqrt(norm_P); 
Jz0((src_ind_x-pts) : (src_ind_x+pts),src_ind_y) = jz0/sqrt(norm_P); 
% Display the structure
figure; visabs(eps_space, xrange, yrange); 
xlabel('x (\mum)'); ylabel('y (\mum)'); 

%% Solve for field distributions

% solver
[Ez, Hx, Hy, omega,Ab] = solve_frequency(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_space, mod_reg, mod_phi, Jz0, Npml);


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
    
    [beta(i), ~, ~, ez(:, i), amp(i)] = beta_sym(L0, wvlen_i, wg_left, eps_clad, eps_wg, mod(i-1-Nsb, 2), hy, pts);
    
    for j = 1:N(1)
        a_sim(i, j) = max(max(abs(Ez_i(j, :)))) / amp(i); 
    end
    
end

