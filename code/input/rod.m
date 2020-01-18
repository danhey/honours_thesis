clear all; close all; clc;

Nsb = 1;  % Total number of frequencies is 2*Nsb+1

%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: nm

xrange = [0 11];  % x boundaries in L0
yrange = [0 5];  % y boundaries in L0
% N = [151 51];  % [Nx Ny]
N = [400 400];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

wvlen0 = (1/0.346); % Input wavelength
wvlen1 = (1/0.352); % Target transition wavelength


hx = diff(xrange) / (N(1)); 
hy = diff(yrange) / (N(2)); 

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
eps_clad = 2.25; 

eps_space = eps_clad*ones(N);



%% Waveguide dimensions
wg_upper = 2.34; 
wg_lower = 2.06; 

eps_mainrod = 8.9;
rod_radius = 0.2;
ring_radius_inner = 2.9;
ring_radius_outer = 3.18;
spacing = L0;
position_start = 0.5;
period = L0;

for i = 1:xrange(2)
    for j = 1:xrange(2)
        l = position_start + (i-1);
        m = position_start + (j-1);
        if i==3 && j == 3
            pos_x = l-0.04;
           pos_y = m;
            within_ring = @(x,y) (x-l+0.04).^2 + (y-m).^2 < 0.12^2;
            eps_space = assign_space(eps_space, xrange, yrange, within_ring, 7.4);
        elseif i==9 && j ==3
            within_ring = @(x,y) (x-l).^2 + (y-m).^2 < 0.56^2;
            eps_space = assign_space(eps_space, xrange, yrange, within_ring, 9.1);
        elseif i==7 && j ==3
            square_x = l;
            square_y = m;
            within_ring = @(x,y) x<(l+0.33) & x>(l-0.33) & y>(m-0.31) & y<(m+0.31);
            eps_space = assign_space(eps_space, xrange, yrange, within_ring, 8.9);
        else 
            within_ring = @(x,y) (x-l).^2 + (y-m).^2 < 0.2^2;
            eps_space = assign_space(eps_space, xrange, yrange, within_ring, 8.9);
        end
    end
end

%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

delta_max = 0.002; % Relative permittivity modulation strength
delta = delta_max * eps0 * 8.9; 
l = square_x;
m = square_y;
%I
quad1 = @(x,y) x > (l-0.33) & x < (l) & y > (m) & y <(m+0.31);
mod_reg = assign_space(mod_reg, xrange, yrange, quad1, delta); 
mod_phi = assign_space(mod_phi, xrange, yrange, quad1, 0); 
%II
quad2 = @(x,y) x > (l) & x < (l+0.33) & y > (m) & y <(m+0.31);
mod_reg = assign_space(mod_reg, xrange, yrange, quad2, delta); 
mod_phi = assign_space(mod_phi, xrange, yrange, quad2, pi/2); 
%III
quad3 = @(x,y) x > (l-0.33) & x < (l) & y > (m-0.31) & y <(m);
mod_reg = assign_space(mod_reg, xrange, yrange, quad3, delta); 
mod_phi = assign_space(mod_phi, xrange, yrange, quad3, pi/2); 
%IV
quad4 = @(x,y) x > (l) & x < (l+0.33) & y > (m-0.31) & y <(m);
mod_reg = assign_space(mod_reg, xrange, yrange, quad4, delta); 
mod_phi = assign_space(mod_phi, xrange, yrange, quad4, 0); 



%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 1; 

mode = 1; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, wg_upper, eps_clad, eps_wg, mode, hy, pts); 

% Location of the source
src_x = pos_x-0.12; 
src_y = pos_y;

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = 6.1164; 
Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0/sqrt(norm_P); 

% Display the structure.
figure; visabs(eps_space, xrange, yrange); 
xlabel('x (\mum)'); ylabel('y (\mum)'); 

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

figure();
visreal(Ez{1}+Ez{2}+Ez{3},xrange,yrange);
