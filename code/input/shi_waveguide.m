clear all; close all; clc;

%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: nm

xrange = [0 10];  % x boundaries in L0
yrange = [-2.5 2.5];  % y boundaries in L0
Npml = [20 20];  % [Nx_pml, Ny_pml]

wvlen0 = 1.55; % Input wavelength
wvlen1 = 1.50; % Target transition wavelength

hx = 0.05;
hy = 0.05;

N = [diff(xrange)/hx diff(yrange)/hy];

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
delta = 0.8*(eps_wg*eps0); 

mod_x = [1.5, 9.001]; % Length of the modulation

within_mod1 = @(x, y) y > 0 & y < wg_upper & x > mod_x(1) & x < mod_x(2); 
within_mod2 = @(x, y) y > wg_lower & y < 0 & x > mod_x(1) & x < mod_x(2); 

mod_reg = assign_space(mod_reg, xrange, yrange, within_mod1, delta); 
mod_reg = assign_space(mod_reg, xrange, yrange, within_mod2, -delta); 


%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 40; 

mode = 0; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, wg_upper, eps_clad, eps_wg, mode, hy, pts); 

% Location of the source
src_x = 0.9999; 
src_y = -0.0001; 

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = 6.1164; 
Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0/sqrt(norm_P); 

% Display the structure
figure; visabs(eps_space, xrange, yrange); 
xlabel('x (\mum)'); ylabel('y (\mum)'); 

%% Solve for field distributions
Nsb = 1;  % Total number of frequencies is 2*Nsb+1
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

% Solve for the coupling constant for adjacent sidebands
for i = 2:2*Nsb+1
    kappa(i-1) = omega0/8 * sum(hy * ez(:, i)' .* ez(:, i-1)' .* mod_reg(100, src_ind_y-pts : src_ind_y+pts)); 
end

%% Calculate the power amplitude from coupled mode theory
if (Nsb > 0)
    a0 = zeros(2*Nsb+1, 1); 
    a0(Nsb+1) = 1; 
    hz = 0.001; 
    zf = diff(mod_x); 

    % Use cmt_coupling to numerically compute the theoretical modal
    % amplitudes
    [a, z] = cmt_coupling(beta, kappa, kappa, a0, hz, zf); 

    % Plot the power amplitude
    figure; 
    for m = 1:2*Nsb+1
        subplot(2*Nsb+1, 1, m); 
        hold on; 
        plot(x, a_sim(m, :)); 
        plot((z+mod_x(1)), abs(a(m, :))); 
        hold off; 
    end
    
end

legend('MF-FDFD', 'CMT')

%%
figure;
visreal(Ez{1}+Ez{2}+Ez{3},xrange,yrange)

