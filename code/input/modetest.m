clear all; close all; clc;

Nsb = 0;  % Total number of frequencies is 2*Nsb+1
src_x = 1; 
%% Set up the domain parameters.
% L0 = 1e-6;  % length unit: nm
L0 = 1e-6;  % length unit: um

xrange = [0 60];  % x boundaries in L0
yrange = [0 5];  % y boundaries in L0
% N = [151 51];  % [Nx Ny]
Npml = [20 20];  % [Nx_pml, Ny_pml]

hx = 0.03;
hy = 0.03;

N = [ceil(diff(xrange)/hx) ceil(diff(yrange)/hy)];

src_max1 = 1i;
src_max2 = 1; 

%%
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

wvlen1 = 1/0.129; % Input wavelength
wvlen0 = 1/0.197; % Target transition wavelength

omega0 = 2*pi*c0 / wvlen0; 

Omega = 2*pi*c0* (1/wvlen1 - 1/wvlen0); 

wvlen_m1 = 2*pi*c0/(omega0 - Omega); 

%% Set up the permittivity.
eps_wg = 12.25; 
eps_clad = 1;
eps_space = ones(N);



%% Waveguide dimensions
wg_upper = 3.05;
wg_lower = 1.95;

within_wg = @(x, y) y > wg_lower & y < wg_upper;

eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);

%% Modulation region
mod_reg = zeros(N); % This matrix contains the modulation strength
mod_phi = zeros(N); % This matrix contains the modulation phase

delta_max = 1; % Relative permittivity modulation strength
delta = .1*(eps_wg*eps0); 

modLeft = @(x, y) x> 5.2 & x<24.2 & y > 2.5 & y < wg_upper;
modRight = @(x, y) x> 35.8 & x<54.8 & y > 2.5 & y < wg_upper;

% 

%% Set up the source
%%%%%%%%%%%%%%%%%% Single mode source %%%%%%%%%%%%%%%%
Jz0 = zeros(N);

pts = 50; 

mode = 1; 

% Use the function beta_sym() to find the modal profile and beta
[beta0, kz0, alpha0, jz0, A0, y] = beta_sym(L0, wvlen0, (wg_upper-wg_lower), 1, eps_wg, 1, hy, pts); 

% Location of the source

src_y = 2.5000; 

src_ind_x = round((src_x-xrange(1)) / diff(xrange) * N(1)) + 1; 
src_ind_y = round((src_y-yrange(1)) / diff(yrange) * N(2)) + 1; 

norm_P = 1; 
Jz0(src_ind_x, (src_ind_y-pts) : (src_ind_y+pts)) = jz0/sqrt(norm_P); 

% Display the structure
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
    
    [beta(i), ~, ~, ez(:, i), amp(i)] = beta_sym(L0, wvlen_i, (wg_upper-wg_lower), eps_clad, eps_wg, mod(i-1-Nsb, 2), hy, pts);
    
    for j = 1:N(1)
        a_sim(i, j) = max(max(abs(Ez_i(j, :)))) / amp(i); 
    end
    
end

if (Nsb > 0)
    a0 = zeros(2*Nsb+1, 1); 
    a0(Nsb+1) = 1; 
    hz = 0.001; 

    % Use cmt_coupling to numerically compute the theoretical modal
    % amplitudes

    % Plot the power amplitude
    figure; 
    for m = 1:2*Nsb+1
        subplot(2*Nsb+1, 1, m); 
        hold on; 
        plot(x, a_sim(m, :)); 
        hold off; 
    end
    
end

figure();
visreal(Ez{1}+Ez{2}+Ez{3}, xrange, yrange);

% figure(4); hold on;
% title('Hx');
% for i=1:2*Nsb+1
%     
%     subplot(2*Nsb+1,1,i)
%     visreal(Hx{i},xrange,yrange)
% end
% hold off;
% 
% figure(5);
% for i=1:2*Nsb+1
%     
%     subplot(2*Nsb+1,1,i)
%     visreal(Hy{i},xrange,yrange)
% end
