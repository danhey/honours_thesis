function [Ez, Hx, Hy, omega,A] = solveTM_coupled_sb(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_r, mod_reg, mod_phi, Jz0, Npml)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Ez, Hx, Hy: Nx-by-Ny arrays of H- and E-field components with 2*Nsb+1
% cells

% dL: [dx dy] in L0
% A: system matrix of A x = b
% omega: angular frequency for given wvlen

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./(N);  % [dx dy]

M = prod(N); 

omega0 = 2*pi*c0/wvlen0;  % angular frequency in rad/sec
n_sb = -Nsb : 1 : Nsb; 

omega = omega0 + Omega*n_sb; 
wvlen = 2*pi*c0./omega; 

%% Deal with the s_factor
Sxf = cell(2*Nsb+1, 1); 
Sxb = cell(2*Nsb+1, 1); 
Syf = cell(2*Nsb+1, 1); 
Syb = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb + 1)
    [Sxf{i}, Sxb{i}, Syf{i}, Syb{i}] = S_create(L0, wvlen(i), xrange, yrange, N, Npml); 
end

%% Set up the permittivity and permeability in the domain.


eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
eps_z = bwdmean_w(eps0 * eps_r, 'x');
eps_z = bwdmean_w(eps_z, 'y'); 
mu_z = mu0 .* ones(N);

% Reshape epsilon into 1D array 
vector_eps_x = reshape(eps_x, M, 1); 
vector_eps_y = reshape(eps_y, M, 1); 
vector_eps_z = reshape(eps_z, M, 1); 
vector_eps = reshape(eps0 * eps_r, M, 1); % <-- I ended up defining Ez to coincide with epsilon

% Reshape phi
vec_phi = reshape(mod_phi, M, 1); 

vector_mu_z = reshape(mu_z, M, 1); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_x = spdiags(vector_eps_x, 0, M, M); 
T_eps_y = spdiags(vector_eps_y, 0, M, M); 
T_eps_z = spdiags(vector_eps_z, 0, M, M); 
T_eps = spdiags(vector_eps, 0, M, M); 

T_mu_z = spdiags(vector_mu_z, 0, M, M); 

T_phi = spdiags(exp(1i*vec_phi), 0, M, M); 

%% Construct derivate matrices
tic
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N); 
time_Dws = toc

%% Reshape Jz0 into a vector
jz0 = reshape(Jz0, M, 1); 

%% Construct A matrix and b vector (Not the most efficient way of doing it)

tic
A_i = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb + 1)
%     A_i{i} = Sxb{i}*Dxb * mu0^-1 * Sxf{i}*Dxf + Syb{i}*Dyb * mu0^-1 * Syf{i}*Dyf + omega(i)^2*T_eps_z; 
    A_i{i} = Sxb{i}*Dxb * mu0^-1 * Sxf{i}*Dxf + Syb{i}*Dyb * mu0^-1 * Syf{i}*Dyf + omega(i)^2*T_eps; 

end
time_A_operators = toc


b0 = 1i * omega(Nsb+1) * jz0;

b = zeros(M * (2*Nsb+1), 1); 

b((Nsb*M)+1 : (Nsb+1)*M, 1) = b0; 

%% Account for coupling (Deinitely not the most efficient way)
tic 
if (Nsb > 0)
    delta_vec = reshape(mod_reg, M, 1); 
    T_delta = spdiags(delta_vec, 0, M, M); 

    Cn_np = cell(2*Nsb, 1); 
    Cn_nm = cell(2*Nsb, 1); 

    for i = 1 : (2*Nsb)
        Cn_np{i} = 1/2 * omega(i)^2 * T_delta * conj(T_phi); 
        Cn_nm{i} = 1/2 * omega(i+1)^2 * T_delta * T_phi; 
    end

    A = sparse((2*Nsb+1)*M, (2*Nsb+1)*M); 

    A(1:M, 1:M) = A_i{1}; 
    A(1:M, M+1:2*M) = Cn_np{1}; 


    A(2*Nsb*M+1 : (2*Nsb+1)*M, 2*Nsb*M+1 : (2*Nsb+1)*M) = A_i{(2*Nsb+1)}; 
    A(2*Nsb*M+1 : (2*Nsb+1)*M, (2*Nsb-1)*M+1 : 2*Nsb*M) = Cn_nm{2*Nsb}; 

    % This is very very inefficient, haha
    for i = 2 : 2*Nsb
        A((i-1)*M+1 : i*M, (i-1)*M+1 : i*M) = A_i{i}; 
        A((i-1)*M+1 : i*M, (i-2)*M+1 : (i-1)*M) = Cn_nm{i-1}; 
        A((i-1)*M+1 : i*M, i*M+1 : (i+1)*M) = Cn_np{i}; 
    end

else
    A(1:M, 1:M) = A_i{1}; 
end

time_setup = toc


%% Solve the equation.
tic 
if all(b==0)
	ez = zeros(size(b));
else
	ez = A\b;
end

time_solve = toc



%% Extract field values
ez_i = cell(2*Nsb+1, 1); 
hx_i = cell(2*Nsb+1, 1); 
hy_i = cell(2*Nsb+1, 1); 

Ez = cell(2*Nsb+1, 1); 
Hx = cell(2*Nsb+1, 1); 
Hy = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb+1)
    ez_i{i} = ez((i-1)*M + 1 : i*M, 1); 
    
    hx_i{i} = -1/(1i*omega(i)) * mu0^-1 * Syf{i}*Dyf * ez_i{i}; 
    hy_i{i} = 1/(1i*omega(i)) * mu0^-1 * Sxf{i}*Dxf * ez_i{i}; 
    
    
    Ez{i} = reshape(ez_i{i}, N); 
    Hx{i} = reshape(hx_i{i}, N); 
    Hy{i} = reshape(hy_i{i}, N); 
end

end
