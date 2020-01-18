function [A,b] = condition(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_r, mod_reg, mod_phi, Jz0, Npml)

constants                           % Import constants file

N = size(eps_r);                    % [N(1) N(2)]
L = [diff(xrange) diff(yrange)];    % [Lx Ly]
dL = L./(N);                        % [dx dy]
M = prod(N); 

omega0 = 2*pi*c0/wvlen0;            % Angular frequency
n_sb = -Nsb : 1 : Nsb;              % Sideband array
omega = omega0 + Omega*n_sb;        % Sideband frequency array
wvlen = 2*pi*c0./omega;             % Wavelength frequency array
    
%% Implement PML equations
Sxf = cell(2*Nsb+1, 1); 
Sxb = cell(2*Nsb+1, 1); 
Syf = cell(2*Nsb+1, 1); 
Syb = cell(2*Nsb+1, 1); 
for i = 1 : (2*Nsb + 1)
    [Sxf{i}, Sxb{i}, Syf{i}, Syb{i}] = S_create(L0, wvlen(i), xrange, yrange, N, Npml); 
end

%% Setup space
% Reshape epsilon and phi into 1D array 
vector_eps = reshape(e0 * eps_r, M, 1); 
vec_phi = reshape(mod_phi, M, 1); 
% Create T diags
T_eps = spdiags(vector_eps, 0, M, M); 
T_phi = spdiags(exp(1i*vec_phi), 0, M, M);

%% Build derivative matrices
for n = 0 : 1 : N(2)-1
    block_y(1 + n*N(1) : (n+1)*N(1) , 1 + n*N(1) : (n+1)*N(1)) = -speye(N(1)); 
    if (n < N(2) - 1)
        block_y(1 + n*N(1) : (n+1)*N(1) , 1 + (n+1)*N(1) : (n+2)*N(1)) = speye(N(1)); 
    else
        block_y(1 + n*N(1) : (n+1)*N(1) , 1 : N(1)) = speye(N(1)); 
    end
end
Dyf(1 : (1)*(N(1)*N(2)) , 1:(1)*(N(1)*N(2))) = 1/dL(2) * block_y;
Dyb = -transpose(Dyf); 
block_x = -speye(N(1)) + circshift(speye(N(1)), -1); 
for n = 0 : 1 : (N(2)*1 - 1)
    Dxf(1 + n*N(1) : (n+1)*N(1), 1 + n*N(1) : (n+1)*N(1)) = 1/dL(1) * block_x; 
end
Dxb = -transpose(Dxf); 

%% Reshape Jz0 into a vector
jz0 = reshape(Jz0, M, 1);  

%% Construct A matrix and b vector 
A_i = cell(2*Nsb+1, 1); 
for i = 1 : (2*Nsb + 1)
    A_i{i} = Sxb{i}*Dxb * u0^-1 * Sxf{i}*Dxf + Syb{i}*Dyb * u0^-1 * Syf{i}*Dyf + omega(i)^2*T_eps; 
end
b0 = 1i * omega(Nsb+1) * jz0;
b = zeros(M * (2*Nsb+1), 1); 
b((Nsb*M)+1 : (Nsb+1)*M, 1) = b0; 

%% Sideband coupling
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
    for i = 2 : 2*Nsb
        A((i-1)*M+1 : i*M, (i-1)*M+1 : i*M) = A_i{i}; 
        A((i-1)*M+1 : i*M, (i-2)*M+1 : (i-1)*M) = Cn_nm{i-1}; 
        A((i-1)*M+1 : i*M, i*M+1 : (i+1)*M) = Cn_np{i}; 
    end
else
    A(1:M, 1:M) = A_i{1}; 
end
end
