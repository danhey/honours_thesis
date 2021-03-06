function [Sx_f, Sx_b, Sy_f, Sy_b] = S_create(L0, wvlen, xrange, yrange, N, Npml)
%S_CREATE Summary of this function goes here
%   Detailed explanation goes here


%% Set up the domain parameters.
e0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
u0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(e0*u0);  % speed of light in vacuum in L0/sec

% L = [diff(xrange) diff(yrange)];  % [Lx Ly]
% dL = L./N;  % [dx dy]

M = prod(N); 

omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor
% Create the sfactor in each direction and for 'f' and 'b'
s_vector_x_f = create_sfactor(xrange, 'f', omega, e0, u0, N(1), Npml(1)); 
s_vector_x_b = create_sfactor(xrange, 'b', omega, e0, u0, N(1), Npml(1)); 
s_vector_y_f = create_sfactor(yrange, 'f', omega, e0, u0, N(2), Npml(2)); 
s_vector_y_b = create_sfactor(yrange, 'b', omega, e0, u0, N(2), Npml(2)); 


% Fill the 2D space with layers of appropriate s-factors
Sx_f_2D = zeros(N); 
Sx_b_2D = zeros(N); 
Sy_f_2D = zeros(N); 
Sy_b_2D = zeros(N); 

for j = 1:N(2)
    Sx_f_2D(:, j) = s_vector_x_f .^-1;  
    Sx_b_2D(:, j) = s_vector_x_b .^-1; 
end

for i = 1:N(1)
    Sy_f_2D(i, :) = s_vector_y_f .^-1; 
    Sy_b_2D(i, :) = s_vector_y_b .^-1; 
end
% Reshape the 2D s-factors into a 1D s-array
Sx_f_vec = reshape(Sx_f_2D, M, 1); 
Sx_b_vec = reshape(Sx_b_2D, M, 1); 
Sy_f_vec = reshape(Sy_f_2D, M, 1); 
Sy_b_vec = reshape(Sy_b_2D, M, 1); 

% Construct the 1D total s-array into a diagonal matrix
Sx_f = spdiags(Sx_f_vec, 0, M, M); 
Sx_b = spdiags(Sx_b_vec, 0, M, M); 
Sy_f = spdiags(Sy_f_vec, 0, M, M); 
Sy_b = spdiags(Sy_b_vec, 0, M, M); 

end

