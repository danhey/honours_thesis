function [a, x] = cmt_coupling(beta, kp_vec, km_vec, a0, hx, xf)
%% This function implements the equation A.15 to describe power exchange under multiple modes in a modulated waveguide, given inputs of propagation constant, and coupling coefficients between modes
M = length(beta);   % Number of propagation constants
x = 0:hx:xf;        % Define x space with hz step
N = length(x);      % Obtain length of x
a = zeros(M, N);	% CMT results, for M couplings
del_beta = diff(beta); 
kp_x = @(x) kp_vec .* exp(-1i * del_beta * x);	% Coupling up 
km_x = @(x) km_vec .* exp(+1i * del_beta * x);	% Coupling down
K = @(x) diag(kp_x(x), 1) + diag(km_x(x), -1); 	% K matrix (Eqtn A.14)
a(:, 1) = a0; 
% Finite difference scheme
for i = 2:N
    a(:, i) = a(:, i-1) - 1i*K(x(i)) * hx * a(:, i-1); % Eqtn A.15
end
end