function [beta, kz, alpha, Ey, A, y_vec] = beta_sym(L0, wvlen, d, eps1, eps2, m, hy, pts)
%BETA_SYM Summary of this function goes here
%   Detailed explanation goes here


%% 
eps0 = 8.854e-12 * L0; 
mu0 = 4*pi*1e-7 * L0; 

c0 = 1/sqrt(eps0*mu0); 
omega = 2*pi*c0 / wvlen; 

%%
kz_max = omega .* sqrt(mu0 * eps0 * (eps2 - eps1)); 

alpha_z1 = @(kz) sqrt(kz_max^2 - kz.^2);

alpha_z_even = @(kz) kz .* tan(kz * d); 
alpha_z_odd = @(kz) kz .* tan(kz * d + pi/2); 

trans_eqn_even = @(kz) alpha_z1(kz) - alpha_z_even(kz); 
trans_eqn_odd = @(kz) alpha_z1(kz) - alpha_z_odd(kz); 


if (mod(m, 2) == 0)
    % Even mode
    kz = eqn_solve(trans_eqn_even, (m/2)*pi/d, pi/2/d + (m/2)*pi/d); 
else
    % Odd mode
    kz = eqn_solve(trans_eqn_odd, pi/2/d + (m-1)*pi/2/d, pi/d + (m-1)*pi/2/d);    
end

%%
beta = sqrt(omega^2 * mu0*eps0 * eps2 - kz^2); 
alpha = sqrt(beta^2 - omega^2 * mu0*eps0 * eps1); 

%%
ylim = pts * hy; 
% y_vec = [-ylim : hy : ylim]'; 
y_vec = [-ylim+hy/2 : hy : ylim+hy/2]'; 


Ey_even = @(y) (y > d).*cos(kz*d) .* exp(-alpha * (y-d)) + ...
          (y < -d).*cos(kz*d) .* exp(alpha * (y+d)) + ...
          (abs(y)<d).*cos(kz*y); 
      
Ey_odd = @(y) (y > d).*sin(kz*d) .* exp(-alpha * (y-d)) + ...
          (y < -d).*sin(kz*d) * (-1) .* exp(alpha * (y+d)) + ...
          (abs(y)<=d).*sin(kz*y); 


% Ey_even = @(y) (y > d).*cos(kz*d) .* exp(-alpha * y) + ...
%           (y < -d).*cos(kz*d) .* exp(alpha * y); 
%       
% Ey_odd = @(y) (y > d).*sin(kz*d) .* exp(-alpha * y) + ...
%           (y < -d).*sin(kz*d) * (-1) .* exp(alpha * y); 

if (mod(m, 2) == 0)
    Ey = Ey_even(y_vec); 
else
    Ey = Ey_odd(y_vec); 
end

%% Normalization
mode_integral = sum(hy * Ey.^2); 
A = sqrt(2 * omega * mu0 / (beta * mode_integral)); 

Ey = A*Ey; 

end

