function [Sx, Sy] = poynting(Ez, Hx, Hy)
%% Input Parameters
% Hz, Ex, Ey: 2D arrays of H- and E-field components


%% Output Parameters
% Sx, Sy: 2D array of x- and y-components of Poynting vector

% Hz_av_x = bwdmean_w(Hz, 'x');
Ez_av_y = bwdmean_w(Ez, 'y');
Sx = -1/2 * real(Ez_av_y .* conj(Hy)); 

% Hz_av_y = bwdmean_w(Hz, 'y'); 
Ez_av_x = bwdmean_w(Ez, 'x'); 
Sy = -1/2 * real(Ez_av_x .* conj(Hx)); 

end