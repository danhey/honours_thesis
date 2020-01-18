function [P] = flux_y(L0, Sx, Sy, ind_x, ind_y, hx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

P = 1*L0*hx * sum(Sy(ind_x(1):ind_x(2), ind_y)); 

end

