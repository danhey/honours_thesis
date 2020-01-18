function [P] = flux_x(L0, Sx, Sy, ind_x, ind_y, hy)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

P = hy * sum(Sx(ind_x, ind_y(1) : ind_y(2))); 

end

