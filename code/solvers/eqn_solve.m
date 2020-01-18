function [x_val] = eqn_solve(func, x1, x2)
%EQN_SOLVE Summary of this function goes here
%   Detailed explanation goes here

res = 3000; 
x = linspace(x1, x2, res); 

err = abs(func(x)); 

[~, ind] = min(err); 

x_val = x(ind); 

end

