clear all; close all; clc; 

%%
L0 = 1e-6; 
wvlen = 1.5; 
d = 0.4; 
eps1 = 1; 
eps2 = 5; 
m = 5; 
hy = 0.025; 
pts = 30; 


% L0 = 1; 
% wvlen = 1.5e-6; 
% d = 0.4e-6; 
% eps1 = 1; 
% eps2 = 5; 
% m = 1; 
% hy = 0.025e-6; 
% pts = 30; 


[beta, kz, alpha, Ez, A, y] = beta_sym(L0, wvlen, d, eps1, eps2, m, hy, pts); 

beta
alpha
plot(y, Ez)