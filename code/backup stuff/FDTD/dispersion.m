clear all;
%Dispersion calc
L0 = 1e-6;
wvlen = linspace(0.1,1,1000);
for mode = 0:2
    i = 1;
    for j = wvlen
        [beta, kz0, ~, ~, ~, ~] = beta_sym(L0, j, (1.1), 1, 12.25, mode, 0.025, 30); 
        kz(i,mode+1) = beta;
        i = i+1;
    end
end
omega = 1./wvlen;
plot(omega,kz);
