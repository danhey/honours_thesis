 clc; 

%%
L0 = 1e-6;
omega = 0.000:0.001:1.5;
wvlen = 1./omega;
d = (1)/2; 
eps1 = 1; 
eps2 = 4; 
m = 0; 
hy = 0.025; 
pts = 100; 

for m = 0:1:4
    for i = 1:length(wvlen)
    [beta, kz, alpha, Ez, A, y] = beta_sym(L0, wvlen(i), d, eps1, eps2, m, hy, pts); 
    %k vals
    TEST1(i,m+1) = beta/(2*pi);
    % omega vals
    TEST2(i,m+1) = 1./wvlen(i);
    end
end

plot(TEST1,TEST2)
xlabel('k')
ylabel('w')