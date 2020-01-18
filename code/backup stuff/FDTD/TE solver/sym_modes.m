function err=sym_modes(x);
global V

%dispersion equations
y=sqrt(V^2-x^2);

%Mode equations
err=abs(x*tan(x)-y);

