function err=asym_modes(hd);
global V

%dispersion equations
qd=sqrt(V^2-hd^2);

%Mode equations
err=abs(-hd*cot(hd)-qd);