function animate_steadystate(Ez,xrange,yrange)
% This function takes an input steady state field component and produces
% an animation of the field propagation. The solution is exact.

n_phi  = 100;
phase = linspace(0,2*pi,n_phi);
% Draw frames
for nphi = 1 : n_phi
    % Add phase to Ez
    f = Ez*exp(1i*phase(nphi));
    visreal(f,xrange,yrange);
    %mesh(real(f));
    drawnow();
end