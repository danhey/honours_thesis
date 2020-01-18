function[Ez_src,Hx_src,neff,AZ,ind] = ezmode(URxx,URyy,ERzz,dxp)
% EZMODE     Calculate the Fundamental Mode of a Slab Waveguide
%           for the Ez Mode
%
% [Ez_src,Hx_src,neff,AZ,ind] = fmode(URxx,URyy,ERzz,dxp)
%
% dxp is the normalized grid resolution
% dxp = k0*dx
% DETERMINE NUMBER OF POINTS ON GRID
Nx = length(ERzz);
% CONSTRUCT DIAGONAL MATERIAL MATRICES
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
ERzz = diag(sparse(ERzz(:)));
% BUILD DERIVATIVE OPERATORS
DHX = spdiags(-ones(Nx,1)/dxp,-1,sparse(Nx,Nx));
DHX = spdiags(ones(Nx,1)/dxp,0,DHX);
DEX = spdiags(-ones(Nx,1)/dxp,0,sparse(Nx,Nx));
DEX = spdiags(ones(Nx,1)/dxp,1,DEX);
% SOLVE EIGEN-VALUE PROBLEM
A = full(ERzz + DHX/URyy*DEX);
B = full(inv(URxx));
[AZ,NEFF] = eig(A,B);
NEFF = sqrt(diag(NEFF));
% FIND FUNDMENTAL MODE
[neff,ind] = max(real(NEFF));
Ez_src = AZ(:,ind);
% COMPUTE Hx_src
Hx_src = -neff*(URxx\Ez_src);