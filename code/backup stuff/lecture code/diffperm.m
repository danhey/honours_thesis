%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "2D FDTD solution for open boundary in a spatial domain with
%                two different permittivities "
% 
% Objective of the program is to solve for the Maxwell's equation for a TM 
% wave containing the xy-plane polarized magnetic field having components Hy
% and Hx and z-polarized electric field Ez. The fields are updated at every 
% timestep, in a unitless space, where all physical parameters of free space 
% are normalized to 1. The update is done using standard update equations 
% obtained from the difference form of Maxwell's curl equations. The field 
% points are defined in a grid described by Yee's algorithm. The H fields are 
% defined at every half coordinate of spacesteps. More precisely, the Hx part
% is defined at every half y-coordinate and full x-coordinate and the Hy part
% is defined at every half x-coordinate and full y-coordinate and E fields 
% i.e the Ez part is defined at every full x and full y-coordinate points. 
% Also, the time update is done using Leapfrog time-stepping. Here, H-fields
% i.e Hx and Hy are updated every half time-step and E fields i.e Ez are 
% updated every full time-step. This is shown by three alternating for-loop 
% updates in groups of two (for H) and one (for E) spanning entire spatial 
% grid inside a main for-loop for time update spanning the entire time-grid.
%
% A source of electric field is defined at the center (vertical line) of the 
% spatial domain, which is hard , in that it does not change its value due to 
% interference from external fields i.e in other words, the source is a 
% perfect electric conductor. The form of the source can be varied using 
% the variables sine, gaussian and impulse. The source is available in four 
% standard forms- Gaussian and Sinusoidal forms. The permittivity on the two 
% sides of the vertical line source is different to show difference in
% traveling velocity and wavelength in the two materials.
%
% The color scaled plot of Ez field over the entire spatial domain is shown 
% at every time step. The simulation can be ended by closing this plot window 
% or by waiting till all the time step updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clearing variables in memory and Matlab command screen
clear all;
clc;

% Grid Dimension in x (xdim) and y (ydim) directions
xdim=500;
ydim=500;

% Total no of time steps
time_tot=200;

% Position of the source (center of the field domain)
xsource=250;

% Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
epsilon(1:xdim/2,1:ydim)=2.25*epsilon0;
% Left part of domain containing glass of n=1.5
mu=mu0*ones(xdim,ydim);

% Initializing electric and magnetic conductivity matrices
sigma=4e-4*ones(xdim,ydim);
sigma_star=4e-4*ones(xdim,ydim);

%Choice of nature of source
gaussian=1;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
%Choose any one as 1 and other as 0. 


% Update loop begins
for n=1:1:time_tot
    
    % Setting spatial update limits
    n1=1;
    n2=xdim-1;
    n11=1;
    n21=ydim-1;
        
    %Vector update instead of for-loop for Hy and Hx fields
    Hy(n1:n2,n11:n21)=((mu(n1:n2,n11:n21)-0.5*deltat*sigma_star(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21))).*Hy(n1:n2,n11:n21)+(deltat/delta)*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=((mu(n1:n2,n11:n21)-0.5*deltat*sigma_star(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21))).*Hx(n1:n2,n11:n21)-(deltat/delta)*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21));
    
    
    %Vector update instead of for-loop for Ez field
    Ez(n1+1:n2+1,n11+1:n21+1)=((epsilon(n1+1:n2+1,n11+1:n21+1)-0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))).*Ez(n1+1:n2+1,n11+1:n21+1)+(deltat/delta)*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1)-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1));
    
    % Smoothening ends of the line source
    Ez(1:xdim,1)=Ez(1:xdim,2);
    Ez(1:xdim,ydim)=Ez(1:xdim,ydim-1);
    
    % Source conditions
    %if sine
    if sine==1
        tstart=1;
        N_lambda=c/(frequency*delta);
        Ez(xsource,1:ydim)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
    end
    %if gaussian
    if gaussian==1
        if n<=42
            Ez(xsource,1:ydim)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
        else
            Ez(xsource,1:ydim)=0;
        end
    end
        
    %Movie type colour scaled image plot of Ez
    imagesc(delta*(1e+6:1e+6:xdim*1e+6),(delta*(1e+6:1e+6:ydim*1e+6))',Ez',[-1,1]);colorbar;
    title(['\fontsize{20}Color-scaled image plot of Ez with \epsilon_l_e_f_t = 2.25 and \epsilon_r_i_g_h_t = 1 at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    