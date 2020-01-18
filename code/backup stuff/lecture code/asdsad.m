%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumantra Chaudhuri, M.S. Electrical Engg.
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
% "Silica ridge waveguide simulation using 2D FDTD"
% 
%              In this program, a 20 micron x 8 micron platform is taken
% with a parallel horizontal ridge waveguide (made of silica i.e n=1.5) of 
% 1 micron width and 20 microns length lying on the centre of the 8 micron 
% width of the platform. The region surrounding the waveguide are filled with 
% air of refractive index, n=1. 
%
%             The smallest dimension in this case is the width of the waveguide 
% at the center i.e 1 micron which is again divided into 'factor' number
% of space steps in the 2D FDTD program below and all dimensions in the
% structure including source wavelength are modeled using this smallest 
% dimension iin terms of numer of space steps.
%
%             An option is given in the program to have Mur’s absorbing boundary 
% condition or perfectly matched layer boundary condition. The 2D FDTD codes
% for the two cases are given in a detailed fashion in the previous programs
% and they are being reused here.
%
%             A color scaled plot of the Ez wave travelling through the
% waveguides starting from the left end of the central waveguide as a line
% source along the width of the central waveguide. The source is hard and
% sinusoidal in nature and its free space wavelength is given as 2 microns.
% The epsilon profile of the platform is shown in the background of the plot 
% to give a feel for the location of the travelling wave. The simulation can 
% be stopped by closing this plot window or by waiting till the plots for
% all the time steps are shown.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Clearing variables in memory and Matlab command screen
clear all;
clc;

%Boundary condition selection, any one can be selected by making it 1
pml=1;
abc=0;

factor=6;
%'factor' :- the number of times the smallest dimension (viz. 0.25 micron in 
%             the given spec) is divided to get one space step (delta)
% Increasing value of factor increases resolution of results but if that is done, 
% to see full result, value of time_tot (total time steps) should be increased. 

%Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 0.25 micron/factor; and factor 
% can be changed)
delta=0.25e-6/factor;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

%Total no of time steps
time_tot=1000;

% Grid Dimension in x (xdim) and y (ydim) directions
ydim=32*factor;%The domain is 30*factor space steps or 32*0.25=8 microns long
xdim=80*factor;%The domain is 80*factor space steps or 80*0.25=20 microns wide

%Free-space wavelength in microns
wav=1;

%Index of the three waveguides
index=3.5;

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);


% Defining of the permittivity profile of the region:-

% Specifying central waveguide location and dimensions i.e full length 
% of the platform viz 20 microns and lying from 14*0.25=3.5 microns to 
% 18*0.25=4.5 microns along the width with a dimension of 1 micron.
%epsilon(:,14*factor+1:18*factor)=index*index*epsilon0;
epsilon(200:250,:)=index*index*epsilon0;
% 2D FDTD update for PML boundary as used in previous program
if pml==1
    
    % Initialization of field matrices
    Ez=zeros(xdim,ydim);
    Ezx=zeros(xdim,ydim);
    Ezy=zeros(xdim,ydim);
    Hy=zeros(xdim,ydim);
    Hx=zeros(xdim,ydim);
    
    % Initializing electric conductivity matrices in x and y directions
    sigmax=zeros(xdim,ydim);
    sigmay=zeros(xdim,ydim);
    
    
    %Perfectly matched layer boundary design
    %Reference:-http://dougneubauer.com/wp-content/uploads/wdata/yee2dpml1/yee2d_c.txt
    %(An adaptation of 2-D FDTD TE code of Dr. Susan Hagness)
    
    %Boundary width of PML in all directions
    bound_width=factor*4;
    
    %Order of polynomial on which sigma is modeled
    gradingorder=6;
    
    %Required reflection co-efficient
    refl_coeff=1e-6;
    
    %Polynomial model for sigma
    sigmamax=(-log10(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width*delta);
    boundfact1=((epsilon(xdim/2,bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact2=((epsilon(xdim/2,ydim-bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact3=((epsilon(bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact4=((epsilon(xdim-bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    x=0:1:bound_width;
    for i=1:1:xdim
        sigmax(i,bound_width+1:-1:1)=boundfact1*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
        sigmax(i,ydim-bound_width:1:ydim)=boundfact2*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
    end
    for i=1:1:ydim
        sigmay(bound_width+1:-1:1,i)=boundfact3*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
        sigmay(xdim-bound_width:1:xdim,i)=boundfact4*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
    end
    
    %Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
    %This is also split into x and y directions in Berenger's model
    sigma_starx=(sigmax.*mu)./epsilon;
    sigma_stary=(sigmay.*mu)./epsilon;
    
    %Multiplication factor matrices for H matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    G=((mu-0.5*deltat*sigma_starx)./(mu+0.5*deltat*sigma_starx));
    H=(deltat/delta)./(mu+0.5*deltat*sigma_starx);
    A=((mu-0.5*deltat*sigma_stary)./(mu+0.5*deltat*sigma_stary));
    B=(deltat/delta)./(mu+0.5*deltat*sigma_stary);
    
    %Multiplication factor matrices for E matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    C=((epsilon-0.5*deltat*sigmax)./(epsilon+0.5*deltat*sigmax));
    D=(deltat/delta)./(epsilon+0.5*deltat*sigmax);
    E=((epsilon-0.5*deltat*sigmay)./(epsilon+0.5*deltat*sigmay));
    F=(deltat/delta)./(epsilon+0.5*deltat*sigmay);
    
    
    NFREQ = 1000;
FREQ  = linspace(0,2*(c/wav),NFREQ);  
K     = exp(-i*2*pi*deltat.*FREQ);

Eref = zeros(1,NFREQ);
% Etrn = zeros(1,NFREQ);
    
    % Update loop begins
    for n=1:1:time_tot
        
        %matrix update instead of for-loop for Hy and Hx fields
        Hy(1:xdim-1,1:ydim-1)=A(1:xdim-1,1:ydim-1).*Hy(1:xdim-1,1:ydim-1)+B(1:xdim-1,1:ydim-1).*(Ezx(2:xdim,1:ydim-1)-Ezx(1:xdim-1,1:ydim-1)+Ezy(2:xdim,1:ydim-1)-Ezy(1:xdim-1,1:ydim-1));
        Hx(1:xdim-1,1:ydim-1)=G(1:xdim-1,1:ydim-1).*Hx(1:xdim-1,1:ydim-1)-H(1:xdim-1,1:ydim-1).*(Ezx(1:xdim-1,2:ydim)-Ezx(1:xdim-1,1:ydim-1)+Ezy(1:xdim-1,2:ydim)-Ezy(1:xdim-1,1:ydim-1));
        
        %matrix update instead of for-loop for Ez field
        Ezx(2:xdim,2:ydim)=C(2:xdim,2:ydim).*Ezx(2:xdim,2:ydim)+D(2:xdim,2:ydim).*(-Hx(2:xdim,2:ydim)+Hx(2:xdim,1:ydim-1));
        Ezy(2:xdim,2:ydim)=E(2:xdim,2:ydim).*Ezy(2:xdim,2:ydim)+F(2:xdim,2:ydim).*(Hy(2:xdim,2:ydim)-Hy(1:xdim-1,2:ydim));
        
        % Source condition incorporating given free space wavelength 'wav'
        % and having a location at the left end of central waveguide just
        % after PML boundary
        tstart=1;
        N_lambda=wav*1e-6/delta;
        Ezx(bound_width,14*factor+1:18*factor)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        Ezy(bound_width,14*factor+1:18*factor)=0.5*sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
                
        Ez=Ezx+Ezy;
        
        %Movie type colour scaled image plot of Ez
        if (mod(n,50)==0)
        h=imagesc(Ez);
        drawnow();
        end
        Eref2(n) = Ez(400,100);
        time(n) = n*deltat;
        for i=1:NFREQ
            Eref(1,i) = Eref(1,i) + (K(i)^(n))*Ez(400,100);
        end
    end
end
% Etrn = Etrn .* (deltat);
% Eref = Eref .* (deltat);

L = length(time);
tr = linspace(min(time),max(time),L);
vr = resample(Eref2,tr);
Ts = tr(2)-tr(1);
Fs = 1/Ts;                                          % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
FTvr = fft(vr)/L;                                   % Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
Iv = 1:length(Fv);

figure(2);
plot(Fv, abs(FTvr(Iv)).*2)
