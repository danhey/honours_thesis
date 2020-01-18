%% Initialise MATLAB
clc; clear all;


% Import constants
constants

%% Define grid

xrange = [0 5];             % x boundaries in L0
yrange = [0 10];            % y boundaries in L0
Npml = [20 20];             % [Nx_pml, Ny_pml]
NPML = [20 20 20 20];
dx = 0.05;                  % Resolution along x
dy = 0.05;                  % Resolution along y

% Calculate matrix size based on given resolution
N = [ceil(diff(xrange)/dx) ceil(diff(yrange)/dy)];


%% Source position (L0 units)
src_x = ceil(2/dx):ceil(3/dx);
src_y = ceil(8/dy);

%% Define problem

% Source characteristics
wavelength = 1;


frequency = c0/wavelength;
omeg = 2*pi/wavelength;


%% Build device
                      
% Pre-allocate material space
mu_space_xx = ones(N);
mu_space_yy = ones(N);
eps_space = ones(N);

%% Define material permittivity
eps_wg = 16;
eps_clad = 1;

%% Waveguide dimensions

within_wg = @(x, y) y > (4.5) & y < 5.5;
eps_space = assign_space(eps_space, xrange, yrange, within_wg, eps_wg);



%% Compute source
% Courant factor time step
dmin = min([dx dy]);
dt = dmin/(2*c0);


%% Compute PML
% 2x grid formulation
Nx2 = 2*N(1);
Ny2 = 2*N(2);

% COMPUTE sigx
sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end

% COMPUTE sigy
sigy = zeros(Nx2,Ny2);
for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

%% Compute update coefficients
% These are outside the main loop so that they're not calculated at every
% time-step

% Hx coefficients
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./mu_space_xx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./mu_space_xx ./ mHx0;

% Hy coefficients
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./mu_space_yy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./mu_space_yy ./ mHy0;

% Dz coefficients
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

% Ez coefficient
mEz1 = 1./eps_space;

%% Pre-allocate variable memory
% Electromagnetic fields (TE)
Hx = zeros(N); Hy = zeros(N);
Dz = zeros(N); Ez = zeros(N);

% Curls
CEx = zeros(N); CEy = zeros(N); CHz = zeros(N);

% Integrations
ICEx = zeros(N); ICEy = zeros(N); IDz = zeros(N);


%frequency = c0/5e9;
%% Minor pre-loop optimisations
f = 1i*2*pi*frequency;
Nx = N(1); Ny = N(2);


%% Fourier stuff
Fs = 1/dt;
Nyquist = Fs/2; %Nyquist frequency
NFREQ = 5000;
FREQ  = linspace(2.2e+14,3.8e+14,NFREQ);  
K     = exp(-i*2*pi*dt.*FREQ);

REF = zeros(Nx,NFREQ);
TRN = zeros(Nx,NFREQ);
SRC = zeros(1,NFREQ);
% Gaussian pulse
%testForm = zeros(Nx,Ny,NFREQ);
t0 = 4.2e-9;
tau = 2*1.e-9;
N_lambda=c0/(frequency);


% POSITION OF RECORD PLANES
ny_ref = Ny-NPML(3);
ny_trn = 21;

% COMPUTE REFRACTIVE INDICES IN RECORD PLANES
nref = sqrt(eps_space(1,ny_ref)*mu_space_xx(1,ny_ref));
ntrn = sqrt(eps_space(1,ny_trn)*mu_space_xx(1,ny_trn));

%% Calculation loop
for T=1:1500
    source(T) = exp(-(T-1)^2/(1/2*frequency)^2);
    % Curl Ex
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny) - Ez(1:Nx,1:Ny-1))./dy;
    CEx(1:Nx,Ny) = (Ez(1:Nx,1) - Ez(1:Nx,Ny))./dy;
    
    % Curl Ey
    CEy(1:Nx-1,1:Ny) = - (Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))./dx;
    CEy(Nx,1:Ny) = - (Ez(1,1:Ny) - Ez(Nx,1:Ny))./dx;
    
    %TFSF
    CEx(src_x,src_y) = CEx(src_x,src_y) - source(T)/dy;
    
    % Integration terms
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    
    % Update Hx and Hy
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;

    % Curl Hz
    CHz(1,1) = (Hy(1,1) - Hy(N(1),1))/dx - (Hx(1,1) - Hx(1,N(2)))/dy;
    CHz(2:Nx,1) = (Hy(2:Nx,1) - Hy(1:Nx-1,1))./dx  - (Hx(2:Nx,1) - Hx(2:Nx,Ny))./dy;    
    
    CHz(1,2:Ny) = (Hy(1,2:Ny) - Hy(Nx,2:Ny))./dx - (Hx(1,2:Ny) - Hx(1,1:Ny-1))./dy;
    CHz(2:Nx,2:Ny) = (Hy(2:Nx,2:Ny) - Hy(1:Nx-1,2:Ny))./dx - (Hx(2:Nx,2:Ny) - Hx(2:Nx,1:(Ny-1)))./dy;

    %TFSF
    CHz(src_x,src_y) = CHz(src_x,src_y) - source(T)/dy;
    
    % Update Integration Term
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Update Ez
    Ez = mEz1.*Dz;
    
    %Sine source
   % Ez(src_x,src_y)=sin(((2*pi*(c0/(N_lambda))*(T)*dt)));
    %Ez(src_x,src_y)=source(T);
    
    %% Draw field
    % Update field graph every 20 steps. Ironically, graphing the field is
    % the slowest part of the code.
    if (mod(T,50) == 0)
       visreal(Ez,xrange,yrange);
       drawnow();
    end
    
    
    for i=1:NFREQ
       REF(:,i) = REF(:,i) + (K(i)^T)*Ez(:,22)*dt;
       %TRN(:,i) = TRN(:,i) + (K(i).^(T)).*Ez(:,ny_trn);
        SRC(i) = SRC(i) + (K(i).^T).*source(T)*dt;
    end
    
end

for i = 1:NFREQ
    REF2(i) = sum(abs(REF(:,i)).^2);
end

figure(2);
plot(FREQ,REF2);

figure(3);
plot(FREQ,abs(SRC).^2);

%  TRN = TRN .* (dt);
%  REF = REF .* (dt);
%  SRC = SRC .* (dt);
% ref2 = abs(REF./SRC).^2;
% trn2 = abs(TRN./SRC).^2;
% 
% for i = 1:NFREQ
%     REF2(i) = sum(ref2(:,i));
%     TRN2(i) = sum(trn2(:,i));
% end
% 
% figure(2); hold on
% plot(FREQ,REF2)
% plot(FREQ,TRN2)
% plot(FREQ,REF2+TRN2);
% legend('REF','TRN','CON');
% hold off

% for nfreq = 1 : NFREQ
% % Compute Wave Vector Components
% %     lam0 = c0/FREQ(nfreq); %free space wavelength
% %     k0 = 2*pi/lam0; %free space wave number
% %     kyinc = k0*nref; %incident wave vector
% %     m = [-floor(Nx/2):floor(Nx/2)]'; %spatial harmonic orders
% %     kx = - 2*pi*m/Sx; %wave vector expansion
% %     kyR = sqrt((k0*nref)^2 - kx.^2); %ky in reflection region
% %     kyT = sqrt((k0*ntrn)^2 - kx.^2); %ky in transmission region
%     % Compute Reflectance
%     ref = REF(:,nfreq)/SRC(nfreq); %normalize to source
%     ref = fftshift(fft(ref)); %compute spatial harmonics
%     %ref = real(kyR/kyinc) .* abs(ref).^2; %compute diffraction eff.
%     REF2(nfreq) = sum(ref); %compute reflectance
%     % Compute Transmittance
%     trn = TRN(:,nfreq)/SRC(nfreq); %normalize to source
%     trn = fftshift(fft(trn)); %compute spatial harmonics
%    % trn = real(kyT/kyinc) .* abs(trn).^2; %compute diffraction eff.
%     TRN2(nfreq) = sum(trn); %compute transmittance
% end

% TRN = TRN .* (dt);
% REF = REF .* (dt);
% SRC = SRC .* (dt);

% REF2 = abs(REF./SRC).^2;
% TRN2 = abs(TRN./SRC).^2;
% CON = REF2 + TRN2;
% 
% clf;
% 
% % PLOT LINES
% h = plot(FREQ,100*REF,'-r','LineWidth',2);
% hold on;
% plot(FREQ,100*TRN,'-b','LineWidth',2);
% plot(FREQ,100*CON,':k','LineWidth',2);
% hold off;
% 
% % SET VIEW
% xlim([FREQ(1) FREQ(NFREQ)]);
% ylim([0 105]);
% h2 = get(h,'Parent');
% set(h2,'FontSize',14,'LineWidth',2);
% h = legend('Reflectance','Transmittance','Conservation');
% set(h,'Location','NorthOutside');
% 
% % LABEL PLOT
% xlabel('Frequency (GHz)');
% ylabel('%', 'Rotation', 0, 'HorizontalAlignment','right');
% 
% figure(); hold on;
% plot(FREQ,(abs(REF).^2)/length(REF));
% plot(FREQ,(abs(TRN).^2)/length(TRN));
% plot(FREQ,(abs(SRC).^2)/length(SRC));
% legend('Eref','Etrn','src');
% hold off;

% %source transform
% L = length(source);
% tr = linspace(1,2000,L);
% vr = resample(source,tr);
% Ts = tr(2)-tr(1);
% Fs = 1/Ts;                                          % Sampling Frequency
% Fn = Fs/2;                                          % Nyquist Frequency
% FTvr = fft(vr)/L;                                   % Fourier Transform
% Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
% Iv = 1:length(Fv);
% figure(3);
% plot(Fv, abs(FTvr(Iv))*2)