%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculates the modes of a symmetric slab waveguide
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% variables:
%   lr: wavelength
%   n1: cladding index
%   n2: core index
%   lr: wavlength (in micrometers)
%   d: vector with film thicknesses
%   pol: polarization state (1: TE, 2:TM)
%
% All length units are in micrometers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
global V

%waveguide parameters
%required parameters: n1, n2, d, lr
n1=1;   %top region (the cover)
n2=sqrt(12.25);  %waveguide film indices
d=0.22;      %waveguide film thicknesses
lr=(1/0.8879);    %wavelength
ko=2*pi/lr;
pol = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equations:
%
%Dispersion equations
%   h=sqrt((n2*ko)^2-beta^2)
%   q=sqrt(beta^2-(n1*ko)^2)
%combine these 2 equations
%    h^2+q^2 = (n2^2-n1^2)*ko^2
%
%Mode equations
%   h*tan(h*d/2) = q (symmetric TE modes)
%   -h*cot(h*d/2)= q (antisymmetric TE modes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multiply both sides of the equations by d
% call x=hd, y=qd
%
% h*tan(h*d/2) = q  becomes x*tan(x/2)=y2
% -h*cot(h*d/2) = q  becomes -x*cot(x/2)=y3
%
%   h^2+q^2 = (n2^2-n1^2)*ko^2 becomes x^2+y1^2 = d^2*ko^2*(n2^2-n1^2)  = V
%

%normalized frequency
V=d*pi/lr*sqrt(n2^2-n1^2)


%break x up into 2 pi intervals
intervals=ceil(1.01*V)

%(1) plot the three equations
%   Dispersion equation
%   Symmetric modes
%   Antisymmetric modes

figure(1)
clf
hold 'on'
if intervals==1
    x=linspace(0, 1.5*V,101);
    x=sort([x,V]);
    I=find(x<=V);
    y1=sqrt(V^2-x.^2);  %dispersion
    y2=x.*tan(x/2);     %symmetric
    y3=-x.*cot(x/2);    %antisymetric
    plot(x(I),y1(I),'b')
    plot(x,y2,'r')
    plot(x,y3,'g')
    grid
    axis image
    axis([0, max(x), 0 max(x)])
    xlabel('h d')
    ylabel('q d')
    hold 'off'

else
    
    for lp=1:intervals
        %only plot in pi/2 intervals to eliminate positive/negative changes
        x=linspace(pi/2*(lp-1),pi/2*lp,101);
        x=x(2:end);
        I=find(x<V);
        y1=sqrt(V^2-x.^2);  %dispersion
        y2=x.*tan(x);     %symmetric
        y3=-x.*cot(x);    %antisymetric
        plot(x(I),y1(I),'b')
        plot(x,y2,'r')
        plot(x,y3,'g')
    end
    grid
    axis image
    axis([0, pi/2*intervals, 0, pi/2*intervals])
    xlabel('hd/2')
    ylabel('qd/2')
    hold 'off'
    
end


%(2) Calculate the symmetric modes
%number of symmetric modes
%tan(hd/2)=0
%   zero crossings every pi
%   hd every 2*pi up to V
n_sym= floor(V/(pi));
x_sym=[];
for lp_sym=0:n_sym
    tmp=fminbnd('sym_modes',lp_sym*pi,lp_sym*pi+pi/2)
    x_sym=[x_sym,tmp]
end

%(3) Calculate the anti-symmetric modes
%number of antisymmetric modes
%cot(hd/2)=0
%   zero crossings every pi/2
%   x every odd multiple of pi
n_asym= round(V/(pi));

x_asym=[];
if n_asym>0
    for lp_asym=1:n_asym
        tmp=fminbnd('asym_modes',pi/2+(lp_asym-1)*pi,pi+(lp_asym-1)*pi)
        x_asym=[x_asym,tmp]
    end
end

hd=2*sort([x_sym,x_asym]);
beta=sqrt((n2*ko)^2-(hd/d).^2)
neff=beta/ko

%(4) plot the mode fields
figure(2)

x1=linspace(-d/2,d/2,101);
x2=linspace(d/2,2*d,101);
x3=linspace(-2*d,-d/2,101);

nn=0;
while nn<length(beta)
    nn=nn+1;
    subplot(length(beta),1,nn)
    
    if rem(nn,2)==0 %even mode
        'even'
        subplot(length(beta),1,nn)
        h=sqrt((n2*ko)^2-beta(nn).^2);
        q=sqrt(beta(nn).^2-(n1*ko)^2);
        C=sin(h*d/2)*exp(q*d/2);
        E1=sin(h*x1);
        E2=C*exp(-q*x2);
        E3=-C*exp(q*x3);
        
        plot(x1,E1)
        hold 'on'
        plot(x2,E2)
        plot(x3,E3)
        plot([-d/2 -d/2], [-1 1],'k')
        plot([d/2 d/2], [-1 1],'k')
        hold 'off'
        
        
 
    else    %odd mode       
        'odd'
        h=sqrt((n2*ko)^2-beta(nn).^2);
        q=sqrt(beta(nn).^2-(n1*ko)^2);
        C=cos(h*d/2)*exp(q*d/2);
        E1=cos(h*x1);
        E2=C*exp(-q*x2);
        E3=C*exp(q*x3);
        plot(x1,E1)
        hold 'on'
        plot(x2,E2)
        plot(x3,E3)
        plot([-d/2 -d/2], [-1 1],'k')
        plot([d/2 d/2], [-1 1],'k')
        hold 'off'
 
    end
end

    
neff    



        