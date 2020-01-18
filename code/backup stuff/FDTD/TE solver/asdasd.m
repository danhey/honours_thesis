clc,
clear all;
close all;
%wave length represented by lambda
lambda=1550e-9;
k0=2*pi/lambda;
%typical values of refractive index(n1 and n2)
n2=1.6;
n1=1.5;
d=5e-6;
%V number or frequency number
V=sqrt((n2)^2-(n1)^2)*(pi*d/lambda);

%plot the circle 
f=@(u,v) u^2+v^2-V^2;
b1=ezplot(f,[0 20 0 20]);
set(b1,'color','red')
axis square;

g=@(u,v) u*tan(u)-v;
hold on;
b2=ezplot(g,[0 4*pi 0 20]);
set(b2,'color','blue')

m=@(u,v) -v- u*cot(u);
hold on;
b3=ezplot(m,[0 4*pi 0 20]);
title('Graphic Silution of transverse electric--blue curve is v=utan(u) and green curve is v=-ucot(u)')
xlabel('u');
ylabel('v');
set(b3,'color','green')

x_initial=[1 3];
y_initial=[5 3]; 
figure(2);
%symmetric modes
for i=1:2
intersection_point = fsolve(@(root)[f(root(1),root(2));g(root(1),root(2))],[x_initial(i) y_initial(i)] )
h=2*intersection_point(1)/d
q=2*intersection_point(2)/d
prop_const1=sqrt((k0*n2)^2-h^2)
prop_const2=sqrt((k0*n1)^2+q^2)
%calculate C and D
%assume B=1 and A=0
%Ez is represented by E only
C=cos(h*d/2)/exp(-q*d/2);
D=C;
x1=linspace(-3*d/6,(-d/2),1000);
        H1=D*exp(q*x1);
        p1=trapz(x1,H1.^2);
        intensity1=H1.^2;
        
x2=linspace((-d/2),(d/2),1000);
      H2=cos(h*x2);
       p2=trapz(x2,H2.^2);
       intensity2=H2.^2;
       
x3=linspace(d/2,3*d/6,1000);
     H3=C*exp(-q*x3);
      p3=trapz(x3,H3.^2);
      intensity3=H3.^2;
      
    plot(x1,H1,x2,H2,x3,H3);
    title('Representation of Symmetric modes, TE0 and TE2')
    xlabel('x')
    ylabel('Ey(x)')
    
   
     plot(x1,intensity1,x2,intensity2,x3,intensity3);
    title('Representation of intensities of Symmetric modes, TE0 and TE2')
    xlabel('x')
    ylabel('Ey(x)')
    
    p=p1+p2+p3;
    core_frac=p2/p
hold on;
end



x_initial=[2 5];
y_initial=[4 2];
%asymmetric modes
figure(3);
for i=1:2
intersection_point=fsolve(@(root)[f(root(1),root(2));m(root(1),root(2))],[x_initial(i) y_initial(i)] )
h=2*intersection_point(1)/d
q=2*intersection_point(2)/d
prop_const1=sqrt((k0*n2)^2-h^2)
prop_const2=sqrt((k0*n1)^2+q^2)
%calculate C ,D
%assume A=1,B=0
C=sin(h*d/2)/exp(-q*d/2);
D=-C;
x1=linspace(-3*d/5,(-d/2),1000);
        H1=D*exp(q*x1);
        p1=trapz(x1,H1.^2);
        
x2=linspace((-d/2),(d/2),1000);
      H2=sin(h*x2);
       p2=trapz(x2,H2.^2);
       
x3=linspace(d/2,3*d/5,1000);
     H3=C*exp(-q*x3);
      p3=trapz(x3,H3.^2);
      
    plot(x1,H1,x2,H2,x3,H3);
     title('Representation of Antisymmetric modes, TE1 and TE3')
    xlabel('x')
    ylabel('Ey(x)')
    p=p1+p2+p3;
    core_frac=p2/p
hold on;
end