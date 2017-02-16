%%%%%%%%%%%%%%%%%%%%%%%% 
clear all 
%Brechungsindex 
n1=1; %Luft 
n2=1-1e-6; %Schicht 
n3=1-1e-6; %Substrat 
%%%%%%%%%%%%%%%%%%%%%% 
%Rauigkeit 
sigma1=0; %Schicht

sigma2=6e-10; %Substrat 
%%%%%%%%%%%%%%%%%%%%%%%% 
%Schichtdicke 
z2=0; 
%Einfallswinkel 
ai=(0:0.0005:5)*pi/180; 
%Wellenvektorübertrag 
qz=4*pi/1.54*sin(ai); 
%Betrag des Wellenvektors 
k=2*pi/1.54*1e10; 
%z-Komponenten 
kz1=k*sqrt(n1^2-cos(ai).^2); 
kz2=k*sqrt(n2^2-cos(ai).^2); 
kz3=k*sqrt(n3^2-cos(ai).^2); 
%modifizierte Fresnelkoeffizienten 
r12=(kz1-kz2)./(kz1+kz2).*exp(-2*kz1.*kz2*sigma1^2); 
r23=(kz2-kz3)./(kz2+kz3).*exp(-2*kz2.*kz3*sigma2^2);
x2=r23; 
x1=(r12+x2)./(1+r12.*x2); 

%sigma2=6e-10;
r23=(kz2-kz3)./(kz2+kz3).*exp(-2*kz2.*kz3*sigma2^2);
x2=r23; 
x1_new=(r12+x2)./(1+r12.*x2); 


semilogy(qz,abs(x1).^2,qz,abs(x1_new).^2);
xlabel('q_z [A^{-1}]'); 
ylabel('intensity'); 
legend('sigma=0A','sigma=6A');
set au replot 
%%%%%%%%%%%%%%%%%%