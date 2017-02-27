%%%%%%%%%%%%%%%%%%%%%%%% 
clear all 
graphics_toolkit gnuplot
%Brechungsindex 
n1=1; %Luft 
n2=1-1e-6; %Schicht 
n3=3-1e-6; %Substrat 
%%%%%%%%%%%%%%%%%%%%%% 
%Rauigkeit 
sigma1=0; %Schicht

sigma2=0; %Substrat 
%%%%%%%%%%%%%%%%%%%%%%%% 
%Schichtdicke 
z2=500e-10; 
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
r12=(kz1-kz2)./(kz1+kz2); 
r23=(kz2-kz3)./(kz2+kz3);
x2=exp(-2*i*kz2*z2).*r23; 
x1=(r12+x2)./(1+r12.*x2);

semilogy(qz,abs(x1).^2);
xlabel('q_z [A^{-1}]'); 
ylabel('intensity'); 
legend('800A-Schicht');
print('schicht.pdf')
set au 
replot 