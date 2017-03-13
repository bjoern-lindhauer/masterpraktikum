%%%%%%%%%%%%%%%%%%%%%%%% 
clear all 
%Brechungsindex 
n1=1; %Luft 
n2=1-25e-8; %Schicht 
n3=1-75e-8; %Substrat 
%%%%%%%%%%%%%%%%%%%%%% 
%Rauigkeit 
sigma1=55e-11; %Schicht 
sigma2=35e-11; %Substrat 
%%%%%%%%%%%%%%%%%%%%%%%% 
%Schichtdicke 
z2=147e-10; 
%Einfallswinkel 
ai =(0:0.0005:3)*pi /180; 
%Wellenvektoruebertrag 
qz=4*pi /1.54*sin ( ai ) ; 
%Betrag des Wellenvektors 
k=2* pi /1.54 * 1e10 ; 
%z−Komponenten 
kz1=k * sqrt (n1^2- cos ( ai ) .^2) ; 
kz2=k * sqrt (n2^2- cos ( ai ) .^2) ; 
kz3=k * sqrt (n3^2- cos ( ai ) .^2) ; 
%modifizierte Fresnelkoeffizienten 
r12=(kz1-kz2 )./( kz1+kz2 ).*exp( -2* kz1 .*kz2*sigma1^2) ; 
r23=(kz2-kz3 )./( kz2+kz3 ).*exp( -2* kz2 .*kz3*sigma2^2) ; 

x2=exp(-2*i*kz2*z2).*r23 ; 
x1=(r12+x2)./(1+ r12.*x2) ; 

semilogy (qz , 2e7*abs (x1) .^2 , '1' ) ; 
%xlabel ( 'q_z [A^{ − 1}] ' ) ; 
xlabel ( '\ alpha_i ' ) ;
data1 = load('messung_bereinigt.txt' ); 
hold on; 
semilogy ( data1 (: ,1) , data1 (: ,2) , '4' ) ; 
ylabel ( ' intensity ' ); 
legend ( ' Theorie ', 'Messwerte');
 %%%%%%%%%%%%%%%%%%