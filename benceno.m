clear all;
clc;
close all;
%-------------------Orden de elementos: C7H8-H2-C6H6-CH4-C12H10-----------------------

global A alfa Diam dQdL Tref Hf Pt a b c d R PM vis alf
ni=[0.04559 0.2577 5.796*10^-3 0.3806 1.2599*10^-4] ; % condiciones iniciales de entrada en Kmol/s
Ten=922; % Temperatura entrada K
Pen=3.447*10^6; % Presión entrada Pa
Pt=500; % Presión total  Psi
PM=[92.138 2.016 78.114 16.043 154.211]; % Peso Molecular  kg/kmol

alfa=[-1 -1 1 1 0;
       0 1 -2 0 1;
       0 -1 2 0 -1]; % matriz de coeficientes de las reacciones
   
%Cpj=[206.8 28.8 168.44 34.3 326.96]; % kj/kmol.K  (Si se considera que la Cp no varía con la temperatura, activar esta línea)
AA=[148.8 28.84 74.06 34.31 326.96];
a= AA*10^(-3);
B=[32.4 0.00765 32.95 5.469 0];
b= B*10^(-5);
C=[0 0.3288 -25.20 0.3661 0];
c=C*10^(-8);
D=[0 -0.8698 77.57 -11.00 0];
d= D*10^(-12);

Hf=[50100 0 82800 -74600 180000]; % kj/kmol , a 25ºC y 1 atm

Diam= 2.045; % Diámetro del reactor en metros; dar valor
Radio= Diam/2;
A= pi*(Radio^2); % m2

Tref= 298.15; % K  es 25ºC
Y0=[ni,Ten,Pen];

R=8.3144; % kj/kmol.K
vis=2.510*10^-5; % viscosidad  Pa.s
alf=0.96;

dQdL= 0;

error=1;
L=10;
while error>10^-3
[x,y]=ode45(@fbenceno,[0 L],Y0);
Conv=(ni(1)-y(:,1))/ni(1); % Conversión de Tolueno
error= 0.5-Conv(end);
L= L+0.05;
end

for i=1: length(x);
    Qv(i)=sum(y(i,(1:5)))*R*y(i,6)/y(i,7);
end

L
Conv(end)

Sel= (y(:,3)-ni(3))./(ni(1)-y(:,1))*(1/-(-1)); % Selectividad a Benceno
Sel(end)

% Representaciones gráficas

figure(1) % Perfil de Temperatura
plot(x,y(:,6),'r')
xlabel('Longitud  (m)')
ylabel('Temperatura  (K)')

figure(2) % nºmoles de cada elemento
plot(x,y(:,[1:5]))
xlabel('Longitud  (m)')
ylabel(' Moles   (Kmol/s)')
legend('C7H8','H2','C6H6','CH4','C12H10')

figure(3) % Conversión de tolueno
plot(x,Conv)
xlabel('Longitud  (m)')
ylabel('X  tolueno')
 
figure(4) % Selectividad a benceno
plot(x,Sel)
xlabel('Longitud  (m)')
ylabel('Selectividad  C6H6')

figure(5) % Variación de Presión
plot(x,y(:,7))
xlabel('Longitud  (m)')
ylabel('Presión    (Pa)')
