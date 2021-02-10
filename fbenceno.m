
function dydL=fbenceno (L,y) % L es la variable independiente e y es la variable dependiente (5 de los moles, T y P)
nj=y(1:5); % 5 sustancias químicas, flujo molar en kmol/s
T=y(6); % K
P=y(7); % Pa
global A alfa Diam dQdL Tref Hf Pt a b c d R PM vis alf
% Orden de elementos: C7H8-H2-C6H6-CH4-C12H10 

pj=(nj/sum(nj))*Pt; % psi

Qv=sum(nj)*R*T/P; % Caudal volumétrico  m3/s
m=nj'.*PM; % Caudal molar  kg/s 
Re=sum(m)*Diam/A/vis; % Reynolds
F=0.046*Re^(-0.2); % Factor de fricción

Href= Hf*alfa'; % kj/kmol (1x3)
Cpj= ((a+b*(T-273)+ c*(T-273)^2+d*(T-273)^3))*1000; % kj/kmol.K
f=(a*((T-273)-(Tref-273))+b/2*((T-273)^2-(Tref-273)^2)+c/3*((T-273)^3-(Tref-273)^3)+(d/4)*((T-273)^4-(Tref-273)^4))*1000;
H=Href + f*alfa'; % kj/kmol   (1x3 + 1x5 * 5x3 = 1x3)

%   Para las dos reacciones implicadas en el proceso de hidrodesalquilación, su velocidad se calcula mediante las siguientes expresiones presentes en bibliografía:
%	  Reacción principal (1):

%   r_1= 3,6858∙(10)^6∙exp(-(2,5616∙(10)^4)/T)∙p_T∙(p_H)^0,5   

%	  Reacción secundaria, al ser reversible se tienen dos expresiones de velocidad (2,3):
%   r_2= 0,62717∙exp(-(1,5362∙(10)^4)/T)∙(p_B)^2     


%   r_3= 0,08124∙exp(-(1,2237∙(10)^4)/T)∙p_D∙p_H    

%   Donde:
%    r_i : Velocidad de cada una de las reacciones implicadas en lbmol/(min·ft3).
%    T : Temperatura expresada en K.
%    p_j: Presión parcial del componente en psi.

% Ecuaciones cinéticas:
r(1)=3.6858*10^6*exp(-(2.5616*10^4)/T)*pj(1)*pj(2)^0.5*(0.453592/(60*0.3048^3)); % velocidad de reacción 1 en kmol/m3.s (siendo(0.453592/(60*0.3048^3)el factor de conversión)
r(2)=0.62717*exp(-(1.5362*10^4)/T)*pj(3)^2*(0.453592/(60*0.3048^3));             % velocidad de reacción 2 en kmol/m3.s
r(3)=0.08124*exp(-(1.2237*10^4)/T)*pj(5)*pj(2)*(0.453592/(60*0.3048^3));         % velocidad de reacción 3 en kmol/m3.s

dnjdL=A*r*alfa; % (5) BALANCES MOLARES, siendo r de cada una de las reacciones y alfa los coeficientes en cada una de las reacciones
dTdL=(dQdL-A*(r*H'))/(Cpj*nj); % BALANCE DE ENERGÍA
dPdL=(2*F*Qv/Diam + R/alf/P*(dTdL*sum(nj)+T*sum(dnjdL)))/(R*T/(alf*P^2)*sum(nj)-(A^2/(PM*nj))); %BALANCE DE CANTIDAD DE MOVIMIENTO
dydL=[dnjdL,dTdL,dPdL]';
