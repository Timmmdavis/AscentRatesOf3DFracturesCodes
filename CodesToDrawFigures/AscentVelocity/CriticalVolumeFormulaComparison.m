
clear
close all

%Test critical volume formulations
nu=0.45;
mu=13e9;
Kc=12e6;
deltagamma=1500*9.81;

E=(2*mu)*(1+nu);
Eprime=E/(1-nu^2);

%% In order of derivation

%Dahm 2000 - Geophys. J. Int. Shape and velocity
a=(Kc/(sqrt(pi)*deltagamma))^(2/3);
Vd=(2/3)*((1-nu)/mu)*(Kc/sqrt(pi))*a^(5/2); % Eq.A2.5 - note he has a range - see A3 (Vd*0.65:*0.133)

%Salimzadeh et al 2020 GRL - critical vol
Fg=1.57;
Vs=((16/27)*((pi^4*Kc^8)/(Eprime^3*Fg^8*deltagamma^5)))^(1/3);%Eq.7

%Davis et al 2020 GRL - critical vol
V=((1-nu)/(16*mu))*((9*pi^4*Kc^8)/(deltagamma^5))^(1/3);

%Smitterello et al Tectonophysics 2021 - critical vol
a=(Kc/(deltagamma*sqrt(pi)))^(2/3); %2D crit len - See Pollard and Townsend and refs therin...
L=a*2; 
alpha=0.3; %0.22->0.33
Vds=alpha*(1-nu^2)*(deltagamma/E)*L^4; %Eq.11

%Mori and Lecampion 2021 - critical when over one
Bk=(deltagamma*Eprime^(3/5)*V^(3/5))/(Kc^(8/5));

%Garagash and Germanovich Aug 2022 Arxiv - Page 16
Kprime=sqrt(2/pi)*Kc;
Eprime=(1/pi)*(E/(1-nu^2));
Vg=0.408*(Kprime^(8/3)/(Eprime*(deltagamma)^(5/3)));


RatioDavisToDahm=V/Vd
RatioDavisToSalim=V/Vs
RatioDavisToSmitt=V/Vds
RatioDavisToMori=Bk
RatioDavisToGerman=V/Vg