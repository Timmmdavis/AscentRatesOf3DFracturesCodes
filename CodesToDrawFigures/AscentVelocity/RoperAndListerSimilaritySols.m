%Roper solution
clear
close all

%% Draw the constant area similarity solutions of Roper and Lister 2007. 

%X=x/K^(2/3); %Eq.3.1

A=(pi/2)*4;
A0=pi/2;
At=A-A0; %Area of tail

FractureHeight=40; %In Z coord
Z=FractureHeight-2; %Eq.6.6 - see fig.3 also
%Rearrange 6.6
T=((Z^3)*16)/(27*At^2);
%Eqs.6.6 to 6.8 Roper and Lister
H=sqrt(Z/(3*T));
K=((4*T)/At)^(1/4);
%Draw profile
Ztail=linspace(0,Z,1000);
Htail=sqrt(Ztail/(3*T));
Zhead=linspace(0,2,1000);
Hhead=(1/2).*sqrt(Zhead).*(2-Zhead).^(3/2);

figure;hold on
plot(Htail,Ztail)
plot(fliplr(Hhead),Zhead+Z);
scatter(H,Z,'k');
ylabel('Z')
xlabel('H')
title('Roper and Lister Fig.8')
WhiteFigure

% %Dimensions:
%Constants
G=8e9;
nu=0.25;
VolumeIn=1.95;
delta_gamma=1000*9.81;
Kc=2e6;
eta=0.05;
t=100000%e9%154000;

%Checking Roper and Lister is correct on dimensional area of Weertman head
%(not quite - needs dividing by (pi)
c=(Kc/(delta_gamma*sqrt(pi)))^(2/3);
P0=0.5*delta_gamma*c;%Twns and Poll
AA=(pi*(1-nu)*P0*c^2)/(G);%Davis Healy Rivalta - Eq.3 
A0=((pi*Kc^2*(1-nu))/(2*G*delta_gamma))/(pi); %Area of head (weertman crack area) - Roper below Eq.6.5

%Convert to 2D area with approx...
[v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,G,nu,eta);
A=2*c*Dn; %Area we start with
A0=((Kc^2*(1-nu))/(2*G*delta_gamma)); %Area of head (weertman crack area)
at=A-A0; %Area of tail


[c,zpnts,hpnts,z,h]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,G,nu,eta,Kc,t);

figure;hold on
plot(hpnts,zpnts)
scatter(h,z,'k');
ylabel('height from crack base')
xlabel('width - opening D_n')
title('Roper and Lister Fig.8 - dimensional')
WhiteFigure

