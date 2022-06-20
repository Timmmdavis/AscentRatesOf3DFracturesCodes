%% Air in Gelatine - Turbulent?
clear
close all

%% Characterizing the physical properties of gelatin, a classic analog for the brittle elastic crust, insight from numerical modeling
rho_r=1001;
g=9.81;
E=1500; 
nu=0.5; 
Kc=50;

%Air:
eta=1.81e-5;
rho_f=1;
% %Water
% eta=8.90e-4;
% rho_f=1000;

%Weight gradient
delta_gamma=(rho_r-rho_f)*g;
%Shear mod
mu=E/(2*(1+nu));
%Critical volume
[V] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
%Fracture size
[v,Dn,c]=AscentVelocityApproximation(V,delta_gamma,mu,nu,eta);
%Reynolds no:
Re=(Dn*v*rho_f)/eta;
if Re<2900 && Re>2300
    disp('Flow is transitioning')
elseif Re>2900
    disp('Flow is turbulent')
else
    disp('Flow is laminar')
end