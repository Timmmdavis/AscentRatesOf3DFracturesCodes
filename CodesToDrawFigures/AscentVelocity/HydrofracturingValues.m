clear
close all

%With increases...:

%% Maximum
Kc=1e6; %typical rock

SlowerAscent=1;
rho_f=1000;
rho_r=2700;

if SlowerAscent==1
    %Height goes up
    VolumeIn=25;           %m3
    deltagamma=(rho_r-rho_f)*9.81;   
    t=60*60%*24*30*12;     %s
    
    %Less height
    eta=0.01;   %0.001-0.01 Pa
    E=40e9;    %10-40 GPa
    nu=0.25;    %dmlss
else
    %% Minimum
    %Height goes up
    VolumeIn=25;           %m3
    deltagamma=(rho_r-rho_f)*9.81;   
    t=60*60%*24*30*12;     %s
    
    %Less height
    eta=0.001;   %0.001-0.01 Pa
    E=10e9;    %10-40 GPa
    nu=0.25;    %dmlss
end

%Shear mod
mu=E/(2*(1+nu));

[Height,rate,tr] = RateOfAscentAndCrackLength(VolumeIn,t,deltagamma,nu,mu,eta);

%Critical volume
[CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,deltagamma);
if VolumeIn<CriticalVolume
    error('Below critical volume, increase values')
end
%Fracture size
[v,Dn,c]=AscentVelocityApproximation(VolumeIn,deltagamma,mu,nu,eta);
%Reynolds no:
Re=(Dn*v*rho_f)/eta;
if Re<2900 && Re>2300
    disp('Flow is transitioning')
elseif Re>2900
    disp('Flow is turbulent')
    disp(Re)
else
    disp('Flow is laminar')
end

%Solve for time to reach height 'x'
HeightIn=2000;
time=(HeightIn^3*eta*mu*pi^2*(-(VolumeIn^3*deltagamma^3*(nu - 1))/mu)^(1/2))/(3*(VolumeIn^3*deltagamma^3 - VolumeIn^3*deltagamma^3*nu));
timeHrs=time/(60*60);

