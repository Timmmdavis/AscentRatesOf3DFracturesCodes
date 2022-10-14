clear
close all

%With increases...:

%% Maximum
Kc=1e6; %typical rock

SlowerAscent=0
rho_f=2700;
rho_r=2750;

if SlowerAscent==1
    %Height goes up
    VolumeIn=0.014e9;           %m3
    deltagamma=(rho_r-rho_f)*9.81;   
    t=60*60%*24*30*12;     %s
    
    %Less height
    eta=30;   %0.001-0.01 Pa
    E=40e9;    %10-40 GPa
    nu=0.25;    %dmlss
else
    %% Minimum
    %Height goes up
    VolumeIn=0.14e9;           %m3
    deltagamma=(rho_r-rho_f)*9.81;   
    t=60*60%*24*30*12;     %s
    
    %Less height
    eta=10;   %0.001-0.01 Pa
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

%Solve for time to reach height 'x'
HeightIn=2000;
time=(HeightIn^3*eta*mu*pi^2*(-(VolumeIn^3*deltagamma^3*(nu - 1))/mu)^(1/2))/(3*(VolumeIn^3*deltagamma^3 - VolumeIn^3*deltagamma^3*nu));
timeHrs=time/(60*60);


%Integrate rate from 2c up to 24km depth
c=((9*VolumeIn*mu)/(16*deltagamma*(1-nu)))^(1/4);
R=2*c; %Start speed
EndDistance=24*1e3;

FuncToInt = @(z)ComputeRateAtZ(z,eta,mu,VolumeIn,deltagamma,nu);
AvgRate = integral(FuncToInt,R,EndDistance);

v
AvgRate=AvgRate/(EndDistance-R)

function [rate]=ComputeRateAtZ(z,eta,mu,VolumeIn,deltagamma,nu)
    %time at given z
    t=(((z.^6*(eta.^2.*pi.^4.*mu))/(9.*VolumeIn.^(3).*deltagamma.^(3).*(1-nu ))).^(1/6)).^3;
    %rate at given z
    rate=((VolumeIn.^(3).*deltagamma^(3).*(1-nu))./(81.*eta^(2).*pi.^(4).*t.^(4).*mu)).^(1/6);   

end



