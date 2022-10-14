function [vradial,vtilde,tmatch,toughviscflag,tmk] = RadialVsAscentFrontVelocity(Kc,nu,E,eta,t,Q0,delta_gamma)
% RadialVsAscentFrontVelocity defines the time until the ascent speed of a fracture
% begins to dominate. Before this the radial growth due to injection is far
% faster. The aim is to understand when to use the formulation of constant
% injection i.e. Mori and Lecampion (2022) arXiv
%% Inputs
%Kc - fracture toughness
%nu - pr
%E - Youngs mod
%eta - visc
%t - current time
%Q0 - injection rate
%delta_gamma - stress gradient
%% Outputs
%L - fracture length (radial)
%v - fracture front velocity (radial)
%tmatch - the time at which the front velocities should match 
%toughviscflag - zero if viscous radial fracture, one if toughness
mu=E/(2*(1+nu));

VolumeIn=Q0*t;
%% Now  
%Max ascent speed - circular...
[vtilde,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,eta);

%% Detournay 2016: Mechanics of Hydraulic Fractures. Annu. Rev. Fluid Mech. 2016. 48:311–39
%Eq.1
Kprime=4*sqrt(2/pi)*Kc;
Eprime=E/(1-nu^2);
muprime=eta*12;

tmk=sqrt((Eprime^13*muprime^5*Q0^3)/(Kprime^(18)));


if t<tmk
    %% Viscous dominated Detournay 2016
    mscl=0.6955; %Eq.35
    L=mscl*((Eprime*Q0^3*t^4)/(muprime))^(1/9);%Eq.31 - fracture radius
    vradial=mscl*4/9*((Eprime*Q0^3)/(muprime*t^5))^(1/9);%diff(L)/diff(t)

    toughviscflag=0;

    %'Constants' and V
    muprime=eta*12;
    Eprime=E/(1-nu^2);
    mu=E/(2*(1+nu));
    %tmatch=(((Eprime*muprime^8)/Q0^6)*((mscl*mu*pi^2)/(delta_gamma^2*4*(1-nu)))^9)^(1/14);
    tmatch=(((mscl*pi^2)/(delta_gamma^2*4))^9*((2*mu^10*muprime^8)/((1-nu)^10*Q0^6)))^(1/14);
    tmatch2=(((mscl*pi^2)/(delta_gamma^2*8))^9*((Eprime^10*muprime^8)/(Q0^6)))^(1/14);
    test=1;
%     %Using solver to find the time where the two velocites match... (vtilde==v)
%     fun=@(t)tm(t,E,nu,eta,delta_gamma,Q0);
%     [tmatch,~]=fminsearch(fun,tmk);


else
    %% Toughness dominated using Tada book:
    kscl=(9/(pi^2*2))^(1/5); %Eq.36
    L=kscl*((Eprime^2*Q0^2*t^2)/(Kprime^2))^(1/5);%Eq.32 - fracture radius
    vradial=kscl*2/5*((Q0^2*Eprime^2)/(t^3*Kprime^2))^(1/5);%diff(L)/diff(t)
    vradial2=2/5*(9/(pi^2*2)*(Q0^2*Eprime^2)/(t^3*Kprime^2))^(1/5);%diff(L)/diff(t)
    %Now finding the time where the two velocites match... (vtilde==v)
    %t at which vtilde==v
    %tt=((9*8^3*27^5*E^2*eta^5*mu^5*pi^9)/((80^5*Kc^2*Q0^3*delta_gamma^10*(1-nu^2)^2*(1-nu)^5)))^(1/8)
    %simplified form:
    %tmatch=1.4558*((Eprime^2*eta^5*mu^5*pi^9)/((Kc^2*Q0^3*delta_gamma^10*(1-nu)^5)))^(1/8);
    tmatch=1.7313*((eta^5*mu^7*pi^9)/((Kc^2*Q0^3*delta_gamma^10*(1-nu)^7)))^(1/8);

    tmatch2=((4*9*8^3*27^5*eta^5*Eprime^7*pi^9)/((2^7*80^5*Kc^2*Q0^3*delta_gamma^10)))^(1/8);
    tmatch2=0.9440*((eta^5*Eprime^7*pi^9)/((Kc^2*Q0^3*delta_gamma^10)))^(1/8);
    tmatch2=9/4*((3*pi^9)/(5^5)*(eta^5*Eprime^7)/(Kc^2*Q0^3*delta_gamma^10))^(1/8);
    
    toughviscflag=1;
end

end

function [res]=tm(t,E,nu,eta,delta_gamma,Q0)

%'Constants' and V
V=Q0*t;
muprime=eta*12;
Eprime=E/(1-nu^2);
mu=E/(2*(1+nu));

%Ascent: 
uV=((4*(1-nu))/(27*pi^2*mu))*((delta_gamma^2*V)/(eta));

%Radial:
mscl=0.6955; %Eq.35
v=mscl*4/9*((Eprime^(1/3)*Q0)/(muprime^(1/3)*t^(5/3)))^(1/3);%diff(L)/diff(t)

%Ascent: 
uV=((4*(1-nu))/(27*pi^2*mu))*((delta_gamma^2*V)/(eta));

%Zero when thes match
res=abs(v-uV);

%Implict eq
res=abs(((4*(1-nu)*delta_gamma^2)/(pi^2*mscl*mu)*((Q0^6*t^14)/(Eprime*muprime^8))^(1/9)-1));



end