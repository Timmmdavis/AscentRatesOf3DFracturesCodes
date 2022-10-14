%Toughness scaling:
clear
close all
%syms Eprime Q0 t muprime Kprime

% Case 2 - Fig.1 
eta=0.05;
Kc=2e6;
E=20e9;
nu=0.25;
Q0=0.015;
delta_gamma=1000*9.81;
mu=E/(2*(1+nu));
Vinjected=1.95;

nu=0.1;



% %% Oil in Gelatin  - PyFrac sims - Fig 3/4
% nu = 0.5-1e-9 ;                          % Poisson's ratio       [dmlss] 
% G= 276  ;                           % Shear modulus         [pa] 
% E = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
% Eprime = E / (1 - nu ^ 2) ;% plain strain modulus  [pa]
% Kc = 19;                        % fracture toughness    [pa.sqrt(m)]
% fluiddensity=1000-160;                  % [kg/m3]
% rockdensity=1000  ;                  % [kg/m3]
% delta_gamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
% eta=48e-3;                            % fluidviscosity [pa.s] - water=~1.1e-3
% clear
% 
% %% Basaltic Dyke
% nu = 0.25 ;                          % Poisson's ratio       [dmlss] 
% G= 25e9  ;                           % Shear modulus         [pa] 
% E = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
% Eprime = E / (1 - nu ^ 2) ;% plain strain modulus  [pa]
% Kc = 6*1e6;                        % fracture toughness    [pa.sqrt(m)]
% fluiddensity=2950 ;                  % [kg/m3]
% rockdensity=3000  ;                  % [kg/m3]
% delta_gamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
% eta=20;                            % fluidviscosity [pa.s] - water=~1.1e-3
% clear
% 
% %% Water in rock properties
% nu = 0.25 ;                          % Poisson's ratio       [dmlss] 
% G= 8e9  ;                           % Shear modulus         [pa] 
% E = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
% Eprime = E / (1 - nu ^ 2) ;% plain strain modulus  [pa]
% Kc = 2*1e6;                        % fracture toughness    [pa.sqrt(m)]
% fluiddensity=1000 ;                  % [kg/m3]
% rockdensity=3000  ;                  % [kg/m3]
% delta_gamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
% eta=0.005;                            % fluidviscosity [pa.s] - water=~1.1e-3


% mu=E/(2*(1+nu));
% [Vc] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
% Vinjected=Vc*100;
% [vtilde,Dn,c]=AscentVelocityApproximation(Vinjected,delta_gamma,mu,nu,eta);
% Q0=(Vinjected*vtilde)/c;

%Mori and Lecampion
tkk=(Kc^(8/3))/((E/(1-nu^2))*Q0*delta_gamma^(5/3));
%Davis etc
tr=(((3^2*pi^8*mu^5)/((1-nu^5)))*((eta^4)/(delta_gamma^9*Vinjected^3)))^(1/4);

t=Vinjected/Q0;%Not defining yet...
[vradial,vtilde,tmatch,toughviscflag,tmk] = RadialVsAscentFrontVelocity(Kc,nu,E,eta,t,Q0,delta_gamma);

disp('t/tmk')
t/tmk
disp('t/tmatch')
t/tmatch

%Same order of magnitude
n=floor( log10(t/tmatch));
if t<tmatch || n==0
    disp('The fracture will be radial at end of injection')
    [vradial,vtilde,tmatch2,toughviscflag,tmk] = RadialVsAscentFrontVelocity(Kc,nu,E,eta,tmatch,Q0,delta_gamma);
    if tmatch2~=tmatch %Changed regime...visc is now toughness
        [vradial,vtilde,tmatch,toughviscflag,tmk] = RadialVsAscentFrontVelocity(Kc,nu,E,eta,tmatch2,Q0,delta_gamma);
    end
    [vtilde,Dn,c]=AscentVelocityApproximation(Q0*tmatch,delta_gamma,mu,nu,eta);
    if toughviscflag==0
        disp('The fracture is still in the viscous regime at the end of injection')
    else
        disp('The fracture is in a toughness regime by the end of injection')
    end
end


t=linspace(0,tmatch,100);

Kfrontspeed=zeros(size(t));
Ascentspeed=zeros(size(t));
Ratio=zeros(size(t));
tvflag=zeros(size(t));
for i=1:numel(t)

    [vradial,vtilde,tmatch,tflag,tmk] = RadialVsAscentFrontVelocity(Kc,nu,E,eta,t(i),Q0,delta_gamma);

    Kfrontspeed(i)=vradial;
    Ascentspeed(i)=vtilde;

    a=((9*(Q0*t(i))*mu)/(16*delta_gamma*(1-nu)))^(1/4);
    Ratio(i)=vradial./vtilde;%(V*uV)/a;
    tvflag(i)=tflag;
end

figure;
hold on
plot(t,Kfrontspeed./Ascentspeed)
%plot(t,ratio)
ylim([0 2])
xlabel('time (s)')
ylabel('v_r/v_z')

