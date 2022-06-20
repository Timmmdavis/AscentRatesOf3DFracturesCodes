% Water in rock properties
nu = 0.25 ;                          % Poisson's ratio       [dmlss] 
G= 8e9  ;                           % Shear modulus         [pa] 
youngs_mod = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
Eprime = youngs_mod / (1 - nu ^ 2) ;% plain strain modulus  [pa]
Kc = 2*1e6;                        % fracture toughness    [pa.sqrt(m)]
fluiddensity=1000 ;                  % [kg/m3]
rockdensity=3000  ;                  % [kg/m3]
deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
eta=0.005;                            % fluidviscosity [pa.s] - water=~1.1e-3

% Basaltic Dyke
nu = 0.25 ;                          % Poisson's ratio       [dmlss] 
G= 25e9  ;                           % Shear modulus         [pa] 
youngs_mod = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
Eprime = youngs_mod / (1 - nu ^ 2) ;% plain strain modulus  [pa]
Kc = 6*1e6;                        % fracture toughness    [pa.sqrt(m)]
fluiddensity=2950 ;                  % [kg/m3]
rockdensity=3000  ;                  % [kg/m3]
deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
eta=20;                            % fluidviscosity [pa.s] - water=~1.1e-3

% Oil in Gelatin 
nu = 0.5 ;                          % Poisson's ratio       [dmlss] 
G= 276  ;                           % Shear modulus         [pa] 
youngs_mod = (2*G)*(1+nu)    ;       % Young's modulus       [pa] 
Eprime = youngs_mod / (1 - nu ^ 2) ;% plain strain modulus  [pa]
Kc = 19;                        % fracture toughness    [pa.sqrt(m)]
fluiddensity=1000-160;                  % [kg/m3]
rockdensity=1000  ;                  % [kg/m3]
deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
eta=48e-3;                            % fluidviscosity [pa.s] - water=~1.1e-3


[CriticalVolume] = CriticalVolumeDavis2020(nu,G,Kc,deltagamma);
VolumeIn=CriticalVolume*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Some analytical insights %%%%%%%%%%%%%%%%%%%%%

%Max velocity 
Max_v=(4/(27*eta*G*pi^2))*VolumeIn*deltagamma^2*(1-nu);
%Volume in radius - Davis
VolFractureRadius=((9*G*VolumeIn)/(16*deltagamma*(1-nu)))^(1/4)   
timeFor2c=(VolFractureRadius*2)/Max_v;

%Length and rate of injection (based on volume and parameters above)
rate=VolumeIn/timeFor2c%volume/lengthofinjection_s                    [m/s] 

%% Compute dmlss toughness
mu=G;
[v,Dn,c]=AscentVelocityApproximation(VolumeIn,deltagamma,mu,nu,eta);
A=(2*c*Dn);
%Roper and Lister 2.10:
m=mu/(1-nu);
%My height from centre
y0=c*100;
%y0=((9.*VolumeIn.^(3).*delta_gamma.^(3).*t.^(2).*(1-nu ))./...
%        (eta^(2).*pi.^(4).*mu)).^(1/6);  
%Rate - Spnce and turcotte
Q=(9*deltagamma*A^3)/(32*eta*y0^3);
Kdmlss=(2*Kc^4)/(3*eta*Q*m^3)


