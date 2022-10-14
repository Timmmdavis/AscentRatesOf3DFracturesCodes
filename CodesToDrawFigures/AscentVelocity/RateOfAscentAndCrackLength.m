function [h,rate,tr] = RateOfAscentAndCrackLength(V,t,delta_gamma,nu,mu,eta,Kc)
%UNTITLED Summary of this function goes here
%   %Inputs:
%   V=Volume of injected fluid
%   t=time since start 
%   delta_gamma=stress gradient
%   nu=Poisson's ratio
%   mu=shear modulus
%   eta=dynamic viscosity
%   Kc=Fracture toughness - only needed if you uncomment below...
%
%   %Outputs:
%   h - height from initial centre to the top of the crack
%   rate - velocity of the crack's upper tip
%   tr - time it takes fracture height to reach 2c

h=((9.*(1-nu ).*delta_gamma.^(3).*V.^(3).*t.^(2))./...
        (pi.^(4).*mu.*eta^(2))).^(1/6);  

rate=(((1-nu).*delta_gamma^(3).*V.^(3))./(9^2.*pi.^(4).*mu.*eta.^(2).*t.^(4))).^(1/6);

tr=((3^(2)*pi^(8)*mu^(5)*eta^(4))/((1-nu)^(5)*delta_gamma^(9)*V^(3))).^(1/4);


end