function [v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,eta)

c=((9*VolumeIn*mu)/(16*delta_gamma*(1-nu)))^(1/4);
Dn=VolumeIn/(pi*c^2);
v=((Dn^2)/(12*eta))*delta_gamma;
%v2=((4*(1-nu))/(27*pi^2*eta))*((delta_gamma^2*VolumeIn)/(eta))

end
