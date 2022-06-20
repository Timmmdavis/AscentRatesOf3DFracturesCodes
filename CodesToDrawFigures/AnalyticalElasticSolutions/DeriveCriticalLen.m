
clear
close all
syms mu nu p c  V deltag

eq=V==((8*(1-nu))/(3*mu))*p*c^3; 

eq2=solve(eq,p);

K1p=(2/pi)*eq2*sqrt(pi*c);
K1g=-(4/(3*pi))*deltag*c*sqrt(pi*c);

eq3=K1p+K1g==0;
solve(eq3,c)

