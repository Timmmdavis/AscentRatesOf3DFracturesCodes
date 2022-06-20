clear; close all
%Plot the upper tip velocity of cracks from:
%GRL: Gravity Hydraulic Fracturing: A Method to Create Self-Driven Fractures
%and those computed using PyFrac. 
%Note the actual data is not in the supplementary material of the GRL
%paper. I asked the authors for this and used numerical differention to get
%the front velocity. 

M=0;
pfluid=0;
laminar=1;

clear
close all
nu=0.5;

%visc = Pa.s
%delta_gamma = kg.m3
%E = pa
c1=[255, 128, 0]/255;
c2=[0., 102, 51]/255;
c3=[255, 0., 127]/255;
c4=[0 0 0];
c5=[0.6350, 0.0780, 0.1840];
figure;
hold on

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

CurrentDir=pwd();

colin=0;
for i=1:3

    if i==1
        col=c1;
    elseif i==2
        col=c2;
    elseif i==3
        col=c3;   
    elseif i==4
        col=c4;   
    elseif i==5
        col=c5;
    elseif i==9
        col=c5;
    end
    colin=colin+1;

    if i==1
        %subplot(1,3,1);hold on

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M2-1933\data.txt"));
        %M2 1933 - 
        E=1995;
        Eerror=0;%ECOSOL-rgd std Table E.1
        Vol=40;
        delta_gamma=260;
        K1=52.61;
        Kc=NaN;
        eta=1.74e-3;
        rho_f=860;
        rho_r=1120;
        %title('M2 Silicon oil: Experiment 1933',Interpreter='latex')
        ylow=1e-5;
        yhigh=1e-1;
        deltattimelaps=1  ;     
        offset=1;
        

    elseif i==2

        %subplot(1,3,2);hold on

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M50-1945\data.txt"));

        %M50 1945 - 
        
        E=306;
        Eerror=61;%ECOSOL-rgd std Table E.1
        Vol=10;
        delta_gamma=160;
        K1=11.38;
        Kc=K1/0.67;
        eta=48e-3;
        rho_f=960;
        rho_r=1120;
        %title('M50 Silicon oil: Experiment 1945',Interpreter='latex')
        ylow=1e-5;
        yhigh=1e-1;
        deltattimelaps=1;
        offset=0.5;


    elseif i==3
        %subplot(1,3,3);hold on

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M1000-1967\data"));

        %M1000 - 1967 -
        
        E=595;
        Eerror=19;%ECOSOL-rgd std Table E.1
        Vol=10;
        delta_gamma=150;
        K1=15.25;
        Kc=K1/0.64;
        eta=970e-3;
        rho_f=970;
        rho_r=1120;
        %title('M1000 Silicon oil: Experiment 1967',Interpreter='latex')
        ylow=1e-5;
        yhigh=1e-1;
        deltattimelaps=1;
        offset=0.3;

    end


    for jj=1:1000 %arbitary large no
        %Remove rows - where fracture jumps back a step
        good=gradient(z)>0;
        z=z(good);
        Time=Time(good);
        Velocity=Velocity(good);
    end

    mu=E/(2*(1+nu));
    V=Vol;
    t=Time*100;
    %Units
    z=z/100; %cm->m
    Velocity=Velocity/100;%cm/s->m/s
    V=V/1e6;%millilitres -> litres then litres -> m^3
    delta_gamma=(rho_r-rho_f)*9.81;
    %Adjust z
    z=(z-min(z)); 
    %extraang from video - ang from vert
    ExAng=10;

    Velocity=(z./Time);
%     %Remove start of Inj 
%     strt=20;
%     FractureHeightPy=FractureHeightPy(strt:end);
%     Times=Times(strt:end);
%     VelocityPy=VelocityPy(strt:end);

    [v,Dn,c]=AscentVelocityApproximation(V,delta_gamma,mu,nu,eta);
    R=c*2;


    plot(t,Velocity,'color',col,'LineWidth',0.5)%,'-.'







 



end


ttri=logspace(1,7,100);
vtri=((1*ttri).^(-2/3));
plot(ttri,vtri,'-.k');

ttri=logspace(2.1,2.8,100);
vtri=((1*ttri).^(-2/3));
%plot(ttri,vtri,'k');
plot([ttri(1),ttri(end)],[vtri(end) vtri(end)],'k');
plot([ttri(1),ttri(1)],[vtri(1) vtri(end)],'k');
plot(ttri,vtri,'k');
%text(mean(ttri),mean(vtri)*1.1,'$t^{-2/3}$',Interpreter='latex',rotation=-30,HorizontalAlignment='center')
text(min(ttri)*0.8,mean(vtri),'$2$',Interpreter='latex',HorizontalAlignment='center')
text(mean(ttri),min(vtri)*0.7,'$3$',Interpreter='latex',HorizontalAlignment='center')


tt=xlabel('Time since start, $t$ (s)','interpreter','latex');
%t.Color = 'red';
ylabel('Speed of upper tip, $v$ (m/s)','interpreter','latex')

WhiteFigure;

% xlim([7e3 3e5])
% ylim([1e-4 10^-2])

% leg=legend({'Numerical: $v$','Analytical: $\tilde{v}_V$','Analytical: $v_V(t)$'},'Interpreter','latex')%,'III PyFrac, V=0.375 m$^3$','III Salimzadeh et al','Max ascent velocity','Analytical velocity'},'Interpreter','latex');
% leg=legend({'Numerical: $v$','$(2t)^{-2/3}$'},'Interpreter','latex')%,'III PyFrac, V=0.375 m$^3$','III Salimzadeh et al','Max ascent velocity','Analytical velocity'},'Interpreter','latex');
% 
leg=legend({'M2 Silicon oil,  \hspace{0.6em}\hspace{0.5em}  Exp 1933','M50 Silicon oil, \quad Exp 1945','M1000 Silicon oil,  Exp 1967','$t^{-2/3}$'},'Interpreter','latex','location','SE');
%leg=legend({'M50 Silicon oil: \quad Exp 1945','M1000 Silicon oil:  Exp 1967','$t^{-2/3}$'},'Interpreter','latex','location','NE');

%leg.Title.String = 'Case';
%title('Modelled crack ascent velocities')


% xlim([10^1 10^7])
% ylim([10^-6 10^0])
% 
% axis('equal')

xl=[10^1 10^6];
yl=[10^-6 10^-1];
% xlim([10^1 10^7])
% ylim([10^-6 10^0])

decades_equal(gca,xl,yl);                      %# Make the decades equal sizes

