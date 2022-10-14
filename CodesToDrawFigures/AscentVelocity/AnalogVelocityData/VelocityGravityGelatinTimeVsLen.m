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

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M2-1837\data_8Sep22.txt"));
        Time=Time*24;
        Tank=1807; Tankheight=22.30;
        E=845;
        E=1345;
        Eerror=0;%ECOSOL-rgd std Table E.1
        Vol=20;
        delta_gamma=260;
        K1=37.13;
        Kc=K1/0.92;%Kc=0;disp('Guess!')
        eta=1.74e-3;
        rho_f=860;
        rho_r=1120;
        %title('M2 Silicon oil: Experiment 1933',Interpreter='latex')
        ylow=1e-5;
        yhigh=1e-1;
        deltattimelaps=1  ;     
        offset=1;        
        InjectionTime=(7-0.25)*24; %Seconds from video
       

    elseif i==2

        %subplot(1,3,2);hold on

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M50-1945\data.txt"));

        %M50 1945 - 
        Tank=1909; Tankheight=21.92;
        E=306;
        E=426;
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
        loadtime=NaN; %no load placed on exp
        InjectionTime=38-25; %Seconds from video


    elseif i==3
        %subplot(1,3,3);hold on

        [Time, x, z, vx, vz, Velocity] = importDelphineData2(strcat(CurrentDir,"\M1000-1967\data"));

        %M1000 - 1967 -
        Tank=1916; Tankheight=21.53;
        E=595;
        E=944;
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
        loadtime=1514-33;%seconds - 25min14sec in vid - injection at 33 secs
        InjectionTime=126-34; %Seconds from video

end

    for jj=1:10 %arbitary large no
        %Remove rows - where fracture jumps back a step
        good=diff(z)>0;
        good=logical([good; 1]); 
        z=z(good);
        Time=Time(good);
        Velocity=Velocity(good);
    end

    %Kc using the formula in Davis 2020 - critical vol - sec 3.1
    mu=E/(2*(1+nu));
    En=6.66e-4*mu;
    Kc=sqrt(En*(E/(1-nu^2)));

    V=Vol;
    t=Time;
    %Units
    z=z/100; %cm->m
    Velocity=Velocity/100;%cm/s->m/s
    V=V/1e6;%millilitres -> litres then litres -> m^3
    delta_gamma=(rho_r-rho_f)*9.81;
    %Adjust z
    z=(z-min(z));
    t=t-min(t);
    %extraang from video - ang from vert
    ExAng=10;
    Tankheight=Tankheight/100; %cm->m

    Velocity=(z./Time);

    [v,Dn,c]=AscentVelocityApproximation(V,delta_gamma,mu,nu,eta);
    R=c*2;


    plot(t,z,'color',col,'LineWidth',0.7)%,'-.'

   % If you want to see when the injection ended was placed on the tank...
    if ~isnan(InjectionTime)
        [~,idx]=min(abs(Time-InjectionTime));
        h=scatter(t(idx),z(idx),10,'^k','filled');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';

    end

end


ttri=logspace(0,5,100);
vtri=((3e-5*ttri).^(1/3));
plot(ttri,vtri,'-.k');

ttri=logspace(2.4,3.4,100);
vtri=((3e-5*ttri).^(1/3));
plot([ttri(1),ttri(end)],[vtri(end) vtri(end)],'k');
plot([ttri(1),ttri(1)],[vtri(1) vtri(end)],'k');
plot(ttri,vtri,'k');
text(min(ttri)*0.8,mean(vtri),'$1$',Interpreter='latex',HorizontalAlignment='center')
text(mean(ttri)*0.85,max(vtri)*1.2,'$3$',Interpreter='latex',HorizontalAlignment='center')


tt=xlabel('Time since start, $t$ (s)','interpreter','latex');
%t.Color = 'red';
ylabel('Fracture height, $h$ (m)','interpreter','latex')

WhiteFigure;

text(2e0,1.7e-4,{'$\uparrow$', 'Injection'},'color',c3,'Interpreter','latex')
text(5e3,1.1e-1,{'$\uparrow$','Hits top','of tank'},'color',c3,'Interpreter','latex',HorizontalAlignment='center')

leg=legend({'M2 oil,  \hspace{0.6em}\hspace{0.5em}  Exp 1837','M50 oil, \quad Exp 1945','M1000 oil,  Exp 1967','$t^{1/3}$'},'Interpreter','latex','location','NW');

xl=[10^0 10^4];
yl=[10^-2 10^0];
xlim(xl)
ylim(yl)

yticks([0.01 0.1 1])
yticklabels({'0.01','0.1','1'})
axis square