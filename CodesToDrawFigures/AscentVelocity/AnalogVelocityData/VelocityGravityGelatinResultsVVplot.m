%Mx =x has units cSt - centistoke - multiply by fluid density for visc

%     '1836' :  deltattimelaps=120
%     '1837' :  deltattimelaps=48
%     '1847' :  deltattimelaps=150
%     '1856' :  deltattimelaps=120
%     '1857' :  deltattimelaps=180

clear
close all
nu=0.5;

CurrentDir=pwd();


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
    t=Time;
    %Units
    z=z/100; %cm->m
    Velocity=Velocity/100;%cm/s->m/s
    V=V/1e6;%millilitres -> litres then litres -> m^3
    delta_gamma=(rho_r-rho_f)*9.81;
    %Adjust z
    z=(z-min(z)); 
    %extraang from video - ang from vert
    ExAng=10;

    %t=linspace(0,1e4,1e5);


    %Calc velocity
    Velocity2=(z./Time);
    %Velocity=smooth(Velocity,150,'sgolay',55);

    [CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
    [v,Dn,c]=AscentVelocityApproximation(V,delta_gamma,mu,nu,eta);
    R=2*c;
    %Turb?
    Re=(Dn*v*rho_f)/eta;
    disp(Re)
    if Re>1400
        %http://www.ecourses.ou.edu/cgi-bin/ebook.cgi?topic=fl&chap_sec=08.1&page=theory
        %fprintf('TurbFlow! red\n')
        fprintf(2,'\nTurbFlow!\n')
    end
    
    lw=1.5;

% %     %Adjust for tail volume (remove the head)
% %     V=V-CriticalVolume;
    
    [zAN,rate,tr] = RateOfAscentAndCrackLength(V,t,delta_gamma,nu,mu,eta,Kc);

    %Supplied rate
    %rate=((V.^(3).*delta_gamma.^(3).*(1-nu))./(81.*eta^(2).*pi.^(4).*t.^(4).*mu)).^(1/6);
    Vhigh=V+1/1e6;Ehigh=E-Eerror;muhigh=Ehigh/(2*(1+nu));etahigh=eta;
    [zANhigh,ratehigh,tr] = RateOfAscentAndCrackLength(Vhigh,t,delta_gamma,nu,muhigh,etahigh,Kc);
    Vlow=V-1/1e6;Elow=E+Eerror;mulow=Elow/(2*(1+nu));etalow=eta;delta_gammalow=delta_gamma;%*cosd(ExAng);
    [zANlow,ratelow,tr] = RateOfAscentAndCrackLength(Vlow,t,delta_gamma,nu,mulow,etalow,Kc);


    plot(Velocity2./v,rate./v,'color',col,'LineWidth',0.5)%,'-.'


end

xlabel('(Analog experiment velocity)$/\tilde{v}_V$','interpreter','latex')
ylabel('Analytical velocity, $v_V(t)/\tilde{v}_V$','interpreter','latex')
WhiteFigure;



%axis('equal')


%hline(1,'--k')
xx=logspace(-6,0,1e4);
h=plot(xx,xx,'k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=logspace(-6,1,1e4);
h=plot(xx,xx/2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=logspace(-6,1,1e4);
h=plot(xx,xx*2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

% xticks([0:0.1:0.5])
% yticks([0:0.1:0.5])

text(1,1,'One to one','Interpreter','latex','rotation',45)

text(0.00004,0.00016,'Prediction 2x slower','Interpreter','latex','rotation',45)
text(0.00015,0.00004,'Prediction 2x faster','Interpreter','latex','rotation',45)


% text(0.0002,1,strcat(['Injection',char(8594)]),'color',c2,'Interpreter','latex')
% %text(0.0003,0.1,strcat(['Hits top of tank',char(8595)]),'color',c2)
% text(0.018,0.01,strcat([char(8593),' Hits top of tank']),'color',c2,Interpreter='latex')

text(0.00011,1,'Injection $\rightarrow$','color',c2,'Interpreter','latex')
text(0.019,0.018,'$\leftarrow$ Hits top of tank','color',c2,Interpreter='latex')


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xl=[10^-5 10^2];
yl=[10^-5 10^2];
decades_equal(gca,xl,yl);    

% xticks([0:0.1:0.5])
% yticks([0:0.1:0.5])

leg=legend({'M2 Silicon oil:  \hspace{0.5em}\hspace{0.5em}  Exp 1933','M50 Silicon oil: \quad Exp 1945','M1000 Silicon oil:  Exp 1967'},'Interpreter','latex','location','SE');
%leg.Title.String = 'Experiment';    

%title('Gelatin analog vs analytical ascent velocities at time t')


