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
    %eta=eta*10;
    mu=E/(2*(1+nu));
    En=6.66e-4*mu;
    Kc=sqrt(En*(E/(1-nu^2)));

    V=Vol;
    t=Time;%/100;
    %Units
    z=z/100; %cm->m
    Velocity=Velocity/100;%cm/s->m/s
    V=V*0.001*0.001;%/1e6;%millilitres -> litres then litres -> m^3
    delta_gamma=(rho_r-rho_f)*9.81;

    %Adjust z and t so zero is start...
    z=(z-min(z));
    t=t-min(t);
    %extraang from video - ang from vert
    ExAng=10;
    Tankheight=Tankheight/100; %cm->m

    %Calc velocity
    Velocity2=(z./t);

    [CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
    [v,Dn,c]=AscentVelocityApproximation(V,delta_gamma,mu,nu,eta);
    R=2*c;
    
    lw=1.5;

    %V=V-CriticalVolume
    
    [zAN,rate,tr] = RateOfAscentAndCrackLength(V,t,delta_gamma,nu,mu,eta,Kc);


%     %Garagash and Germanovich Aug 2022 Arxiv - Page 16
%     Kprime=sqrt(2/pi)*Kc;
%     Eprime=(1/pi)*(E/(1-nu^2));
%     muprime=pi^2*eta;
%     Vg=0.408*(Kprime^(8/3)/(Eprime*(delta_gamma)^(5/3)));
%     zAN2=2.2*(((V-Vg).^2.*delta_gamma.^(7/3).*t)./(muprime*Kprime^(4/3))).^(1/3);
%     zAN=zAN2;

    %Link so both same z at time t (when z=R)....
    %This avoids the injection rate....
    Rscl=1.5;
    [~,idx]=min(abs(t-InjectionTime*2));
    diffx=zAN(idx)-z(idx);
    zAN=zAN-diffx;

    plot(z./R,zAN./R,'color',col,'LineWidth',0.7)%,'-.
    h=scatter(z(idx)./R,zAN(idx)./R,10,'sk','filled');%,'-.
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

end

xlabel('(Analog experiment height)$/(2a)$','interpreter','latex')
ylabel('Analytical height, $h/(2a)$','interpreter','latex')
WhiteFigure;

lim=7.5;

%hline(1,'--k')
xx=logspace(-5,2,1e4);
h=plot(xx,xx,'k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=logspace(-5,2,1e4);
h=plot(xx,xx/2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=logspace(-5,2,1e4);
h=plot(xx,xx*2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


text(lim*0.65,lim*0.62,'One to one','Interpreter','latex','rotation',45)

strt=lim*0.78;
text(strt,strt/2.2,'Prediction 2x shorter','Interpreter','latex','rotation',27.5, 'horizontalAlignment', 'center')
text(strt/2.2,strt,'Prediction 2x longer','Interpreter','latex','rotation',90-27.5, 'horizontalAlignment', 'center')

text(2.9,5.3e0,{'$\uparrow$', 'Hits top','of tank '},'color',c2, 'horizontalAlignment', 'left',Interpreter='latex')

axis('equal')
xlim([0 lim])
ylim([0 lim])


leg=legend({'M2:  \hspace{0.5em}\hspace{0.5em}  Exp 1933','M50: \quad Exp 1945','M1000:  Exp 1967'},'Interpreter','latex','location','NE');



