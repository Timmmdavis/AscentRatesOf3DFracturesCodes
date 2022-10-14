clear; close all
%Plot the upper tip velocity of cracks from:
%GRL: Gravity Hydraulic Fracturing: A Method to Create Selfâ€?Driven Fractures
%and those computed using PyFrac. 
%Note the actual data is not in the supplementary material of the GRL
%paper. I asked the authors for this and used numerical differention to get
%the front velocity. 


M=0;
pfluid=0;
laminar=1;

%Same for all
g=9.81;
c1=[0, 0.4470, 0.7410];c1=[150, 150, 150]./256;
c2=[0.4660, 0.6740, 0.1880];c2=[255, 51, 51]./256;
c3=[0, 0.75, 0.75];
c4=[0 0 0];
c5=[0.6350, 0.0780, 0.1840];

VVplot=figure;hold on
HHplot=figure;hold on

VLpIn=[2,10,100];%[2,10,100]
vp=[];
colin=0;
for q=1:numel(VLpIn)
    VLp=VLpIn(q);
for j=[1:3]%2:3

    if j==1
        col=c1;
        lw=1;
    elseif j==2
        col=c2;
        lw=0.7;
    elseif j==3
        col=c3;
        lw=0.7;
    end
    colin=colin+1;

    Case=j; %works for first 5 (below are below eq)
    if Case==1
        % Oil in Gelatin 
        nu = 0.5-1e-9 ;                         % Poisson's ratio       [dmlss] 
        G= 276  ;                               % Shear modulus         [pa] 
        E = (2*G)*(1+nu)    ;                   % Young's modulus       [pa] 
        Eprime = E / (1 - nu ^ 2) ;             % plain strain modulus  [pa]
        Kc = 19;                                % fracture toughness    [pa.sqrt(m)]
        fluiddensity=1000-160;                  % [kg/m3]
        rockdensity=1000  ;                     % [kg/m3]
        deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
        eta=48e-3;                              % fluidviscosity [pa.s] - water=~1.1e-3

        if VLp==2
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('OilGelatin2V.csv');
        elseif VLp==10
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('OilGelatin10V.csv');
        elseif VLp==100
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('OilGelatin100V.csv');
        end

    elseif Case==2
        % Basaltic Dyke
        nu = 0.25 ;                          % Poisson's ratio       [dmlss] 
        G= 25e9  ;                           % Shear modulus         [pa] 
        E = (2*G)*(1+nu)    ;                % Young's modulus       [pa] 
        Eprime = E / (1 - nu ^ 2) ;          % plain strain modulus  [pa]
        Kc = 6*1e6;                          % fracture toughness    [pa.sqrt(m)]
        fluiddensity=2950 ;                  % [kg/m3]
        rockdensity=3000  ;                  % [kg/m3]
        deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
        eta=20;                              % fluidviscosity [pa.s] - water=~1.1e-3

        if VLp==2
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('Magma2V.csv');
        elseif VLp==10
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('Magma10V.csv');
        elseif VLp==100
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('Magma100V.csv');
        end


    elseif Case==3
        %Case III
        % Water in rock properties
        nu = 0.25 ;                         % Poisson's ratio       [dmlss] 
        G= 8e9  ;                           % Shear modulus         [pa] 
        E = (2*G)*(1+nu)    ;               % Young's modulus       [pa] 
        Eprime = E / (1 - nu ^ 2) ;         % plain strain modulus  [pa]
        Kc = 2*1e6;                         % fracture toughness    [pa.sqrt(m)]
        fluiddensity=1000 ;                 % [kg/m3]
        rockdensity=3000  ;                 % [kg/m3]
        deltagamma=(rockdensity-fluiddensity)*9.81; % gradient in weight [pa.m^-1]
        eta=0.005;                          % fluidviscosity [pa.s] - water=~1.1e-3

        if VLp==2
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('WaterShale2V.csv');
        elseif VLp==10
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('WaterShale10V.csv');
        elseif VLp==100
            [FractureHeightPy, Times, VelocityPy,FractureHeightPyFromSrc] = importCSVResults2('WaterShale100V.csv');
        end

    end


    %Desired height from src (comment if you want from base of crack)...
    FractureHeightPy=FractureHeightPyFromSrc;

    %% Smooth
    VelocityPy(3:end) = movmean(VelocityPy(3:end),40); 
    
 
    %Get the elastic constants
    [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( E,nu );

    %Compute the critical volume required to ascend
    [CriticalVolume] = CriticalVolumeDavis2020(nu,G,Kc,deltagamma);

    if VLp==2
        VolumeIn=CriticalVolume*2;
    elseif VLp==10
        VolumeIn=CriticalVolume*10;
    elseif VLp==100
        VolumeIn=CriticalVolume*100;
    else 
        error('Check this')
    end
   
    [ ~,E,lambda,nu,G ] = ElasticConstantsCheck( E,nu );

    TimeTakenToReachEnd=Times(end)/3600;%hours

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Some analytical insights %%%%%%%%%%%%%%%%%%%%%
    
    %Length scales
    [v,Dn,c]=AscentVelocityApproximation(VolumeIn,deltagamma,mu,nu,eta);
    R=c*2;

    %Length and rate of injection (based on volume and parameters above)
    timeFor2c=(R)/v;
    InjectionRate=VolumeIn/timeFor2c;%volume/lengthofinjection_s                    [m/s] 

    %Davis et al 2020 GRL - critical vol - head volume
    [CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,deltagamma);
    V=VolumeIn;%-CriticalVolume;    
    
    pis=pi;

    %Substitute vals:
    t=Times;
    
    %Rate of ascent -analytical
    [z,rate,tr] = RateOfAscentAndCrackLength(V,t,deltagamma,nu,mu,eta,Kc);

%     %Garagash and Germanovich Aug 2022 Arxiv - Page 16
%     Kprime=sqrt(2/pi)*Kc;
%     Eprime=(1/pi)*(E/(1-nu^2));
%     muprime=pi^2*eta;
%     Vg=0.408*(Kprime^(8/3)/(Eprime*(deltagamma)^(5/3)));
%     z2=2.2*(((V-Vg).^2.*deltagamma.^(7/3).*t)./(muprime*Kprime^(4/3))).^(1/3);
%     rate2=gradient(z2)./gradient(t);
%     z=z2;rate=rate2;

    %So the vector is correctly plotted
    FractureHeightPy(end)=NaN;
    z(end)=NaN;


    figure(VVplot)
   if VLp==2
        h=plot(VelocityPy/v,rate/v,'color',col,'LineWidth',lw);%,'-.'
    elseif VLp==10
        h=plot(VelocityPy/v,rate/v,'color',col,'LineWidth',lw,LineStyle='--');%,'-.'
    elseif VLp==100
        h=plot(VelocityPy/v,rate/v,'color',col,'LineWidth',lw,LineStyle='-.');%,'-.'
    end
    %Turn off in legend
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    vp=[vp max(VelocityPy)];
    text(0.15,0.2,'$\leftarrow$ Increasing time','rotation',50,'Interpreter','latex')
    

    figure(HHplot)
    if VLp==2
        h=plot(FractureHeightPy/R,z/R,'color',col,'LineWidth',lw);%,'-.'
    elseif VLp==10
        h=plot(FractureHeightPy/R,z/R,'color',col,'LineWidth',lw,LineStyle='--');%,'-.'
    elseif VLp==100
        h=plot(FractureHeightPy/R,z/R,'color',col,'LineWidth',lw,LineStyle='-.');%,'-.'
    end
    %Turn off in legend
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %text(2.5,4,strcat(['increasing time',char(8594)]),'rotation',45)
    text(2.56,3.7,'Increasing time $\rightarrow$','rotation',45,'Interpreter','latex')
    
end
end

xlabel('(PyFrac height)$/(2a)$','interpreter','latex')
ylabel('Analytical height, $h/(2a)$','interpreter','latex')
%set(gca, 'YScale', 'log')
WhiteFigure;

%Maximum length shown
maxL=8;

axis('equal')

%hline(1,'--k')
xx=linspace(0,10,1e4);
h=plot(xx,xx,'k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

h=plot(xx,xx/2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

h=plot(xx,xx*2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=linspace(0,1.5,5);
yy=ones(size(xx))*1.5;
h=plot(xx,yy,':k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
xx=linspace(0,1.5,5);
yy=ones(size(xx))*1.5;
h=plot(yy,xx,':k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([0 maxL])
ylim([0 maxL])

xticks([0:1:maxL])
yticks([0:1:maxL])

text(0.25*maxL*2,0.26*maxL*2,'One to one','Interpreter','latex','rotation',45)
text(0.3*maxL*2,0.16*maxL*2,'Analytical 2x shorter','Interpreter','latex','rotation',27)
text(0.14*maxL*2,0.3*maxL*2,'Analytical 2x longer','Interpreter','latex','rotation',19+45)


text(0.025*maxL*2,0.475*maxL*2,'Oil in gelatine','Interpreter','latex',Color=c1);
text(0.025*maxL*2,0.425*maxL*2,'Basaltic dyke','Interpreter','latex',Color=c2);
text(0.025*maxL*2,0.45*maxL*2,'Water in shale','Interpreter','latex',Color=c3);

plot([-1 -1],[-1 -1],'k')
plot([-1 -1],[-1 -1],'k--')
plot([-1 -1],[-1 -1],'k-.')
leg=legend({'$2 V_c$','$10 V_c$','$100 V_c$'},'Interpreter','latex','location','SE');

figure(VVplot)

xlabel('(PyFrac speed)$/\tilde{v}_V$','interpreter','latex')
ylabel('Analytical speed, $v_V(t)/\tilde{v}_V$','interpreter','latex')
WhiteFigure;



axis('equal')

xx=linspace(0,0.5,1e4);
h=plot(xx,xx,'k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=linspace(0,0.5,1e4);
h=plot(xx,xx/2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xx=linspace(0,0.5,1e4);
h=plot(xx,xx*2,'--k');%,'AutoUpdate','off')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


xlim([0 0.5])
ylim([0 0.5])

xticks([0:0.1:0.5])
yticks([0:0.1:0.5])

text(0.25,0.26,'One to one','Interpreter','latex','rotation',45)

text(0.3,0.16,'Analytical 2x slower','Interpreter','latex','rotation',27)
text(0.14,0.3,'Analytical 2x faster','Interpreter','latex','rotation',19+45)


text(0.025,0.475,'Oil in gelatine','Interpreter','latex',Color=c1);
text(0.025,0.425,'Basaltic dyke','Interpreter','latex',Color=c2);
text(0.025,0.45,'Water in shale','Interpreter','latex',Color=c3);

plot([-1 -1],[-1 -1],'k')
plot([-1 -1],[-1 -1],'k--')
plot([-1 -1],[-1 -1],'k-.')
leg=legend({'$2 V_c$','$10 V_c$','$100 V_c$'},'Interpreter','latex','location','SE');





