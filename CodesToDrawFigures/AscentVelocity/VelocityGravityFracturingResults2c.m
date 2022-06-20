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

%Same for all
nu=0.25;
Kc=2e6;
g=9.81;
c1=[0, 128, 255]/255;
c2=[0, 128, 255]/255;
c3=[0, 0.75, 0.75];
c4=[0 0 0];
c5=[0.6350, 0.0780, 0.1840];

cc1=[0 0 51]/255;
cc2=[0 0 153]/255;
cc3=[0 0 255]/255;
cc4=[102 102 255]/255;
cc5=[153 153 255]/255;

figure;hold on
xlabel('crack length [m]')
ylabel('velocity [m/s]')

for j=[2]%2:3
    if j==1
        col=c1;
    elseif j==2
        col=c2;
    elseif j==3
        col=c3;   
    elseif j==4
        col=c4;   
    elseif j==5
        col=c5;
    elseif j==9
        col=c5;
    end

    Case=j; %works for first 5 (below are below eq)

    if Case==2
        %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC] = importCSVResults2('PyFracResult2_Long.csv');
       
    elseif Case==3
        %Case III
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult3.csv');
        rhon=2000;
        eta=0.005;
        E=20e9;
        InjectedMass=750;
        injectionrate=0.015; %m3SECS
        injectiontime=25; %SECS    
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC] = importCSVResults2('PyFracResult3_Long.csv');
       
   
    end
    [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( E,nu );
    %Comp vars
    delta_gamma=rhon*g;
    VolumeIn=injectionrate*injectiontime;%/rhon;
    Eprime=E/(1-nu^2);
    [ ~,E,lambda,nu,G ] = ElasticConstantsCheck( E,nu );

   
    FractureHeightPy=FractureHeightPySRC;
    TimeTakenToReachEnd=Times(end)/3600;%hours

    %Remove start of Inj 
    strt=20;
    FractureHeightPy=FractureHeightPy(strt:end);
    Times=Times(strt:end);
    VelocityPy=VelocityPy(strt:end);

    [v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,eta);
    R=c*2;
    plot(FractureHeightPy/R,VelocityPy,'color',col,'LineWidth',0.5)%,'-.'
    %plot(FractureHeight/R,Velocity,':','marker','.','color',col);

    %Time to get to 6c:
    [~,idx]=min(abs(FractureHeightPy-(c*6)))
    Times(idx)

    [~,idx]=min(abs(FractureHeightPy-R));
    %hh=scatter(FractureHeightPy(idx)/R,VelocityPy(idx),15,'ok','filled');
    %hh.Annotation.LegendInformation.IconDisplayStyle = 'off';

    injectiontime=VolumeIn/injectionrate;
    [~,idx]=min(abs(Times-injectiontime));
    hh=scatter(FractureHeightPy(idx)/R,VelocityPy(idx),18,'k','^','filled');
    hh.Annotation.LegendInformation.IconDisplayStyle = 'off';

    yline(v,'--','color','k')%,'LineWidth',1) %v*1.35
    %Decay rate: Roper and Lister (2007) LargeFractureToughness:
    m=G/(1-nu);
    t=linspace(0,1e5,1e6);
    at=(2*Dn*c);

    t=Times;

    V=VolumeIn;

    [z,rate,tr] = RateOfAscentAndCrackLength(V,t,delta_gamma,nu,mu,eta,Kc);

%     z(t<tr)=NaN;
%     rate(t<tr)=NaN;

    %Adjust so rate at 2c is 1
    [~,idx1]=min(abs(z-c*2));
    %[~,idx]=min(abs(FractureHeightPy-c*2));
    idx1=1;

    %Simplified
    V=VolumeIn;pis=pi;
   

    [~,idx2]=min(abs(rate-min(VelocityPy))); %Cut off results at end of sim
    idx2=numel(rate)%idx2*2



    

    
    if Case==2 
        string=["10 min";  "6 h"; "15 h"; "30 h" ; "60 h" ];
        testtimes=[600 21600 54000 108000 216000];
        for jj=1:numel(testtimes)

            if jj==1

                col2=cc1;
                elseif jj==2

                col2=cc2;
                elseif jj==3

                col2=cc3;
                elseif jj==4

                col2=cc4;
                elseif jj==5

                col2=cc5;                
            end

            %Find 500s - Py results
            [~,idx]=min(abs(Times-testtimes(jj)));
            hh=scatter(FractureHeightPy(idx)/R,VelocityPy(idx),[],col2,'*');
            hh.Annotation.LegendInformation.IconDisplayStyle = 'off';
            text((FractureHeightPy(idx)-17)/R,VelocityPy(idx)*1.45,string(jj),'color',col2,Interpreter='latex')
        end
        %text(38/R,(10^-1)*1.06,strcat([char(8592),' End of Injection']))
        text(30/R,(10^-1)*1.0,'$\leftarrow$ End of Injection','color','k',Interpreter='latex')
        %text(42/R,0.0069,strcat(['2c',char(8594)]))
    end

    %Plot: adjusting to only show rates after 2c & v
    plot(z(idx1:idx2)./R,rate(idx1:idx2),'-.','color','k');
    %plot(z(1:end)./R,rate(1:end),'k');
    %plot(FractureHeightPy(idx:idx2)./R,rate(idx:idx2),'-.','color',col);



end
t=xlabel('Normalised distance to upper tip, $h/2a$','interpreter','latex');
%t.Color = 'red';
ylabel('Speed of upper tip, $v$ (m/s)','interpreter','latex')
set(gca, 'YScale', 'log')
WhiteFigure;

xlim([0 5.25])
ylim([10^-4 10^0])

leg=legend({'$v$, numerical','$\tilde{v}_V$, analytical','$v_V(t)$, analytical'},'Interpreter','latex')%,'III PyFrac, V=0.375 m$^3$','III Salimzadeh et al','Max ascent velocity','Analytical velocity'},'Interpreter','latex');
%leg.Title.String = 'Case';
%title('Modelled crack ascent velocities')

pbaspect([1 1 1])


