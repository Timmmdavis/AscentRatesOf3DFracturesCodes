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
c1=[0, 0, 0]/255;
c2=[0, 128, 255]/255;
c3=[0, 0.75, 0.75];
%c4=[0 0 0];
c4=[0, 0, 0]/255;
c5=[0, 128, 255]/255;
c6=[0, 0.75, 0.75];

cc1=[0 0 51]/255;
cc2=[0 0 153]/255;
cc3=[0 0 255]/255;
cc4=[102 102 255]/255;
cc5=[153 153 255]/255;



figure;hold on
xlabel('crack length [m]')
ylabel('velocity [m/s]')

for j=[1:6]%2:3
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
    elseif j==6
        col=c6;
    end

    Case=j; %works for first 5 (below are below eq)

    if Case==1
        %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_Long.csv');
        meshsz=30;

    elseif Case==2
         %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_LongNewImplicit_Grid30NewFilter.csv');
        meshsz=30;    

    elseif Case==3
         %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_LongNewExplicit_Grid30NewFilter.csv');
        meshsz=50;        

    elseif Case==4
        %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_LongNewPredCorr_Grid50NewFilter.csv');
        meshsz=50;
       
    elseif Case==5
        %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_LongNewImplicit_Grid50NewFilter.csv');
        meshsz=50;
       
    elseif Case==6
         %Case II
        [FractureHeight, Time, Velocity] = importCSVResults('SalimzadehResult2.csv');
        rhon=1000;
        eta=0.05;
        E=20e9;
        InjectedMass=1950;
        injectionrate=0.015; %m3SECS
        injectiontime=130; %SECS
        [FractureHeightPy, Times, VelocityPy,FractureHeightPySRC,xcoord] = importCSVResults3('PyFracResult2_LongNewExplicit_Grid50NewFilter.csv');
        meshsz=50;

        %Remove where scheme starts to break down...
        FractureHeightPy=FractureHeightPy(1:end-20);
        Times=Times(1:end-20);
        VelocityPy=VelocityPy(1:end-20);
        FractureHeightPySRC=FractureHeightPySRC(1:end-20);



%         %Remove where scheme starts to break down...
%         FractureHeightPy=FractureHeightPy(1:end-20);
%         Times=Times(1:end-20);
%         VelocityPy=VelocityPy(1:end-20);
%         FractureHeightPySRC=FractureHeightPySRC(1:end-20);        

    end
    [ K,E,lambda,nu,mu ] = ElasticConstantsCheck( E,nu );
    %Comp vars
    delta_gamma=rhon*g;
    VolumeIn=injectionrate*injectiontime;%/rhon;
    Eprime=E/(1-nu^2);
    [ ~,E,lambda,nu,G ] = ElasticConstantsCheck( E,nu );

%     n=6;
%     Times=Times(1:n:end);
%     FractureHeightPySRC=FractureHeightPySRC(1:n:end);
%     VelocityPy=VelocityPy(1:n:end);
   
%     %Take average of every two results
%     n=numel(VelocityPy);
%     if floor(n/2)==n/2
%         % code for even
%         VelocityPy=(VelocityPy(1:2:end-1)+VelocityPy(2:2:end))/2;
%         Times=(Times(1:2:end)+Times(2:2:end))/2;
%         FractureHeightPy=(FractureHeightPy(1:2:end)+FractureHeightPy(2:2:end))/2;
%         FractureHeightPySRC=(FractureHeightPySRC(1:2:end)+FractureHeightPySRC(2:2:end))/2;
%     else 
%       % code for odd
%         VelocityPy=(VelocityPy(1:2:end-1-1)+VelocityPy(2:2:end))/2;
%         Times=(Times(1:2:end-1)+Times(2:2:end))/2;
%         FractureHeightPy=(FractureHeightPy(1:2:end-1)+FractureHeightPy(2:2:end))/2;
%         FractureHeightPySRC=(FractureHeightPySRC(1:2:end-1)+FractureHeightPySRC(2:2:end))/2;
%     end
%     %Calc gradient
%     VelocityPy=gradient(FractureHeightPy)./gradient(Times);
%     %FInd upticks in the velocity and remove these
%     Badflag=diff(VelocityPy);
%     Badflag=[0;Badflag];
%     Badflag=Badflag<0;
%     %Removing them
%     VelocityPy=VelocityPy(Badflag);
%     FractureHeightPySRC=FractureHeightPySRC(Badflag);
%     FractureHeightPy=FractureHeightPy(Badflag);
%     Times=Times(Badflag);
% %     Badflag=abs(xcoord)>0;
% %     %Removing them
% %     VelocityPy=VelocityPy(Badflag);
% %     FractureHeightPySRC=FractureHeightPySRC(Badflag);
% %     FractureHeightPy=FractureHeightPy(Badflag);
% %     Times=Times(Badflag);

    VelocityPy(3:end) = movmean(VelocityPy(3:end),40);
    
    
    FractureHeightPy=FractureHeightPySRC;
    TimeTakenToReachEnd=Times(end)/3600;%hours

%     %Remove start of Inj 
%     strt=20;
%     FractureHeightPy=FractureHeightPy(strt:end);
%     Times=Times(strt:end);
%     VelocityPy=VelocityPy(strt:end);

    [v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,eta);
    R=c*2;
    if j>3
    plot(Times,VelocityPy,'color',col,'LineWidth',0.5)%,'-.'
    else
    plot(Times,VelocityPy,'--','color',col,'LineWidth',0.5)%
    end
    %scatter(Times,VelocityPy,'.k');
%     scatter(Times,VelocityPy,5,xcoord,'filled');
%     colorbar
%     clim([mean(xcoord)-1 mean(xcoord)+3])
    %plot(FractureHeight/R,Velocity,':','marker','.','color',col);

    %Time to get to 6c:
    [~,idx]=min(abs(FractureHeightPy-(c*6)))
    Times(idx)

    [~,idx]=min(abs(FractureHeightPy-R));
    %hh=scatter(FractureHeightPy(idx)/R,VelocityPy(idx),15,'ok','filled');
    %hh.Annotation.LegendInformation.IconDisplayStyle = 'off';

    injectiontime=VolumeIn/injectionrate;
    [~,idx]=min(abs(Times-injectiontime));
    hh=scatter(injectiontime,VelocityPy(idx),18,'k','^','filled');
    hh.Annotation.LegendInformation.IconDisplayStyle = 'off';

    %yline(v,'--','color','k')%,'LineWidth',1)
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
        string=["";  "6 h"; "15 h"; "30 h" ; "60 h" ];
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
            hh=scatter(t(idx),VelocityPy(idx),[],col2,'*');
            hh.Annotation.LegendInformation.IconDisplayStyle = 'off';
            text((t(idx)-17),VelocityPy(idx)*1.45,string(jj),'color',col2,Interpreter='latex')
        end
        %text(38/R,(10^-1)*1.06,strcat([char(8592),' End of Injection']))
        %%text(injectiontime*1.1,(10^-1)*1.0,'$\leftarrow$ End of Injection','color','k',Interpreter='latex')
        %text(42/R,0.0069,strcat(['2c',char(8594)]))
    end


 



end


ttri=logspace(4.6,5.5,100);
vtri=((0.4*ttri).^(-2/3));
%[z,vtri,tr] = RateOfAscentAndCrackLength(V,ttri,delta_gamma,nu,mu,eta,Kc);
plot(ttri,vtri,'-.k');

ttri=logspace(4.9,5.3,100);
vtri=((0.4*ttri).^(-2/3));
%[z,vtri,tr] = RateOfAscentAndCrackLength(V,ttri,delta_gamma,nu,mu,eta,Kc);
%plot(ttri,vtri,'k');
plot([ttri(1),ttri(end)],[vtri(end) vtri(end)],'k');
plot([ttri(1),ttri(1)],[vtri(1) vtri(end)],'k');
plot(ttri,vtri,'-k');
%text(mean(ttri),mean(vtri)*1.1,'$t^{-2/3}$',Interpreter='latex',rotation=-30,HorizontalAlignment='center')
text(min(ttri)*0.9,mean(vtri)*1.,'$2$',Interpreter='latex',HorizontalAlignment='center')
text(mean(ttri)*.9,min(vtri)*0.86,'$3$',Interpreter='latex',HorizontalAlignment='center')
 

tt=xlabel('Time since start, $t$ (s)','interpreter','latex');
%t.Color = 'red';
ylabel('Speed of upper tip, $v$ (m/s)','interpreter','latex')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
WhiteFigure;
% 
% xlim([10^2 10^6])
% ylim([10^-4 10^0])

xl=[10^3 10^6];
yl=[10^-4 10^-1];
decades_equal(gca,xl,yl);   

leg=legend({'Numerical: $v$','Analytical, $\tilde{v}_V$','Analytical: $v_V(t)$'},'Interpreter','latex')%,'III PyFrac, V=0.375 m$^3$','III Salimzadeh et al','Max ascent velocity','Analytical velocity'},'Interpreter','latex');
leg=legend({'$v$: Predictor 30$\cdot$30','$v$: Implicit 30$\cdot$30','$v$: Explicit 30$\cdot$30','$v$: Predictor 50$\cdot$50','$v$: Implicit 50$\cdot$50','$v$: Explicit 50$\cdot$50','$t^{-2/3}$'},'Interpreter','latex')%,'III PyFrac, V=0.375 m$^3$','III Salimzadeh et al','Max ascent velocity','Analytical velocity'},'Interpreter','latex');
%leg.Title.String = 'Case';
%title('Modelled crack ascent velocities')

%axis('equal')


