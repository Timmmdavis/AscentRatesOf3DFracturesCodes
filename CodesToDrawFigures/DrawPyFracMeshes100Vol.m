%DrawPyFracMeshes
clear
close all

%Times(seconds)=[500., 12875., 25250., 37625., 50000.];

xlen=71;
ylen=491;

%Get to mesh dir...
%Change this is you change the top level directory name
FolderName='SpeedForResubmissionAdditionalData';

%Get the address of the current working directory
pathstring = pwd;      

%Splitting this into a cell array
if ispc
    parts = strsplit(mfilename('fullpath'), '\'); %Windows      
else
    parts = strsplit(mfilename('fullpath'), '/'); %Mac/Linux
end

%Finding the scripts root directory and its location (n)
if verLessThan('matlab', '9.1') %(below v2016b where 'contains' was introduced)
    [~,n] = find(~cellfun(@isempty,strfind(parts,FolderName)));
else %Cleaner version:
    [~,n] = find(contains(parts,FolderName));     
end

%Adding all folders to the path that are in this dir 
if ispc
    NewStr=strjoin(parts(1,1:n),'\'); %Windows  
else
    NewStr=strjoin(parts(1,1:n),'/'); %Mac/Linux
end

Str2=strcat(NewStr,"\PyFracSimulations\MeshData");
cd(Str2);
figure;hold on
%cmap2 = colormap_cpt('Ccool-warm2');
%colormap(cmap2)
colormap(flipud(gray))

%Constants
mu=8e9;
nu=0.25;
E=(2.*mu).*(1+nu);
delta_gamma=2000*9.81;
Kc=2e6;
eta=0.005;
[CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
VolumeIn=CriticalVolume*100;
%Injection depth
InjD=486.52119372446924; 
%My approx
c=((9*VolumeIn*mu)/(2^4*delta_gamma*(1-nu)))^(1/4);
R=c*2;
p=2*delta_gamma*c/3; %Davis 2020 GRL -> tail is closed K-=0
P0=p;
P1=delta_gamma;

%Weertman length
ch=(Kc/(delta_gamma*sqrt(pi)))^(2/3);
ch=((3*sqrt(pi)*Kc)/(8*delta_gamma))^(2/3);%3D - davis 2020


%Head vol - zero at base and Kc at penny SIDES
clat=((6*sqrt(pi)*(Kc))/(8*delta_gamma))^(2/3);
p=(2*delta_gamma*clat)/3;
[Volume] = PennyCrackConstantPressure_Volume(nu,E,p,clat);

Kprime=sqrt(2/pi)*Kc;
Eprime=(1/pi)*(E/(1-nu^2));
VolumeGara=0.408*(Kprime^(8/3))/(Eprime*delta_gamma^(5/3));


%For cross sections
Scl=20000;
[CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma);
%Convert to 2D area with approx...for cross sections:
[v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,eta);
[zGer,rateGer,heightGer,widthGer]=GaragashAndGermanovich2022Arxiv(Kc,E,nu,eta,delta_gamma,0,VolumeIn);
A=2*c*Dn; %Area we start with
A=A+((Kc^2*(1-nu))/(2*mu*delta_gamma)); %We remove this from the roper and lister func but! we use the whole tail vol in the sol
xprofile=0;
xmv=c*3;

data=load('PyFracLargeWaterShale_Time0p03e3.csv');
x=data(:,1);y=data(:,2);widths=data(:,3);
t=0.03e3;
time(1)=t;
xmvCross=xmv-c*1.4;
[plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmvCross,xprofile,Scl,R);
[cRL,zpnts,hpnts,zRL,hRL]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,mu,nu,eta,Kc,t);
steplength=gradient(zpnts);
AreaDnAn=sum(hpnts*2.*steplength,'omitnan')
AreaDnAn2=AreaDnAn-((Kc^2*(1-nu))/(2*mu*delta_gamma))
plot((hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
plot((-hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
scatter((hRL*Scl+xmvCross)/R,(zRL-InjD-cRL)/R,'.k');
[~,HeadVol(1)]=ReshapeAndDraw(x,y,widths,xlen,0,R,InjD,ch);
text(0/R,-(InjD*1.15)/R,'10 min','HorizontalAlignment', 'center','Interpreter','latex','color',[0 0 51]/255)
text(0/R,-(InjD*1.23)/R,strcat('$V$ = ',num2str(VolumeIn),'m$^3$'),'HorizontalAlignment', 'center',Interpreter='latex')

% %% Weertman profile for first one
% %Half height of weertman crack: Rivalta 2015 dyke review Eq.3
% c2d=(Kc/(delta_gamma*sqrt(pi)))^(2/3);
% zhead=linspace(-c2d,c2d,1000);
% % %Rivalta 2015 dyke review Eq.4 % HALF OPENING!
% hhead=(((1-nu)*Kc)/(2*mu)).*sqrt(c2d./pi).*sqrt(1-(zhead./c2d).^2).*(1+(zhead./c2d));
% plot((hhead*Scl+xmvCross)/R,(zhead+c2d-InjD-cRL)/R,'--k')
% plot((-hhead*Scl+xmvCross)/R,(zhead+c2d-InjD-cRL)/R,'--k')

data=load('PyFracLargeWaterShale_Time0p1e3.csv');
x=data(:,1);y=data(:,2);widths=data(:,3);
t=0.1e3;
time(2)=t;
xmvCross=xmv*2-c*1.5;
[plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmvCross,xprofile,Scl,R);
[cRL,zpnts,hpnts,zRL,hRL]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,mu,nu,eta,Kc,t);
plot((hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
plot((-hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
scatter((hRL*Scl+xmvCross)/R,(zRL-InjD-cRL)/R,'.k');
[~,HeadVol(2)]=ReshapeAndDraw(x,y,widths,xlen,xmv,R,InjD,ch);
text((xmv)/R,-(InjD*1.15)/R,'6 h','HorizontalAlignment', 'center','Interpreter','latex','color',[0 0 153]/255)

data=load('PyFracLargeWaterShale_Time0p6e3.csv');
x=data(:,1);y=data(:,2);widths=data(:,3);
t=0.6e3;
time(3)=t;
xmvCross=xmv*3-c*1.5;
[plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmvCross,xprofile,Scl,R);
[cRL,zpnts,hpnts,zRL,hRL]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,mu,nu,eta,Kc,t);
plot((hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
plot((-hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
scatter((hRL*Scl+xmvCross)/R,(zRL-InjD-cRL)/R,'.k');
[~,HeadVol(3)]=ReshapeAndDraw(x,y,widths,xlen,xmv*2,R,InjD,ch);
text((xmv*2)/R,-(InjD*1.15)/R,'15 h','HorizontalAlignment', 'center','Interpreter','latex','color',[0 0 255]/255)

data=load('PyFracLargeWaterShale_Time1p3e3.csv');
x=data(:,1);y=data(:,2);widths=data(:,3);
t=1.3e3;
time(4)=t;
xmvCross=xmv*4-c*1.5;
[plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmvCross,xprofile,Scl,R);
[cRL,zpnts,hpnts,zRL,hRL]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,mu,nu,eta,Kc,t);
plot((hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
plot((-hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
scatter((hRL*Scl+xmvCross)/R,(zRL-InjD-cRL)/R,'.k');
[~,HeadVol(4)]=ReshapeAndDraw(x,y,widths,xlen,xmv*3,R,InjD,ch);
text((xmv*3)/R,-(InjD*1.15)/R,'30 h','HorizontalAlignment', 'center','Interpreter','latex','color',[102 102 255]/255)

data=load('PyFracLargeWaterShale_Time2p8e3.csv');
x=data(:,1);y=data(:,2);widths=data(:,3);
t=2.8e3;
time(5)=t;
xmvCross=xmv*5-c*1.5;
[plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmvCross,xprofile,Scl,R);
[cRL,zpnts,hpnts,zRL,hRL]=RoperAndListerConstantAreaSimilarity(A,delta_gamma,mu,nu,eta,Kc,t);
plot((hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
plot((-hpnts*Scl+xmvCross)/R,(zpnts-InjD-cRL)/R,'--k')
scatter((hRL*Scl+xmvCross)/R,(zRL-InjD-cRL)/R,'.k');
[~,HeadVol(5)]=ReshapeAndDraw(x,y,widths,xlen,xmv*4,R,InjD,ch);
text((xmv*4)/R,-(InjD*1.15)/R,'60 h','HorizontalAlignment', 'center','Interpreter','latex','color',[153 153 255]/255)

drawbrace([xmv*3.7 (max(zpnts)-InjD-cRL*4)]/R, [xmv*3.7 (max(zpnts)-InjD-cRL*2)]/R, 5, 'Color', 'k');
text((xmv*3.25)/R,(max(zpnts)-InjD-cRL*2.5)/R,{'Fracture', 'head'},'HorizontalAlignment', 'center',Interpreter='latex')

view(0,90)
cbar = colorbar( 'Location', 'NorthOutside' );
ylabel( cbar, 'Aperture $w$ (mm)', Interpreter='latex');

% cpos=get(c,'position');
% set(c,'position',[cpos(1) cpos(2) cpos(3) cpos(4)/2])

% plot([0 50],[50 50],'k')
% text(25,60,'50 m','HorizontalAlignment', 'center')

MveAxis=xmv*5;
plot([MveAxis MveAxis]/R,[0-InjD 300-InjD]/R,'k')
widths=3.5;
step=50;
plot([MveAxis-widths MveAxis+widths]/R,[-InjD -InjD]/R,'k')
%plot([MveAxis MveAxis+widths]/R,[-InjD+50 -InjD+50]/R,'k')
plot([MveAxis MveAxis+widths]/R,[-InjD+100 -InjD+100]/R,'k')
%plot([MveAxis-widths MveAxis+widths]/R,[-InjD+150 -InjD+150]/R,'k')
plot([MveAxis MveAxis+widths]/R,[-InjD+200 -InjD+200]/R,'k')
%plot([MveAxis-widths MveAxis+widths]/R,[-InjD+250 -InjD+250]/R,'k')
plot([MveAxis MveAxis+widths]/R,[-InjD+300 -InjD+300]/R,'k')
offset=5;
text((MveAxis+widths+offset)/R,-InjD/R,'0 m',Interpreter='latex')
%text((MveAxis+widths+offset)/R,(-InjD+50)/R,'50 m')
text((MveAxis+widths+offset)/R,(-InjD+100)/R,'100 m',Interpreter='latex')
%text((MveAxis+widths+offset)/R,(-InjD+150)/R,'150 m')
text((MveAxis+widths+offset)/R,(-InjD+200)/R,'200 m',Interpreter='latex')
%text((MveAxis+widths+offset)/R,(-InjD+250)/R,'250 m')
text((MveAxis+widths+offset)/R,(-InjD+300)/R,'300 m',Interpreter='latex')

step=R;
%plot([MveAxis-widths MveAxis+widths]/R,[-InjD -InjD]/R,'k')
%plot([MveAxis-widths MveAxis]/R,[-InjD+R -InjD+R]/R,'k')
plot([MveAxis-widths MveAxis]/R,[-InjD+2*R -InjD+2*R]/R,'k')
plot([MveAxis-widths MveAxis]/R,[-InjD+4*R -InjD+4*R]/R,'k')
plot([MveAxis-widths MveAxis]/R,[-InjD+8*R -InjD+8*R]/R,'k')
%plot([MveAxis-widths MveAxis+widths]/R,[-InjD+3*R -InjD+3*R]/R,'k')
offset=-30;
%text((MveAxis+widths+offset)/R,-InjD/R,'0 m')
%text((MveAxis+widths+offset)/R,(-InjD+R)/R,'2c')
text((MveAxis+widths+offset)/R,(-InjD+2*R)/R,'$4a$',Interpreter='latex')
text((MveAxis+widths+offset)/R,(-InjD+4*R)/R,'$8a$',Interpreter='latex')
text((MveAxis+widths+offset)/R,(-InjD+8*R)/R,'$16a$',Interpreter='latex')
%plot([MveAxis MveAxis]/R,[0-InjD 8*R-InjD]/R,'k')


axis('equal')
%ylim([-140/R 85/R])
WhiteFigure
axis off


%Draw analytical solution
smpl=1000;
[rx,rz]=meshgrid(linspace(-c,c,smpl),linspace(-c,c,smpl));
step=rz(2)-rz(1);
[th,r]=cart2pol(rx,rz);
rx(r>c)=NaN;
rz(r>c)=NaN;

% [c1,c2,k,k_prime,Kk,Ek,Q]=Assign_c1c2_andEllipticIntegrals(c,c);
% [w_P0,w_P1]=ComputeEllipticalCrackOpeningProfile(c,c,rx,rz,nu,E,k,k_prime,P1,P0,Q,c1,c2,Kk,Ek);
% aperture=(w_P0+w_P1).*2;

%[aperture]=AnalyticalPennyOpeningProfile(rx,rz,VolumeIn,delta_gamma,mu,nu,c,P0);
[Dngamma]=PennyCrackAssymGrad_disp(rx,rz,delta_gamma,mu,nu,c);
[Dnp]=PennyCrackConstPressure_disp(rx,rz,P0,mu,nu,c);
aperture=Dngamma+Dnp;

% [v,Dn,c]=AscentVelocityApproximation(VolumeIn,delta_gamma,mu,nu,0);
% aperture=ones(size(rx))*Dn;
% aperture(r>c)=NaN;

NumIntVol=sum((step*step).*aperture(:),'omitnan');
TadaV=((16*(1-nu^2))/(3*E))*p*c^3;
zheight=10;
mesh(rx/R,(rz+zheight)/R,aperture*1000,'FaceColor' ,'flat')
th=linspace(0,2*pi,100);
[xcirc,ycirc]=pol2cart(th,c);
plot3(xcirc/R,(ycirc+zheight)/R,ones(size(xcirc))+10,'k')
plot3(([0 sind(45)*R/2])/R,([0 sind(45)*-R/2]+zheight)/R,ones(2,1)+10,'k')
text((sind(45)*R/4)/R,(sind(45)*-R/4+zheight*1.7)/R,10,{'$a$'},'Interpreter','latex')
%plot([c*1.4 c*1.55]/R,[+zheight +zheight]/R,'k')
%text( (c*1.65)/R,+zheight/R,{'3D analytical','length scale'},Interpreter='latex')


%Draw analytical solution - cross section
smpl=1000;
rz=linspace(-c,c,smpl)';
rx=zeros(size(rz));
[th,r]=cart2pol(rx,rz);
rx(r>c)=NaN;
rz(r>c)=NaN;
[Dngamma]=PennyCrackAssymGrad_disp(rx,rz,delta_gamma,mu,nu,c);
[Dnp]=PennyCrackConstPressure_disp(rx,rz,P0,mu,nu,c);
aperture=Dngamma+Dnp;

ywall_grad=[rz;flipud(rz)]; %merging pos & neg side
xwall_grad=[aperture/2;flipud(-aperture/2)];

xmv=60;
greymid=0.6;
plt=patch(((xwall_grad*Scl)+xmv)/R,(ywall_grad+10)/R,[greymid,greymid,greymid],'LineStyle','none');%,'linestyle','--');

%Put in box
plot([-c*1.3 c*2.5]/R,[c*1.3+zheight c*1.3+zheight]/R,'k')
plot([-c*1.3 c*2.5]/R,[-c*1.3+zheight -c*1.3+zheight]/R,'k')
plot([c*2.5 c*2.5]/R,[-c*1.3+zheight c*1.3+zheight]/R,'k')
plot([-c*1.3 -c*1.3]/R,[c*1.3+zheight -c*1.3+zheight]/R,'k')
text( (-c*1.3)/R,(c*1.7+zheight)/R,{'a)'},Interpreter='latex')
text( (-c*1.3)/R,(-c*2.5)/R,{'b)'},Interpreter='latex')



caxis([0 Dn*2000])%2* Dn in mm
ylim([(-InjD-R/1.5)/R 1])

figure;hold on
plot(time*0.000277778,HeadVol)
hline(CriticalVolume,'k')
xlabel('time hrs')
ylabel('head vol m^3')
legend('numerical vol','V_c')

function [plt,HeadVol]=ReshapeAndDraw(x,y,widths,xlen,xmv,R,InjDepth,ch)

X=reshape(x,xlen,[])+xmv;
Y=reshape(y,xlen,[]);
Wgrd=reshape(widths,xlen,[]);
Wgrd=Wgrd*1000; %To mm!

flag=Wgrd~=0;
UpperTip=max(Y(flag));
flag2=Y>UpperTip-2*ch & Y<UpperTip;
HeadVol=sum(Wgrd(flag2))/1000;%Back to m's


%Draw the outside.
sz=(max(Y(:))-min(Y(:)))/1000; %m
[xnew,ynew] = meshgrid(min(Y(:)):sz:max(Y(:)),min(X(:)):sz:max(X(:)));
znew = interp2(Y,X,Wgrd,xnew,ynew, 'linear');
val=0.01;
contour3(ynew/R,xnew/R,znew+10,[val val]+10,'k'); 

Wgrd(Wgrd==0)=NaN;
plt=mesh(X/R,Y/R,Wgrd,'FaceColor' ,'flat');
scatter3(xmv/R,-InjDepth/R,10,'.k')

end


function [plt]=ReshapeAndDrawCrossSection(x,y,widths,xlen,xmv,xprofile,Scl,R)
    X=reshape(x,xlen,[]);
    Y=reshape(y,xlen,[]);
    Wgrd=reshape(widths,xlen,[]);
    Wgrd=Wgrd;%*1000; %To mm!
    
    [~,idx]=min(abs(X(:,1)-xprofile));
    
    %Stratified
    zStrat=Y(idx,:)';
    DnStrat=Wgrd(idx,:)';

    if DnStrat(1)==0 %Remove blank points start and end
        for i=1:numel(DnStrat)
            if DnStrat(i)~=0
                value=i-1;
                break
            end
        end
        DnStrat=DnStrat(value:end);
        zStrat=zStrat(value:end);
    end
    DnStrat=flipud(DnStrat);
    zStrat=flipud(zStrat);
    if DnStrat(1)==0 %Remove blank points start and end
        for i=1:numel(DnStrat)
            if DnStrat(i)~=0
                value=i-1;
                break
            end
        end
        DnStrat=DnStrat(value:end);
        zStrat=zStrat(value:end);
    end
    DnStrat=flipud(DnStrat);
    zStrat=flipud(zStrat);    

    steplength=zStrat(2)-zStrat(1);
    AreaDnNum=sum(DnStrat.*steplength,'omitnan')

    ywall_grad=[zStrat;flipud(zStrat)]; %merging pos & neg side
    xwall_grad=[DnStrat/2;flipud(-DnStrat/2)];
    
    greymid=0.6;
    plt=patch(((xwall_grad*Scl)+xmv)/R,(ywall_grad)/R,[greymid,greymid,greymid],'LineStyle','none');%,'linestyle','--');
    
    %plt=plot(Wgrd(idx,:),Y(idx,:))
%     WhiteFigure;
%     xlabel("Opening")
%     ylabel('depth')
%     xlim([-0.5 0.5])
%     ylim([-0.5e3 0.5e3])


end

