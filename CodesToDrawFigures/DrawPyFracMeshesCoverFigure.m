%DrawPyFracMeshes
clear
close all

%Times(seconds)=[500., 12875., 25250., 37625., 50000.];

xlen=71;
ylen=71;
xmv=65;

%Get to mesh dir...
%Change this is you change the top level directory name
FolderName='VelocityPaperUpdatedCodes';

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
VolumeIn=1.95;
delta_gamma=1000*9.81;
Kc=2e6;
%Injection depth
InjD=93.75; 
%My approx
c=((9*VolumeIn*mu)/(2^4*delta_gamma*(1-nu)))^(1/4);
R=c*2;
p=2*delta_gamma*c/3; %Davis 2020 GRL -> tail is closed K-=0
P0=p;
P1=delta_gamma;


load('PyFracSim2_Time600')
ReshapeAndDraw(x,y,widths,xlen,0,R)
%text(0/R,-125/R,'10 min','HorizontalAlignment', 'center')
%text(0/R,-140/R,'V = 1.95 m^3','HorizontalAlignment', 'center')

load('PyFracSim2_Time10800')
ReshapeAndDraw(x,y,widths,xlen,xmv,R)
%text((xmv)/R,-125/R,'3 h','HorizontalAlignment', 'center')

load('PyFracSim2_Time21600')
ReshapeAndDraw(x,y,widths,xlen,xmv*2,R)
%text((xmv*2)/R,-125/R,'6 h','HorizontalAlignment', 'center')

load('PyFracSim2_Time32400')
ReshapeAndDraw(x,y,widths,xlen,xmv*3,R)
%text((xmv*3)/R,-125/R,'9 h','HorizontalAlignment', 'center')

load('PyFracSim2_Time54000')
ReshapeAndDraw(x,y,widths,xlen,xmv*4,R)
%text((xmv*4)/R,-125/R,'15 h','HorizontalAlignment', 'center')


view(0,90)

set(gcf,'Position',[10 10 600 500])


axis('equal')
ylim([-140/R 85/R])
WhiteFigure
axis off


%Draw analytical solution
smpl=1000;
[rx,rz]=meshgrid(linspace(-c,c,smpl),linspace(-c,c,smpl));
step=rz(2)-rz(1);
[th,r]=cart2pol(rx,rz);
rx(r>c)=NaN;
rz(r>c)=NaN;

caxis([0 1.5])

axis('equal')
pos = get(gcf, 'Position');
pos(3)=pos(4)*1.2;
set(gcf, 'Position', pos) 

function [plt]=ReshapeAndDraw(x,y,widths,xlen,xmv,R)

X=reshape(x,xlen,[])+xmv;
Y=reshape(y,xlen,[]);
Wgrd=reshape(widths,xlen,[]);
Wgrd=Wgrd*1000; %To mm!

%Draw the outside.
sz=(max(Y(:))-min(Y(:)))/1000; %m
[xnew,ynew] = meshgrid(min(Y(:)):sz:max(Y(:)),min(X(:)):sz:max(X(:)));
znew = interp2(Y,X,Wgrd,xnew,ynew, 'linear');
val=0.01;
contour3(ynew/R,xnew/R,znew+10,[val val]+10,'k'); 

Wgrd(Wgrd==0)=NaN;
plt=mesh(X/R,Y/R,Wgrd,'FaceColor' ,'flat');
%scatter3(xmv/R,-93.75/R,10,'.k')

end

