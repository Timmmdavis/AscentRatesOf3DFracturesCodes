%Spence and Turcotte Fig.1
clear
close all
figure; hold on
h0=1;%ref width
l0=1;%ref len
A=pi*h0*l0;
g=1;%9.81;
deltarho=1;%1000
t=[5 10 20 50];
eta=1;

y=linspace(-1,6,10000);
step=(y(2)-y(1));
ybar=y/l0;
for i=1:numel(t)
tbar=(h0^2*g*deltarho*t(i))/(l0*eta);   %eq.4  - dmlss time
hbar=sqrt((1/tbar)+(ybar/tbar));        %eq.12 - dmlss width
lbar=((3*pi)/4)^(2/3)*tbar^(1/3)-1;     %eq.14 - dmlss height of 'tip'
htip=sqrt((1/tbar)+(lbar/tbar));        %eq.12 - dmless width of 'tip'
area=sum(hbar(ybar<lbar).*step)*2         %Area of fracture
plot(hbar(ybar<lbar),ybar(ybar<lbar));
scatter(htip,lbar,'.k')

eta=1;
at=A;
%Roper and Lister Eq.6.7 - dimensional
z2=((9.*at.^2.*deltarho.*g.*t)./(16.*eta)).^(1/3);
htip2=sqrt((eta.*z2)./(deltarho.*g.*t));
scatter(htip2,z2-1,'.b') %Roper and Lister do base of crack at 0

end
xlim([0 1])
ylim([-1 5.5])
WhiteFigure;
xlabel('$\overline{h}$',Interpreter='latex')
ylabel('$\overline{y}$',Interpreter='latex')

