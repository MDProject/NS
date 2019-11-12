clear all;close all; % 1D
k1=0.179;
k2=0.356;
k3=0.536;%0.714;
HMIN=15;
HMAX=20;
HNUM=100;
H=linspace(HMIN,HMAX,HNUM);
func1=@(z) sin(k1.*z).*sin(k2.*z);
func2=@(z) sin(k1.*z).*sin(k3.*z);
func3=@(z) sin(k2.*z).*sin(k3.*z);
for i=1:1:HNUM
    Int1(i)=integral(func1,-H(i),H(i),'AbsTol',1e-12);
    Int2(i)=integral(func2,-H(i),H(i),'AbsTol',1e-12);
    Int3(i)=integral(func3,-H(i),H(i),'AbsTol',1e-12);
end
figure(1)
plot(H,Int1,'-or',H,Int2,'-ob',H,Int3,'-oc','linewidth',0.8);
hold on;
plot(zeros(HNUM,1)+16.48,linspace(-10,10,HNUM),'-k','linewidth',2.2);
hold on;
plot(16.48,0,'.k','MarkerSize' ,56);
lgd = legend('$Mode1\ k_z = 0.179$','$Mode2\ k_z = 0.355$','$Mode3\ k_z = 0.536$','$HD Boundary$');
set(lgd,'Interpreter','latex');
tl = title('$Eigenmodes\ Mutual\ Orthogonality$');
set(tl,'Interpreter','latex');
xlabel({'$Channel\ coordinate\ z\ (in\ unit\ of\ \sigma)$'},'fontsize',50,'Interpreter','latex');ylabel({'$Eigenmodes\ mutual\ integral\ \int v^(1)_zv^(2)_zdz\ $'},'fontsize',50,'Interpreter','latex');
set(gca,'fontsize', 25);
