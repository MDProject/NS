clear all;
KMIN=0.1; % 0.225   0.45    0.7
KMAX=0.8; % 0.301 0.497   0.738
dk=0.001;
NUM=(KMAX-KMIN)/dk+1;
base='D:\\WinSCP\\LJ_fluid\\eigenmode\\kx2_0\\'; 
savePath='D:\\WinSCP\\LJ_fluid\\eigenmode\\kx2_0\\tau\\';
ifPlot=1;
ifDecayTime=0;
time=[];
for i=1:1:NUM
    k=KMIN+(i-1)*dk;
    path=[base,num2str(k,'%.3f')];
    fp=fopen(path,'r');
    [A,count]=fscanf(fp,'%f');
    arrayLength=count/2;
    for j=1:1:arrayLength
        if i==1
            time(j)=A((j-1)*2+1);
        end
        logCorrelationVX(i,j)=log(A((j-1)*2+2));
    end
    fclose(fp);
end
% plot decay time tau versus kz      0.232710566932577	0.465421133865155	0.698131700797732
kz=linspace(KMIN,KMAX,NUM);
startIdx=5;endIdx=20;    % 5 20(1)
for i=1:1:NUM
    x=time(startIdx:endIdx);
    y=logCorrelationVX(i,startIdx:endIdx);
    coef=polyfit(x,y,1);
    tau(i)=-1/coef(1);
end
idx=1;
for n=2:1:NUM-1
    if tau(n)>=tau(n-1)&&tau(n)>=tau(n+1)
        tauMax(idx)=tau(n);
        kz_eigen(idx)= KMIN+(n-1)*dk;
        idx=idx+1;
    end
end
k1 = 0.180;
tau1 = 1.347209507704577e+01;
k2 = 0.355;
tau2 = 3.397867804595660e+00;
k3 = 0.536;
tau3 = 1.380239155467137e+00;
figure('visible','on');
plot(kz,tau,'or','linewidth',2.2);
hold on;
plot(k1,tau1,'.k',k2,tau2,'.k',k3,tau3,'.k','MarkerSize' ,56);
hold on;
plot(kz,0.427./kz./kz,'-b','linewidth',2);
set(gca,'fontsize', 50);xlabel({'$Eigen\ wavevector\ k_z\ (in\ unit\ of\ \frac{1}{\sigma})$'},'fontsize',50,'Interpreter','latex');ylabel({'$Decay\ time\ \tau (in unit of \sqrt{m\sigma}$'},'fontsize',50,'Interpreter','latex');
% [0.177870375545877,0.356718667515213,0.537149444520186,0.719312927334858]

k1=0.178;
k2=0.357;
k3=0.5365;%0.714;
HMIN=14.5;
HMAX=20.5;
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
p1 = subplot(1,2,2);
plot(H,Int1,'-or',H,Int2,'-ob',H,Int3,'-oc','linewidth',0.8);
hold on;
plot(zeros(HNUM,1)+16.48,linspace(-10,10,HNUM),'-k','linewidth',2.2);
hold on;
ylim([-3 3]);
plot(16.48,0,'.k','MarkerSize' ,56);
%lgd = legend('$Mode1\ k_z = 0.179$','$Mode2\ k_z = 0.355$','$Mode3\ k_z = 0.536$','$HD Boundary$','Location','southeast');
%set(lgd,'Interpreter','latex');
tl = title(''); %Eigenmodes\ Mutual\ Orthogonality
set(tl,'Interpreter','latex');
xlabel({'$z\ (\sigma)$'},'fontsize',70,'Interpreter','latex');ylabel({'$\int v^(1)_z(z)v^(2)_z(z)dz\ $'},'fontsize',70,'Interpreter','latex');
set(gca,'fontsize', 40);
p2 = subplot(1,2,1);
plot(kz,tau,'or','linewidth',2.2);
hold on;
plot(k1,tau1,'.k',k2,tau2,'.k',k3,tau3,'.k','MarkerSize' ,56);
hold on;
plot(kz,0.427./kz./kz,'-b','linewidth',2);
xlabel({'$\ k_z\ (1/\sigma)$'},'fontsize',70,'Interpreter','latex');ylabel({'$Decay\ time\ \tau (\sqrt{m\sigma}$'},'fontsize',70,'Interpreter','latex');
set(gca,'fontsize', 40);
% [0.177870375545877,0.356718667515213,0.537149444520186,0.719312927334858]