clear all;
L=17.48;
H=6;
% L=13.5;
% H=13.5;
ls=[1.2]; % only unknown 2.2
kx=pi/L;%4*pi/L;
n=3;  
% slip length dispersion relation       x = kz*H
ky=zeros(length(ls),n); % anti-symmetric
kz=zeros(length(ls),n); % symmetric modes
for k=1:1:length(ls)
    func=@(x) kx*H*tanh(kx*H)+x*tan(x)+ls(k)*H*(x^2/H/H+kx^2);
    for i=1:1:n
        x0=[i*pi];
        options = optimoptions('fsolve','MaxIterations',2000,'MaxFunctionEvaluations',500,'FunctionTolerance',10^(-40));
        x=fsolve(func,x0,options);
        ky(k,i)=x/H;
    end
end
for k = 1:1:length(ls)
    func=@(x) kx*H/tanh(kx*H)-x/tan(x)+ls(k)*H*(x^2/H/H+kx^2);
    for i=1:1:n
        x0=[(i+0.5-1/3)*pi];
        options = optimoptions('fsolve','MaxIterations',2000,'MaxFunctionEvaluations',500,'FunctionTolerance',10^(-40));
        x=fsolve(func,x0,options);
        kz(k,i)=x/H;
    end
end
% \rho = 1725/(13.5*13.5*1.5*8)=0.789
% (ls,n the root)
rho=0.805;
eta = 1.9634;
R=rho/eta;
tau=1./((kx^2+ky.^2)/R);

%% Plot velocity field
savePath='D:\\MatlabR2016a\\LJ_fluid\\WritePaper\\Figure\\Antisym';
x=linspace(-L,L,25);
y=linspace(-H,H,25);
[X,Y]=meshgrid(x,y);
iwR=kx^2+ky.^2;
% lambda=-kx*cos(ky*H)/iwR/cosh(kx*H);
% ***** Set Figure Format *****
figFormat = 'png';
figAppend=['.',figFormat];
% *****                                *****
for i=1:1:n
    figure('visible', 'off');
    Vx=vpa((-kx*cos(ky(i)*H)/(kx^2+ky(i).^2)/cosh(kx*H)*sinh(kx*Y)-ky(i)/(kx^2+ky(i).^2)*sin(ky(i)*Y)).*cos(kx*X));
    Vy=vpa((-kx*cos(ky(i)*H)/(kx^2+ky(i).^2)/cosh(kx*H)*cosh(kx*Y)+kx/(kx^2+ky(i).^2)*cos(ky(i)*Y)).*sin(kx*X));
    quiver(X,Y,Vx,Vy,'linewidth',1.1);set(gca,'fontsize', 20);
    %title('Anti-Symmetric Hydrodynamic Eigenmodes','fontsize',16);
    xlabel({'$Coordinate\ x\ (\mu m)$'},'fontsize',17,'Interpreter','latex');ylabel({'$Coordinate\ z\ (\mu m)$'},'fontsize',17,'Interpreter','latex');
    save=[savePath,'Analy_',num2str(double(iwR(i)),'%.5f'),figAppend];
    param = ['-d',figFormat];
    print(gcf, param, '-r1500' , save);
end


%% Plot velocity field
savePath='D:\\MatlabR2016a\\LJ_fluid\\WritePaper\\Figure\\Sym';
x=linspace(-L,L,25);
y=linspace(-H,H,25);
[X,Z]=meshgrid(x,y);
iwR=kx^2+kz.^2;
% lambda=-kx*cos(ky*H)/iwR/cosh(kx*H);
% ***** Set Figure Format *****
figFormat = 'png';
figAppend=['.',figFormat];
% *****                                *****
for i=1:1:n
    figure('visible', 'off');
    Vx = vpa((-kx/(kx^2+kz(i)^2)*sin(kz(i)*H)/sinh(kx*H)*cosh(kx*Z)+kz(i)/(kx^2+kz(i)^2)*cos(kz(i)*Z)).*sin(kx*X));
    Vz = vpa((kx/(kx^2+kz(i)^2)*sin(kz(i)*H)/sinh(kx*H)*sinh(kx*Z)-kx/(kx^2+kz(i)^2)*sin(kz(i)*Z)).*cos(kx*X));
    quiver(X,Z,Vx,Vz,'linewidth',1.1);set(gca,'fontsize', 20);
    %title('Symmetric Hydrodynamic Eigenmodes','fontsize',16);
    xlabel({'$Coordinate\ x\ (\mu m)$'},'fontsize',17,'Interpreter','latex');ylabel({'$Coordinate\ z\ (\mu m)$'},'fontsize',17,'Interpreter','latex');
    save=[savePath,'Analy_',num2str(double(iwR(i)),'%.5f'),figAppend];
    param = ['-d',figFormat];
    print(gcf, param, '-r1500' , save);
end