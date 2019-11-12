ls = 1.2;
H = 16.48;
L = 17.48;
T= 2.6;
kT = 3.6;
rho = 0.805;
kz = 0.179;
lambda = (1.95*0.179^2)/rho;
dt = 0.02;
% \tilde{C}(t) for 1D eigenmode kz = 0.179, C(t) extracted from MD
 C2 = kT/rho/(H-sin(2*kz*H)/2/kz) / (L * 2*T*2);
 C = sqrt(C2);
N1 = 45000; 
N2 = 300; % C(t) point number
Cg = zeros(1,N1);
rg = normrnd(0,1,1,N1);
g1 = exp(-lambda*dt);
g2 = sqrt(1-g1^2);
for i = 1:1:N1
	if i == 1
		Cg(i) = rg(1);
	else
		Cg(i) = g1*Cg(i-1) +  g2*rg(i);
	end
end
Cg = C*Cg;
figure(1)
p1 = subplot(1,2,1);
plot(linspace(0,(N1-1)*dt,N1),Cg,'or','linewidth',0.7);
%tl = title(p1,{'$Time\ series\ \tilde{A}(t)\ $'},'fontsize', 20);
%set(tl,'Interpreter','latex');
xlabel({'$t (\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',80,'Interpreter','latex');ylabel({'$\tilde{A}(t) $'},'fontsize',80,'Interpreter','latex');
ylim([-0.15 0.15]); xlim([0 1000]);
set(gca,'fontsize', 40);
p2 = subplot(1,2,2);
plot(linspace(0,(N2-1)*dt,N2),-lambda*linspace(0,(N2-1)*dt,N2),'-ob','linewidth',1.1);
%tl = title(p2,{'$Logarithm\ of\ \tilde{A}(t)$'},'fontsize', 20);
%set(tl,'Interpreter','latex');
xlabel({'$t(\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',80,'Interpreter','latex');ylabel({'$\ln(<\tilde{A}(0)\tilde{A}(t)>/<\tilde{A}(0)\tilde{A}(0)>)$'},'fontsize',80,'Interpreter','latex');
set(gca,'fontsize', 40);
