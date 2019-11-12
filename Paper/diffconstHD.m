diffusion_z_10600 = diffusion(1,:,10600);
%plot(z0,diffusion_z_400.*rho./Rho ,'-oy',z0,diffusion_z_800.*rho./Rho ,'-ob',z0,diffusion_z_1600.*rho./Rho ,'-og',z0,diffusion_z_3200.*rho./Rho,'-*k',z0,diffusion_z_8786.*rho./Rho ,'-or','linewidth',2.2);
%legend('n=400','n=800','n=1600','n=3200','n=8786');
plot(z0(2:length(z0)),diffusion_z_10600(2:length(z0)).*rho./Rho(2:length(z0)) ,'.-r','linewidth',3.3,'MarkerSize',40);
hold on;
diffusion_z_md = [0.218,0.208,0.238,0.237,0.222,0.224,0.217,0.235,0.201];
z0_md = [0,2,4,6,8,10,12,14,16];
plot(z0_md,diffusion_z_md,'-db','linewidth',2.8,'MarkerSize',18);
leg = legend('HMs','MD');
pos=get(leg,'pos');
w=pos(3); % 原宽度
ratio=5; % 调整比例
pos(3)=w*ratio;
pos(1)=pos(1)-w*(ratio-1);
set(leg,'pos',pos)
%title('Diffusion constant distribution by HD modes composition');
ylim([0.1 0.3]);
set(gca,'fontsize', 40);xlabel({'$z\ (\sigma)$'},'fontsize',50,'Interpreter','latex');ylabel({'$D(z)\ (\sqrt{\varepsilon\sigma^2/m})$'},'fontsize',50,'Interpreter','latex');
% save('D:\\WinSCP\\LJ_fluid\\InstantMagnitude\\New folder\\Result_L_250','-v7.3');