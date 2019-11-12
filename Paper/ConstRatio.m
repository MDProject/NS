filename1 = 'D:\\WinSCP\\LJ_fluid\\InstantMagnitude\\New folder\\Result_L_250.mat';
load(filename1,'lambdai_vec_decrease','L','H');
tau0 = 0.025;

Larray = [250,300,400,500];
beta_max = [10571,12471,16626,20781];
plot(Larray,beta_max,'.r','linewidth',3,'MarkerSize' ,60);
hold on;
plot(Larray,(20781-16626)/100.*Larray+6,'-b','linewidth',3);
set(gca,'fontsize', 40);
xlabel({'$L (\sigma)$'},'fontsize',50,'Interpreter','latex');ylabel({'$\beta_{max}$'},'fontsize',50,'Interpreter','latex');
legend('numerical values','asymptotic line');


rho = 0.805;
eta = 1.95;
R = rho/eta;
tau0 = 0.025;
ratio_theory = R/8/tau0/pi;
ratio = [0.6415,0.6305,0.6305,0.6305];
Larray = [250,300,400,500];
plot(Larray,ratio,'.r','MarkerSize' ,70);
hold on;
line([200,600],[ratio_theory,ratio_theory],'linestyle','--','linewidth',3);
legend('numerical values','theoretical value');
xlim([200 600]);ylim([0.5 0.7]);
set(gca,'fontsize', 40);
xlabel({'$L (\sigma)$'},'fontsize',50,'Interpreter','latex');ylabel({'$M$'},'fontsize',50,'Interpreter','latex');



