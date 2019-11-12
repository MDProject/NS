kz = 0.179;
file1 = 'D:\\WinSCP\\LJ_fluid\\InstantMagnitude\\0.179';
%file2 = 'D:\\WinSCP\\LJ_fluid\\eigenmode\\kx_0\\0.179';
file2 = 'D:\\WinSCP\\LJ_fluid\\BACKUP\\LJ_fluid\\DATA\\Correlation\\kx2_0\\0.179';
fp1 = fopen(file1,'r');
fp2 = fopen(file2,'r');
[A,N] = fscanf(fp1,'%f');

%%% norm of kz = 0.179 mode
ls = 1.2;
H = 16.48;
L = 17.48;
T = 2.6;
natom = 5115;
% the norm of v^2 in code is (H-sin(2*k*H)/2/k) * (L * T * 4)
%%%

for i = 1:1:N
    if mod(i,2) == 0
        Coef(fix(i/2)) = A(i);
    else
        time1(fix(i+1)/2) = A(i);
    end
end
Coef = Coef * H * L * T  * 8 / natom / (H-sin(2*kz*H)/2/kz) / (L * T * 4);
mean(Coef.^2)

[A,N] = fscanf(fp2,'%f');
for i = 1:1:N
    if mod(i,2) == 0
        Cor(fix(i/2)) = A(i);
    else
        time2(fix(i+1)/2) = A(i);
    end
end

figure(1)
p1 = subplot(1,2,1);
plot(time1,Coef,'ob','linewidth',1.1);
%tl = title(p1,{'$Projection\ coefficient\ A(t)$'},'fontsize', 5);
%set(tl,'Interpreter','latex');
xlabel({'$t (\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',80,'Interpreter','latex');ylabel({'$A(t) $'},'fontsize',80,'Interpreter','latex');
ylim([-0.15 0.15]); xlim([0 1000]);
set(gca,'fontsize', 40);
p2 = subplot(1,2,2);
lnCor = (log(Cor));
plot(time2,lnCor,'oc','linewidth',2.4);
%tl = title(p2,{'$Logarithm\ of\ autocorrelation\ \ln(<A(0)A(t)>/<A(0)A(0)>)$'},'fontsize', 6);
%set(tl,'Interpreter','latex');
xlabel({'$t(\sqrt{m\sigma^2/\varepsilon})$'},'fontsize',80,'Interpreter','latex');ylabel({'$\ln(<A(0)A(t)>/<A(0)A(0)>)$'},'fontsize',80,'Interpreter','latex');
set(gca,'fontsize', 40);


