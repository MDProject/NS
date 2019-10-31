function [ ] = plot_slip_length(eigen_func,L,H,Nx,Nz,ls1,ls2,nidx)

DX = 2*L/(Nx-1);
DZ = 2*H/(Nz-1);

phi = eigen_func(:,nidx);
Phi = zeros(Nx,Nz);

for j = 2:1:Nz-1
    for i = 1:1:Nx-1
        Phi(i,j) = phi((j-2)*(Nx-1)+i);
    end
end
Phi(Nx,:) = Phi(1,:);
R1 = (DZ-2*ls1)./(DZ+2*ls1);
R2 = (DZ-2*ls2)./(DZ+2*ls2);
Phi_down = R1.*Phi(:,2);
Phi_up = R2.*Phi(:,Nz-1);

% \partial^2 vx / \partial z^2
% lower boundary, 2nd order
dvx_dz1 = (Phi_down+Phi(:,2)-2*Phi(:,1))/DZ/DZ;
dvx_dz2 = (Phi_up+Phi(:,Nz-1)-2*Phi(:,Nz))/DZ/DZ;
vx_1 = (Phi(:,2) - Phi_down)/2/DZ;
vx_2 = (Phi_up-Phi(:,Nz-1))/2/DZ;

ls1_fdm = vx_1./dvx_dz1;
ls2_fdm = -vx_2./dvx_dz2;

x = linspace(-L,L,Nx);
figure(1)
plot(x,ls1_fdm,'-or','linewidth',2.2);
hold on;
plot(x,ls1,'-.b','linewidth',2.2,'MarkerSize' ,26);
hold on;
set(gca,'fontsize', 40);
figure(2)
plot(x,ls2_fdm,'-or','linewidth',2.2);
hold on;
plot(x,ls2,'-.b','linewidth',2.2,'MarkerSize' ,26);
hold on;
set(gca,'fontsize', 40);

end

