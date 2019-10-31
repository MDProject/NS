% Vx,Vz,and Phi obey the same standard order
% boundary condition is set according to 2nd order accuracy

function [Vx,Vz] = reshape_phi(eigen_func,L,H,Nx,Nz,ls1,ls2,nidx)

DX = 2*L/(Nx-1);
DZ = 2*H/(Nz-1);

phi = eigen_func(:,nidx);
Vx = zeros(Nx,Nz);
Vz = zeros(Nx,Nz);
Phi = zeros(Nx,Nz);

for j = 2:1:Nz-1
    for i = 1:1:Nx-1
        Phi(i,j) = phi((j-2)*(Nx-1)+i);
    end
end
Phi(Nx,:) = Phi(1,:);

% using centered difference to calculate the velocity field
% calculate ghost points
R1 = (DZ-2*ls1)./(DZ+2*ls1);
R2 = (DZ-2*ls2)./(DZ+2*ls2);
Phi_down = R1.*Phi(:,2);
Phi_up = R2.*Phi(:,Nz-1);

% vx = \partial\phi/ \partial z, 4th order accuracy
for j = 3:1:Nz-2
    Vx(:,j) = (-Phi(:,j+2)+8*Phi(:,j+1)-8*Phi(:,j-1)+Phi(:,j-2))/12/DZ;
end
Vx(:,2) = (-Phi(:,4)+8*Phi(:,3)-8*Phi(:,1)+Phi_down)/12/DZ;
Vx(:,Nz-1) = (-Phi_up+8*Phi(:,Nz)-8*Phi(:,Nz-2)+Phi(:,Nz-3))/12/DZ;
% following 2nd order 
Vx(:,1) = (Phi(:,2) - Phi_down)/2/DZ;
Vx(:,Nz) = (Phi_up-Phi(:,Nz-1))/2/DZ;

% vz = -\partial\phi / \partial x
for i = 1:1:Nx
    if i == 1
        Vz(1,:) = (-Phi(3,:)+8*Phi(2,:)-8*Phi(Nx-1,:)+Phi(Nx-2,:))/12/DX;
    elseif i == 2
        Vz(2,:) = (-Phi(4,:)+8*Phi(3,:)-8*Phi(1,:)+Phi(Nx-1,:))/12/DX;
    elseif i == Nx-1
        Vz(Nx-1,:) = (-Phi(2,:)+8*Phi(Nx,:)-8*Phi(Nx-2,:)+Phi(Nx-3,:))/12/DX;
    elseif i == Nx
        Vz(Nx,:) = (-Phi(3,:)+8*Phi(2,:)-8*Phi(Nx-1,:)+Phi(Nx-2,:))/12/DX;
    else
        Vz(i,:) = (-Phi(i+2,:)+8*Phi(i+1,:)-8*Phi(i-1,:)+Phi(i-2,:))/12/DX;
    end
end
Vz = -Vz;

        


