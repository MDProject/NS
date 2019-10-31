function [] = plot_velocity(Vx,Vz,L,H,scaling)

[Nx,Nz] = size(Vx);
nx = floor(Nx/scaling);
nz = floor(Nz/scaling);
step_x = floor(Nx/nx);
step_z = floor(Nz/nz);

x = linspace(-L,L,Nx);
z = linspace(-H,H,Nz);
[X,Z] = meshgrid(x,z);
Vx = Vx';
Vz = Vz'; % transpose is important!

% data points are too dense, scaling the arrows
X = X(1:step_z:Nz,1:step_x:Nx);
Z = Z(1:step_z:Nz,1:step_x:Nx);
Ux = Vx(1:step_z:Nz,1:step_x:Nx);
Uz = Vz(1:step_z:Nz,1:step_x:Nx);
quiver(X,Z,Ux,Uz); 

end