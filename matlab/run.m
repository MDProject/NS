ReadMatrix
opts.tol = 1e-10;
[EigenFuncs,EigenValues] = eigs(MatA,MatB,6,'sm',opts);
[Vx,Vz] = reshape_phi(EigenFuncs,L,H,Nx,Nz,ls1,ls2,2);
plot_velocity(Vx,Vz,L,H,10);