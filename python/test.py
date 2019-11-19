from scipy import sparse
import numpy as np
from scipy.sparse.linalg import eigs
from DataManager import *

file_mat_bh = "D:\\WinSCP\\PDE\\Matrix\\biharmonic_py.txt"
file_mat_lap = "D:\\WinSCP\\PDE\\Matrix\\laplacian_py.txt"
file_setting = "D:\\WinSCP\\PDE\\Matrix\\setting.txt"
file_vec_slipLength = "D:\\WinSCP\\PDE\\Matrix\\ls.txt"

[L, H, Nx, Nz, MatA, MatB, ls1, ls2] = ReadData(file_mat_bh, file_mat_lap, file_setting, file_vec_slipLength)

[EigVal, EigVec] = eigs(MatA, 3, MatB, which='SM')

[Vx, Vz] = extract_velocity_field(EigVec, L, H, Nx, Nz, ls1, ls2, 2)

plot_velocity_field(Vx, Vz, L, H, 2)

print("asd")


