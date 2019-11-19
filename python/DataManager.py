from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt


def ReadData(file_mat_bh, file_mat_lap, file_setting, file_vec_slipLength):

    mat_bh = np.loadtxt(file_mat_bh)
    row = mat_bh[:, 0]
    col = mat_bh[:, 1]
    val = mat_bh[:, 2]
    row = row.astype(np.int32)
    col = col.astype(np.int32)
    MatA = sparse.csr_matrix((val, (row, col)))

    mat_lap = np.loadtxt(file_mat_lap)
    row = mat_lap[:, 0]
    col = mat_lap[:, 1]
    val = mat_lap[:, 2]
    row = row.astype(np.int32)
    col = col.astype(np.int32)
    MatB = sparse.csr_matrix((val, (row, col)))

    vec_setting = np.loadtxt(file_setting)
    L = vec_setting[0]
    H = vec_setting[1]
    Nx = vec_setting[2].astype(np.int32)
    Nz = vec_setting[3].astype(np.int32)

    # left side column ls1 for lower boundary and right side column ls2 for upper boundary
    vec_ls = np.loadtxt(file_vec_slipLength)
    ls1 = vec_ls[:, 0]
    ls2 = vec_ls[:, 1]
    return [L, H, Nx, Nz, MatA, MatB, ls1, ls2]

def extract_velocity_field(EigenVec, L, H, Nx, Nz, ls1, ls2, nidx):

    DX = 2 * L / (Nx - 1)
    DZ = 2 * H / (Nz - 1)
    phi = EigenVec[:, nidx]
    Vx = np.zeros((Nx + 1, Nz + 1))
    Vz = np.zeros((Nx + 1, Nz + 1))
    Phi = np.zeros((Nx + 1, Nz + 1))

    for j in range(2, Nz):
        for i in range(1, Nx):
            Phi[i, j] = phi[(j-2)*(Nx-1)+i-1]
    Phi[Nx, :] = Phi[1, :]  # periodic boundary

    R1 = (DZ-2*ls1)/(DZ+2*ls1)
    R2 = (DZ-2*ls2)/(DZ+2*ls2)
    R1 = np.append(0, R1)
    R2 = np.append(0, R2)
    Phi_down = np.multiply(R1, Phi[:, 2])
    Phi_up = np.multiply(R2, Phi[:, Nz-1])

    # 4th order accuracy
    for j in range(3, Nz - 1):
        Vx[:, j] = (-Phi[:, j + 2] + 8*Phi[:, j + 1] - 8*Phi[:, j - 1] + Phi[:, j - 2])/12/DZ
    Vx[:, 2] = (-Phi[:, 4] + 8*Phi[:, 3] - 8*Phi[:, 1] + Phi_down)/12/DZ
    Vx[:, Nz - 1] = (-Phi_up + 8*Phi[:, Nz] - 8*Phi[:, Nz - 2] + Phi[:, Nz - 3])/12/DZ
    # 2nd order accuracy
    Vx[:, 1] = (Phi[:, 2] - Phi_down)/2/DZ
    Vx[:, Nz] = (Phi_up - Phi[:, Nz - 1])/2/DZ

    for i in range(Nx + 1):
        if i == 1:
            Vz[1, :] = (-Phi[3, :] + 8*Phi[2, :] - 8*Phi[Nx - 1, :] + Phi[Nx - 2, :])/12/DX
        elif i == 2:
            Vz[2, :] = (-Phi[4, :] + 8*Phi[3, :] - 8*Phi[1, :] + Phi[Nx - 1, :])/12/DX
        elif i == Nx - 1:
            Vz[Nx - 1, :] = (-Phi[2, :] + 8 * Phi[Nx, :] - 8 * Phi[Nx - 2, :] + Phi[Nx - 3, :]) / 12 / DX
        elif i == Nx:
            Vz[Nx, :] = Vz[2, :] = (-Phi[3, :] + 8*Phi[2, :] - 8*Phi[Nx - 1, :] + Phi[Nx - 2, :])/12/DX
        else:
            Vz[i, :] = Vz[2, :] = (-Phi[i + 2, :] + 8*Phi[i + 1, :] - 8*Phi[i - 1, :] + Phi[i - 2, :])/12/DX
    Vz = -Vz
    return [Vx, Vz]

def plot_velocity_field(Vx, Vz, L, H, scaling):

    Nx = Vx.shape[0]
    Nz = Vx.shape[1]
    nx = Nx/scaling
    nz = Nz/scaling
    step_x = Nx/nx
    step_z = Nz/nz

    x = np.linspace(-L, L, Nx)
    z = np.linspace(-H, H, Nz)
    [X, Z] = np.meshgrid(x, z)
    fig, ax = plt.subplots()
    q = ax.quiver(X, Z, Vx.transpose(), Vz.transpose())
    plt.show()










