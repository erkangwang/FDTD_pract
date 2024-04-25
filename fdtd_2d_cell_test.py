import numpy as np


# 2D FDTD: have to consider derivatives d/dx d/dy, Yet d/dz=0
# Hz mode and Ez mode are de-coupled. need separate loops. Now I only finished the Ez mode
# array variables are stored in 2D arrays of size Nx*Ny
def fdtd_cell_2D_test(t_space_ind, src_ind):

    e0 = 8.854 * 10 ** -12
    mu0 = 4 * np.pi * 10 ** -7
    c0 = np.sqrt(1/(e0*mu0))  # m/s
    c=1
    lamda = 500 * 10 ** -9  # m
    f = lamda / c0
    Nx = 101
    Ny = 101
    Yee_Grid_X_Space_ind = np.linspace(0, Nx - 1, num=Nx, dtype=int)
    Yee_Grid_Y_Space_ind = np.linspace(0, Ny - 1, num=Ny, dtype=int)
    dx = 1
    dy = 1
    dt = min(dx, dy) / np.sqrt(2) / c
    Yee_Grid_X_Space, Yee_Grid_Y_Space = np.meshgrid(Yee_Grid_X_Space_ind*dx, Yee_Grid_Y_Space_ind*dy)

    E_z = np.zeros((Nx, Ny))
    H_x = np.zeros((Nx, Ny-1))
    H_y = np.zeros((Nx-1, Ny))

    #######
    # main FDTD loop
    tau = dt * 20
    g = np.exp(-(dt*(t_space_ind-50)/tau)**2)
    E_z_ini=E_z
    E_z_all_time_list = E_z_ini[np.newaxis, :, :]

    for t in t_space_ind:

        # Update magnetic fields at time step n+1/2
        H_x = H_x - dt / dy * (E_z[:, 1:] - E_z[:, :-1])
        H_y = H_y + dt / dx * (E_z[1:, :] - E_z[:-1, :])

        # Update electric field at time step n+1
        diff_H_x = dt / dy * (H_x[1:-1, 1:] - H_x[1:-1, :-1])
        diff_H_y = dt / dx * (H_y[1:, 1:-1] - H_y[:-1, 1:-1])
        E_z[1:-1, 1:-1] = E_z[1:-1, 1:-1] + (diff_H_y - diff_H_x)

        E_z[src_ind][src_ind] = E_z[src_ind][src_ind] + g[t]
        E_z_all_time_list = np.concatenate((E_z_all_time_list, E_z[np.newaxis, :, :]), axis=0)


        # E_y_all_time = np.array(E_y_all_time.T, E_y.T)
    ###########################################
    # out_array = np.append(Yee_Grid_Space, E_y, axis=0)
    # out_array = np.array([out_array[:101], out_array[101:]])
    out_array_list = E_z_all_time_list
    # out_array = np.array([Yee_Grid_Space.T, E_y.T])
    return out_array_list, g