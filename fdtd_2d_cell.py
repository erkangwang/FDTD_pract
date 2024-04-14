import numpy as np


# 2D FDTD: have to consider derivatives d/dx d/dy, Yet d/dz=0
# Hz mode and Ez mode are de-coupled. need separate loops. Now I only finished the Ez mode
# array variables are stored in 2D arrays of size Nx*Ny
def fdtd_cell_2D(t_space_ind, src_ind):
    c0 = 3 * 10 ** 8  # m/s
    delta_t = 10 ** -8  # s
    Nx = 101
    Ny = 101
    Yee_Grid_X_Space_ind = np.linspace(0, Nx - 1, num=Nx, dtype=int)
    Yee_Grid_Y_Space_ind = np.linspace(0, Ny - 1, num=Ny, dtype=int)
    dx = 0.00005  # m
    dy = 0.00005  # m
    Yee_Grid_X_Space = Yee_Grid_X_Space_ind * dx
    Yee_Grid_Y_Space = Yee_Grid_Y_Space_ind * dy


    CE_x = np.zeros((Nx, Ny)) # curl of H field
    CE_y = np.zeros((Nx, Ny))  # curl of H field
    CH_z = np.zeros((Nx, Ny))  # curl of E field
    E_z = np.zeros((Nx, Ny))
    D_z = np.zeros((Nx, Ny))
    H_x = np.zeros((Nx, Ny))
    H_y = np.zeros((Nx, Ny))

    CH_x = np.zeros((Nx, Ny)) # curl of H field
    CH_y = np.zeros((Nx, Ny))  # curl of H field
    CE_z = np.zeros((Nx, Ny))  # curl of E field
    E_x = np.zeros((Nx, Ny))
    D_x = np.zeros((Nx, Ny))
    E_y = np.zeros((Nx, Ny))
    D_y = np.zeros((Nx, Ny))
    H_z = np.zeros((Nx, Ny))

    # update coefficients
    epsilon_zz = np.ones((Nx, Ny))  # array
    mu_xx = np.ones((Nx, Ny))  # array
    scale_factor = 5e-6
    m_Ey = (c0 * delta_t) / epsilon_zz * scale_factor
    m_Hx = - (c0 * delta_t) / mu_xx * scale_factor
    m_Hy = - (c0 * delta_t) / mu_xx * scale_factor
    #######
    # main FDTD loop
    tau = delta_t * 3
    g = np.exp(-(delta_t * (t_space_ind - 20) / tau) ** 2)
    E_z_ini=E_z
    E_z_all_time_list = E_z_ini[np.newaxis, :, :]

    for t in t_space_ind:

        #compute curls fo Ez mode
        for nx in Yee_Grid_X_Space_ind[0:Nx - 1]:
            for ny in Yee_Grid_Y_Space_ind[0:Ny - 2]:
                CE_x[nx][ny] = (E_z[nx][ny + 1] - E_z[nx][ny]) / dy
            CE_x[nx][Ny-1] = (0 - E_z[nx][Ny-1]) / dy
        
        for ny in Yee_Grid_Y_Space_ind[0:Ny - 1]:
            for nx in Yee_Grid_X_Space_ind[0:Nx - 2]:
                CE_y[nx][ny] = - (E_z[nx + 1][ny] - E_z[nx][ny]) / dx
            CE_y[Nx-1][ny] = - (0 - E_z[Nx-1][ny]) / dx
        
        CH_z[0][0] = (H_y[0][0] - 0) / dx - (H_x[0][0] - 0) / dy
        for nx in Yee_Grid_X_Space_ind[1:Nx - 1]:
            CH_z[nx][0] = (H_y[nx][0] - H_y[nx-1][0]) / dx - (H_x[nx][0] - 0) / dy

        for ny in Yee_Grid_Y_Space_ind[1:Ny - 1]:
            CH_z[0][ny] = (H_y[0][ny] - 0) / dx - (H_x[0][ny] - H_x[0][ny - 1]) / dy
            for nx in Yee_Grid_X_Space_ind[1:Nx - 1]:
                CH_z[nx][ny] = (H_y[nx][ny] - H_y[nx-1][ny]) / dx - (H_x[nx][ny] - H_x[nx][ny-1]) / dy

#        # compute curls fo Hz mode
#        for nx in Yee_Grid_X_Space_ind[0:Nx - 1]:
#            CH_x[nx][0] = (0 - H_z[nx][0]) / dy
#            for ny in Yee_Grid_Y_Space_ind[1:Ny - 1]:
#                CH_x[nx][ny] = (H_z[nx][ny + 1] - H_z[nx][ny]) / dy
#
#
#        for ny in Yee_Grid_Y_Space_ind[0:Ny - 1]:
#            CH_y[0][ny] = (0 - H_z[0][ny]) / dx
#            for nx in Yee_Grid_X_Space_ind[1:Nx - 1]:
#                CH_y[nx][ny] = - (H_z[nx + 1][ny] - H_z[nx][ny]) / dx
#
#
#        for ny in Yee_Grid_Y_Space_ind[0:Ny - 2]:
#            for nx in Yee_Grid_X_Space_ind[0:Nx - 2]:
#                CE_z[nx][ny] = (E_y[nx][ny] - E_y[nx - 1][ny]) / dx - (E_x[nx][ny] - E_x[nx][ny - 1]) / dy
#        CE_z[Nx - 1][ny] = (0 - E_y[Nx - 1][ny]) / dx - (0 - E_x[Nx - 1][ny]) / dy
#        for nx in Yee_Grid_X_Space_ind[0:Nx - 2]:
#            CH_z[nx][Ny - 1] = (E_y[nx][Ny - 1] - E_y[nx - 1][Ny - 1]) / dx - (E_x[nx][Ny - 1] - E_x[nx][Ny - 2]) / dy
#        CE_z[Nx - 1][Ny - 1] = (0 - E_y[Nx - 1][Ny - 1]) / dx - (0 - E_x[Nx - 1][Ny - 1]) / dy

        #updat Hx Hy and Dz field
        for ny in Yee_Grid_Y_Space_ind[0:Ny - 1]:
            for nx in Yee_Grid_X_Space_ind[0:Nx - 1]:
                H_x[nx][ny] = H_x[nx][ny] + m_Hx[nx][ny] * CE_x[nx][ny]
                H_y[nx][ny] = H_y[nx][ny] + m_Hy[nx][ny] * CE_y[nx][ny]
                #D_z[nx][ny] = D_z[nx][ny] + c0 * delta_t * CH_z[nx][ny]
                #E_z[nx][ny] = D_z[nx][ny] / epsilon_zz[nx][ny] * scale_factor
                E_z[nx][ny] = E_z[nx][ny] + m_Ey[nx][ny] * CH_z[nx][ny]

        # inject source
        E_z[src_ind][src_ind] = E_z[src_ind][src_ind] + g[t]
        norm_factor = np.amax(E_z)
        E_z=E_z/norm_factor
        #E_z_all_time_list = np.append(E_z_all_time_list, E_z, axis=None)

        E_z_all_time_list = np.concatenate((E_z_all_time_list, E_z[np.newaxis, :, :]), axis=0)
        #E_z_all_time_list.append(E_z)

        # E_y_all_time = np.array(E_y_all_time.T, E_y.T)
    ###########################################
    # out_array = np.append(Yee_Grid_Space, E_y, axis=0)
    # out_array = np.array([out_array[:101], out_array[101:]])
    out_array_list = E_z_all_time_list
    # out_array = np.array([Yee_Grid_Space.T, E_y.T])
    return out_array_list