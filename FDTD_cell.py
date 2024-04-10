import numpy as np

#1D FDTD: d/dx=d/dy=0 can simplify the E and H field updating in time and space
# array variables are stored in 1D arrays of length Nz
def fdtd_cell_1D(t_space_ind, src_ind):
    c0 = 3*10**8 # m/s
    delta_t = 10**-8 # s
    Nz = 101
    Yee_Grid_Space_ind = np.linspace(0, Nz-1, num=Nz, dtype=int)
    dz=0.01 # m
    Yee_Grid_Space=Yee_Grid_Space_ind*dz
    H_x = np.zeros(Nz)
    E_y = np.zeros(Nz)
    # update coefficients
    epsilon_yy = np.ones(Nz)   #array
    mu_xx = np.ones(Nz)   #array
    scale_factor = 0.00005
    m_Ey = (c0*delta_t)/epsilon_yy*scale_factor
    m_Hx = (c0*delta_t)/mu_xx*scale_factor
    #######
    #main FDTD loop
    E1 = E_y[0]
    tau = delta_t*3
    g = np.exp(-(delta_t*(t_space_ind-10)/tau)**2)
    E_y_all_time = Yee_Grid_Space
    for t in t_space_ind:

        H1 = H_x[0]
        for nz in Yee_Grid_Space_ind[0:Nz-2]:
            #1st H_x is at t+delta_t/2 time point, at nz th cell
            #2nd H_x is at t-delta_t/2 time point, at nz th cell
            #all E fields are at t time point, either nz+1 or nz th cell
            H_x[nz] = H_x[nz] + m_Hx[nz]*(E_y[nz+1] - E_y[nz])/dz
        H_x[Nz-1] = H_x[Nz-1] + m_Hx[Nz-1]*(E1 - E_y[Nz-1])/dz #out-of-boundary conditions

        E1 = E_y[0]
        E_y[0] = E_y[0] + m_Ey[0]*(H_x[0] - H1)/dz #out-of-boundary conditions
        for nz in Yee_Grid_Space_ind[1:Nz-1]: #Yee_Grid_Space_ind[1:]
            #1st E_y is at t+delta_t time point, at nz th cell
            #2nd E_y is at t time point, at nz th cell
            #all H fields are at t+delta_t/2 time point, either nz or nz-1 th cell
            E_y[nz] = E_y[nz] + m_Ey[nz]*(H_x[nz] - H_x[nz-1])/dz

        
        # inject source
        E_y[src_ind] = E_y[src_ind] + g[t]
        E_y_all_time= np.append(E_y_all_time, E_y, axis=0)
        #E_y_all_time = np.array(E_y_all_time.T, E_y.T)
    ###########################################
    #out_array = np.append(Yee_Grid_Space, E_y, axis=0)
    #out_array = np.array([out_array[:101], out_array[101:]])
    out_array = E_y_all_time
    #out_array = np.array([Yee_Grid_Space.T, E_y.T])
    return out_array

