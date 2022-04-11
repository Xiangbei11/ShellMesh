import numpy as np 


def _create_n_array(n, dim):
    if isinstance(n, int):
        n = np.repeat(n, dim)
        num_points = np.sum(n) - (dim - 1)
        return n, int(num_points)

    elif (isinstance(n, np.ndarray) and n.size == dim):
        num_points = np.sum(n) - (dim - 1)
        return n, int(num_points)


def create_ffd(control_points, nxp, nyp, nzp):
    zdim = control_points.shape[2] - 1
    ydim = control_points.shape[1] 
    xdim = control_points.shape[0]

    nzp, num_points_z = _create_n_array(nzp, zdim)

    bot_top_lines = np.zeros((xdim, ydim, num_points_z, 3))
    for i in range(xdim):
        for j in range(ydim):
            for k in range(zdim):
                if k == 0: 
                    bot_top_lines[i, j, 0:nzp[k], :] = np.linspace(control_points[i,j,k,:], control_points[i,j,k+1,:], nzp[k])
               
                elif k >= 1: 
                    temp = np.zeros((xdim, ydim, nzp[k], 3))

                    # print('temp: ', temp)

                    temp[i, j, :, :] = np.linspace(control_points[i,j,k,:], control_points[i,j,k+1,:], nzp[k])
                    
                    # print('temp val: ', temp)

                    # print('k: ', k)
                    # print('temp: ', temp[i, j, :, :])
                    # print('temp shape: ', temp[i,j,:,:].shape)
                    # print('bot_top_lines: ',  bot_top_lines[i, j, (k*nzp-k):((k+1)*nzp - k), :])
                    # print('bot_top_lines shape: ', bot_top_lines[i, j, (k*nzp-k):((k+1)*nzp - k), :].shape)
                    # print('bot_top: ', bot_top_lines[i, j, :, :])
                    # print('bot_top shape: ', bot_top_lines[i, j, :, :].shape)
                    
                    bot_top_lines[i, j, (np.sum(nzp[:k]) - k) : (np.sum( nzp[:(k+1)] ) - k), :] = temp[i, j, :, :]

    # return bot_top_lines
                

    ydim = control_points.shape[1] - 1
    xdim = control_points.shape[0]

    nyp, num_points_y = _create_n_array(nyp, ydim)

    # print('ydim: ', ydim)
    # print('bot_top shape: ', bot_top_lines.shape)


    section_mesh = np.zeros((xdim, int(num_points_y), bot_top_lines.shape[2], 3))
    # print('section_mesh: ', section_mesh)
    for i in range(xdim):
        for j in range(ydim):
            if j == 0: 
                section_mesh[i, 0:nyp[j], :, :] = np.linspace(bot_top_lines[i,j,:,:], bot_top_lines[i,j+1,:,:], nyp[j]) 
            elif j >= 1: 
                # print('temp: ', temp)
                # print('temp shape: ', temp.shape)
                
                # print('temp: ', np.linspace(bot_top_lines[i,j,:,:], bot_top_lines[i,j+1,:,:], nyp))
                # print('temp.shape: ', np.linspace(bot_top_lines[i,j,:,:], bot_top_lines[i,j+1,:,:], nyp).shape)
                temp = np.zeros((xdim, nyp[j], bot_top_lines.shape[2], 3))
                # print('temp: ', temp)
                # print('temp shape: ', temp.shape)

                temp[i, :, :, :] = np.linspace(bot_top_lines[i,j,:,:], bot_top_lines[i,j+1,:,:], nyp[j])
                # print('temp val: ', temp)
                # print('temp val shape: ', temp.shape)

                # print('section mesh: ', section_mesh)
                # print('section mesh shape: ', section_mesh.shape)

                # print('temp: ', temp[i, :, :, :])
                # print('temp shape: ', temp[i, :, :, :].shape)
                
                
                section_mesh[i, (np.sum(nyp[:j]) - j):(np.sum(nyp[:(j+1)] ) - j), :, :] = temp[i, :, :, :]
    # return bot_top_lines, section_mesh

    # return bot_top_lines, section_mesh


    xdim = control_points.shape[0] - 1
    nxp, num_points_x = _create_n_array(nxp, xdim)


    ffd_control_points = np.zeros((int(num_points_x), section_mesh.shape[1], section_mesh.shape[2], 3))

    for i in range(xdim):
        if i == 0: 
            ffd_control_points[0:nxp[i],:,:,:] = np.linspace(section_mesh[i,:,:,:], section_mesh[i+1,:,:,:], nxp[i]) 

        elif i >= 1: 
            temp = np.zeros((nxp[i], section_mesh.shape[1], section_mesh.shape[2], 3))
            print('temp: ', temp)
            print('temp shape: ', temp.shape)

            temp = np.linspace(section_mesh[i,:,:,:], section_mesh[i+1,:,:,:], nxp[i])
            print('temp val: ', temp)
            print('temp val shape: ', temp.shape)

            ffd_control_points[(np.sum(nxp[:i]) - i):(np.sum(nxp[:(i+1)] ) - i), :, :, :] = temp

    return ffd_control_points

    
    # xdim = (control_points.shape[0] - 1) * nxp
    # ydim = (control_points.shape[1] - 1) * nyp
    # zdim = (control_points.shape[2] - 1) * nzp
    # ffd_control_points = np.zeros((xdim, ydim, zdim, 3))
    # shape = (nxp, nyp, nzp)

    # if control_points.shape[-1] == 3:
    #     for i in range(control_points.shape[0] - 1):
    #         for j in range(control_points.shape[1] - 1):
    #             for k in range(control_points.shape[2] - 1):


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    nxp = 5
    nyp = 5
    nzp = 5

    point000 = np.array([170. ,0. ,100.])
    point001 = np.array([170., 0., 170.])

    point002 = np.array([170., 0., 240.])
    point003 = np.array([170., 0., 300.])

    point010 = np.array([170., 100., 100.])
    point011 = np.array([170., 100., 170.])
    
    point012 = np.array([130., 100., 240.])
    point013 = np.array([130., 100., 300.])

    point020 = np.array([130., 200., 100.])
    point021 = np.array([130., 200., 170.])
    point022 = np.array([130., 200., 240.])
    point023 = np.array([130., 200., 300.])

    point100 = np.array([240. ,0. ,100.])
    point101 = np.array([240. ,0. ,170.])
    point102 = np.array([240. ,0. ,240.])
    point103 = np.array([240. ,0. ,300.])

    point200 = np.array([50. ,0. ,100.])
    point201 = np.array([50. ,0. ,170.])
    point202 = np.array([50. ,0. ,240.])
    point203 = np.array([50. ,0. ,300.])

    point110 = np.array([240. ,100. ,100.])
    point111 = np.array([240. ,100. ,170.])
    point112 = np.array([200. ,100. ,240.])
    point113 = np.array([200. ,100. ,300.])
    
    point210 = np.array([300. ,100. ,100.])
    point211 = np.array([300. ,100. ,170.])
    point212 = np.array([200. ,100. ,240.])
    point213 = np.array([200. ,100. ,300.])
    
    point220 = np.array([500. ,200. ,100.])
    point221 = np.array([500. ,200. ,170.])
    point222 = np.array([500. ,200. ,240.])
    point223 = np.array([200. ,200. ,300.])

    control_points = np.zeros((2, 2, 2, 3))

    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001
    # control_points[0,0,2,:] = point002
    
    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011
    # control_points[0,1,2,:] = point012

    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
    # control_points[1,0,2,:] = point102    
  
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111
    # control_points[1,1,2,:] = point112

    # control_points[1,1,0,:] = point110
    # control_points[1,1,1,:] = point111
    # control_points[1,1,2,:] = point112


    # control_points[2,0,0,:] = point200
    # control_points[2,0,1,:] = point201
    # control_points[2,0,2,:] = point202

    # control_points[2,1,0,:] = point210
    # control_points[2,1,1,:] = point211
    # control_points[2,1,2,:] = point212

    # control_points[2,2,0,:] = point220
    # control_points[2,2,1,:] = point221
    # control_points[2,2,2,:] = point222

    # control_points[0,1,3,:] = point013

    # control_points[0,2,0,:] = point020
    # control_points[0,2,1,:] = point021
    # control_points[0,2,2,:] = point022
    # control_points[0,2,3,:] = point023



    # control_points[1,0,3,:] = point103
    

    # control_points[1,1,3,:] = point113


    # control_points[1,2,0,:] = point120
    # control_points[1,2,1,:] = point121
    # control_points[1,2,2,:] = point122
    # control_points[1,2,3,:] = point123


    ffd_control_points = create_ffd(control_points, nxp, nyp, nzp)

    ''' 
    PLOT
    '''
    
    ffd_control_points1 = np.reshape(ffd_control_points, (ffd_control_points.shape[0] * ffd_control_points.shape[1] * ffd_control_points.shape[2], 3))
    # section_mesh1 = np.reshape(section_mesh, ((section_mesh.shape[0] * section_mesh.shape[1] * section_mesh.shape[2], 3)))

    vp_init = Plotter()
    vps1 = Points(ffd_control_points1, r=8, c = 'red')
    # # # vps2 = Points(section_mesh1, r=8, c = 'red')


    vp_init.show(vps1, 'Bspline Volume', axes=1, viewup="z", interactive = True)
    # vp_init.show(vps2, 'Section Mesh', axes=1, viewup="z", interactive = True)

    # print(control_points)
    # print(control_points.shape)
