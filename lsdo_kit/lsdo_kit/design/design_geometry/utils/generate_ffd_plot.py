from vedo import Points, Plotter, LegendBox
import numpy as np 

def generate_ffd_plot(original_points, modified_points):

    pts_shape = original_points.shape 
    nxp = pts_shape[0]
    nyp = pts_shape[1]
    nzp = pts_shape[2]

    org_pts_reshape = np.reshape(original_points, (nxp*nyp*nzp, 3))
    mod_pts_reshape = np.reshape(modified_points, (nxp*nyp*nzp, 3))

    # print('Original Points: ', original_points[0,:,:,:])
    # print('Modified Points: ', modified_points[0,:,:,:])

    vp_init = Plotter()
    vps = []
    vps1 = Points(org_pts_reshape, r=8, c = 'blue')
    vps2 = Points(mod_pts_reshape, r=8, c='red')
    vps.append(vps1)
    vps.append(vps2)

    vp_init.show(vps, 'FFD Changes', axes=1, viewup="z", interactive = True)
