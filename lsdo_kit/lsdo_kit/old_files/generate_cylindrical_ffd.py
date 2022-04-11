import numpy as np 

'''
This function assumes that the control point the user is feeding in are 
8 values to define the corners of a rectangular prism that would be able 
to contain the entity.

This script then takes those points, and creates a cylinder. 

x-axis: The height direction
y-axis: The direction from which the diameter of the circle will be determined
z-axis: We are assuming a cylinder with circles at the end, therefore same as y-axis
'''

def calc_center_of_rectangle(control_points):
    y_center_1 = ((control_points[0,0,0,:]) + (control_points[0,1,0,:])) / 2 
    z_center_1 = ((control_points[0,0,0,:]) + (control_points[0,0,1,:])) / 2

    y_center_2 = ((control_points[0,0,1,:]) + (control_points[0,1,1,:])) / 2 
    z_center_2 = ((control_points[0,1,0,:]) + (control_points[0,1,1,:])) / 2

    y_center = (y_center_1 + y_center_2) / 2
    z_center = (z_center_1 + z_center_2) / 2

    return y_center, z_center

def max_y_radius(control_points):
    temp1 = np.sqrt((control_points[0,0,0,:] - control_points[0,1,0,:]) ** 2)
    temp2 = np.sqrt((control_points[1,0,1,:] - control_points[1,1,1,:]) ** 2)
    max_rad_y_1 = max(temp1[1], temp2[1])/2

    temp1 = np.sqrt((control_points[1,0,0,:] - control_points[1,1,0,:]) ** 2)
    temp2 = np.sqrt((control_points[1,0,1,:] - control_points[1,1,1,:]) ** 2)
    max_rad_y_2 = max(temp1[1], temp2[1])/2

    return max_rad_y_1, max_rad_y_2

def create_ffd(control_points, nxp, nyp, ncp):

    # nxp refers to the amount of points that make up the height of the cylinder
    # nyp refers to the amount of points along the radius of the end of the cylinder
    # ncp refers to the amount of points in one of the circles located at the ends of the cylinder

    t = np.linspace(0, 2*np.pi, ncp, endpoint=False)

    max_rad_circle_1, max_rad_circle_2 = max_y_radius(control_points)

    y_center, z_center = calc_center_of_rectangle(control_points)

    x_largest_radius_circle_1 = max_rad_circle_1 * np.cos(t)
    z_largest_radius_circle_1 = max_rad_circle_1 * np.sin(t)

    x_largest_radius_circle_2 = max_rad_circle_2 * np.cos(t)
    z_largest_radius_circle_2 = max_rad_circle_2 * np.sin(t)

    # Create the inner points to the circle on each end of the cylinder
    # x_circle_1_inner_points = np.linspace()
    # z_circle_1_inner_points = 
    # print(max_rad_y_1)
    # print(max_rad_y_2)

    return x_largest_radius_circle_1, z_largest_radius_circle_1, y_center, z_center




if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox


    nxp = 10
    nyp = 5
    ncp = 10

    point000 = np.array([170. ,0. ,100.])
    point001 = np.array([170., 0., 170.])

    point010 = np.array([170., 100., 100.])
    point011 = np.array([170., 100., 170.])

    point100 = np.array([240. ,0. ,100.])
    point101 = np.array([240. ,0. ,170.])

    point110 = np.array([240. ,100. ,100.])
    point111 = np.array([240. ,100. ,170.])

    control_points = np.zeros((2, 2, 2, 3))

    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001
    
    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011

    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
  
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    x_circle_1, z_circle_1, y_center, z_center = create_ffd(control_points, nxp, nyp, ncp)

    print(y_center)

    print(x_circle_1)

    plt.scatter(x_circle_1, z_circle_1)
    plt.scatter(y_center[1], y_center[2])
    plt.show()


    '''
    Have the user provide a two centers, two radii, and a length for the cylinder. 

    If they provide one center, then they can provide one radius 

    If they provide two centers, they can still only provide one radius, or they can provide 2 radii
    '''