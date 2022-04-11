import numpy as np

def calculate_rotation_mat(v1, v2):
    '''
    This does NOT work when the vectors point in exactly oppposite directions from each other!!
    Add check, and work around for this! 

    Link to the calculation:
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    '''
    
    v = np.cross(v1,v2)
    s = np.linalg.norm(v)
    c = np.dot(v1, v2)

    v_skew = [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]

    rot_mat = np.eye(3) + v_skew + np.linalg.matrix_power(v_skew,2)*(1/(1+c))

    return rot_mat
