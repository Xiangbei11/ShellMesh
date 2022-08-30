from options_dictionary import OptionsDictionary

import numpy as np
import scipy.sparse as sps

class Member(object):
    def __init__(self, **kwargs):
        ''' Add member to the aircraft model and update the list of members. mem_name dictionary
        ----------
        'name' : str 
            Name of the member.
        'points_to_be_projected': np.array
            with shape (2, num_points, 3) 
        'projection_direction': np.array
            with shape (2,3)/(2, num_points, 3)
        '''
        self.options = OptionsDictionary()
        self.options.declare('name',types=str)
        self.options.declare('id', types=int)
        self.options.declare('node_indices', types=list, default=[])
        self.options.declare('tri_connectivity', types=np.ndarray)
        self.options.declare('constrained_node_indices', types=list, default=[])
        self.options.declare('u_v_vec', types=np.ndarray)
        self.options.declare('mapping', types=sps.csc.csc_matrix)
        self.options.declare('constrained_edges', types=np.ndarray)
        self.options.declare('constrained_boundary_node_indices', types=list, default=[])
        self.options.declare('coordinates', types=np.ndarray)
        self.options.declare('opt_coordinates', types=np.ndarray)
        self.options.declare('opt_connectivity', types=np.ndarray)
        self.options.declare('opt_mapping', types=sps.csc.csc_matrix)
        self.options.declare('ctrl_pointset_list', types=list, default=[])
        self.options.update(kwargs)        
        # 
        # self.declare('curves',types=np.ndarray,default=np.array([0,0,0,0]))
        # 
        # self.declare('dirts',types=np.ndarray,default=np.array([]))
        # self.declare('weight',types=np.ndarray,default=np.array([]))
        # self.declare('r',types=float,default=0.)
        # self.declare('cc',types=np.ndarray,default=np.array([]))
        # self.declare('type',types=int,default=0)
        # self.declare('pt_names',default=[])
        # self.declare('pt_dict',default={'starting_upper_point': np.empty([3],dtype=np.float32),
        #                                                 'ending_upper_point': np.empty([3],dtype=np.float32),
        #                                                 'starting_lower_point': np.empty([3],dtype=np.float32),
        #                                                 'ending_lower_point': np.empty([3],dtype=np.float32)})
        
        # self.declare('DT_tris',types=np.ndarray)
        # self.declare('DT_coords',types=np.ndarray)
        # self.declare('DT_fixed',types=np.ndarray)
        # self.declare('projected',types=list,default=[])
        # self.declare('quads',types=np.ndarray)
        # self.declare('coords',types=np.ndarray)
        # self.update(kwargs)

        # self.submembers = []


if __name__ == '__main__':
    test = Member()
    print(test)