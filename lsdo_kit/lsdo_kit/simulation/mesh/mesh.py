import numpy as np

class Mesh:
    '''
    Class for representing meshes to be fed into analysis components.
    '''

    def __init__(self, name, pointset_list=[]) -> None:
        self.name = name
        self.pointset_list = []
        self.num_mesh_points = 0

        for pointset in pointset_list:
            self.add_pointset(pointset)

    def add_pointset(self, pointset, name=None):
        '''
        Create/add on to the single instance of evaluated_pointsets
        '''
        if name is not None:
            pointset.name = name

        elif pointset.name is None:
            name = f'pointset_{pointset.pointset_id}'


        self.pointset_list.append(pointset)
        self.num_mesh_points += np.cumprod(pointset.shape)[-2]
