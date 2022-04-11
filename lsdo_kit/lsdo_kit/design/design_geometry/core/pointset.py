import numpy as np
from typing import List


class PointSet:
    def __init__(self, pointset_id, output_starting_ind, shape, parent_pointers, relative_map, absolute_map,
                offset, offset_absolute_map, physical_coordinates, permutation_matrix, name=None) -> None:
        '''
        pointset_id: int
        shape : np.ndarray # shape=(n,3) or (nu, nv 3), or (nu, nv, nw, 3)
        parent_pointers : List  # list of IDs to parent pointsets or BSplineEntities
        output_starting_ind : int, where this pointset subset starts in the entire pointset vector
        relative_map : np.ndarray, map from parent pointset (geometry control points) to a subset of all pointsets
        absolute_map : np.ndarray, map from parent pointset (geomeetry control points) to all pointsets
        offset : np.ndarray
        offset_absolute_map : np.ndarray
        physical_coordinates : np.ndarray
        '''

        self.pointset_id = pointset_id
        if name is None:
            self.name = f'pointset_{self.pointset_id}'
        else:
            self.name = name
            
        self.output_starting_ind = output_starting_ind
        self.shape = shape # shape=(n,3) or (nu, nv 3), or (nu, nv, nw, 3)
        self.parent_pointers = parent_pointers  # list of IDs to parent pointsets or BSplineEntities
        self.relative_map = relative_map
        self.absolute_map = absolute_map
        self.offset = offset
        self.offset_absolute_map = offset_absolute_map
        self.physical_coordinates = physical_coordinates
        self.permutation_matrix = permutation_matrix

    def translate(self, x, y, z):
        # TODO Check to see if this makes sense
        self.offset += np.array(x, y, z)

    def reshape(self, shape):
        # TODO Check to see if this makes sense
        self.shape = shape

    def reorder(self, einsum_string):
        # TODO implement this (do things like transpose, but for more dimensions)
        print("Sorry! This has not been implemented yet!")
        pass

    def extract(self, location):
        # TODO implement this. Get point/pointset given parametric coordinates
        print("Sorry! This has not been implemented yet!")
        pass

    def __add__(self, array):
        # TODO Check this!
        self.translate(array[0], array[1], array[2])

    def __mult__(self, scaling_value):
        # TODO Check this!
        pass
    
    def __sub__(self, array):
        self.translate(-array[0], -array[1], -array[2])


    def __equals__(self, other):
        # TODO implement this (and recall what we are trying to do with this)
        print("Sorry! This has not been implemented yet!")
        pass

    
class ProjectedPointSet(PointSet):
    pass

class DerivedPointSet(PointSet):
    pass

class EvaluatedPointSets(PointSet):     
    pass