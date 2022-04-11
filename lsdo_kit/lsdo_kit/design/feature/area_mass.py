import numpy as np

class AreaMass(object):
    def __init__(self, name, pointset):
        super().__init__(name, pointset)

    def add_mass(self, mass):
        if len(mass) == 1:
            self.mass = mass * np.ones(self.pointset.shape)

        elif mass.shape == self.pointset.shape[:1]:
            self.mass = mass
        
        else:
            raise Exception('Mass must either be a single number or equal the shape of the pointset')

            