import numpy as np

class PointMass(object):
    def __init__(self, name):
        self.name = name
        self.point_mass_dict = {}


    def add_pointmass(self, point, mass):
        self.point_mass_list[point] = mass

    # def add_mass(self, mass):
    #     if len(mass) == 1:
    #         self.mass = mass

    # def add_point(self, point):
    #     self.point = point
            