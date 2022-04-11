import lsdo
import csdl
from csdl.utils.parameters import Parameters


class Base(object):

    def __init__(self, name, **kwargs):
        self.name = name
        self.parameters = Parameters()
        self.base_initialize()
        self.initialize()
        self.parameters.update(kwargs)


class Simulation(Base):

    def base_initialize(self, kwargs):
        self.parameters.declare('num_conditions', types=int)

    def dummy_method(self):
        temp=self.parameters['num_conditions']

class CruiseSimulation(Simulation):

    def initialize(self):
        self.parameters.declare('velocity', types=[str, float])
        self.parameters.declare('altitude')

    def dummy_vel_fnc(self):
        temp = self.parameters['velocity']


Cr = CruiseSimulation('asda')
