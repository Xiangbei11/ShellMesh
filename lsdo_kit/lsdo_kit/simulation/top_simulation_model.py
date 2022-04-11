import csdl
from lsdo_kit.simulation.simulation_model import SimulationModel


class TopSimulationModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('simulation_dict')
        self.parameters.declare('design')

    def define(self):
        simulation_dict = self.parameters['simulation_dict']
        design = self.parameters['design']

        for simulation_name, simulation_obj in simulation_dict.items():
            # print('SIMULATION OBJ: ', simulation_obj)
            # print('SIMULATION NAME: ', simulation_name)
            sim_model = SimulationModel(simulation_name=simulation_name, simulation_obj=simulation_obj, design=design)
            self.add(sim_model, name=simulation_name, promotes=[])

                
                

