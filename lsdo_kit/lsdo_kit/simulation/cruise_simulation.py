from lsdo_kit.simulation.simulation import Simulation

class CruiseSimulation(Simulation):
    def __init__(self, name, velocity=None, density=None):
        
        # self.operating_conditions_dict = {
        #     'velocity' : velocity,
        #     'denisty'  : density,
        # }
        super().__init__(name)


    # def set_value(self, velocity : float = None , density : float = None):
    #     if velocity is not None:
    #         self.operating_conditions_dict['velocity'] = velocity
    #     if density is not None:
    #         self.operating_conditions_dict['density'] = density

    # def connect(self, velocity : str = None, density : str = None):
    #     if velocity is not None:
    #         self.operating_conditions_dict['velocity'] = velocity

    #     if density is not None:
    #         self.operating_conditions_dict['density'] = density        


if __name__ == "__main__":
    Cr = CruiseSimulation('cruise')
    Cr.set_value(velocity=50)
    Cr.connect(density='density_model.output')

    for attribute, value in Cr.__dict__.items():
        print('attribute: ', attribute)
        print('type attribute: ', type(attribute))
        print('value: ', value)



# class CruiseSimulation2(Simulation):
#     def initialize(self):
#         self.parameters.declare('velocity', types=[str, float], allow_none)
#         self.parameters.declare('density', types=[str, float], allow_none)