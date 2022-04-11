import csdl
import numpy as np

from lsdo_kit.design.feature.area_mass import AreaMass

class ForcesMomentsModel(csdl.Model):
    def initialize(self):
        self.declare.parameters['feature_dict']

    def define(self):
        feature_dict = self.parameters['feature_dict']
        
        moments_list = []
        mass_list = []

        for feature_name, feature_obj in feature_dict.items():
            if type(feature_obj) is AreaMass:
                pass

            else:
                mass = self.create_input(f'{feature_name}_mass', val=feature_obj.mass)
                cg   = self.create_input(f'{feature_name}_cg',   val=feature_obj.pointset)
                moment = mass * cg 
                mass_list.append(mass)
                moments_list.append(moment)

        total_mass = 0
        for mass_obj in mass_list:
            total_mass = total_mass + mass_obj
        
        total_moment = 0
        for moment_obj in moments_list:
            total_moment = total_moment + moment_obj

        total_cg = total_moment / total_cg

        self.register_output('total_mass', total_mass)
        self.register_output('total_moment', total_moment)
        self.register_output('total_cg', total_cg)

        