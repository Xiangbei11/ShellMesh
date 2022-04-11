import csdl
from lsdo_kit.design.design_geometry.design_geometry_model import DesignGeometryModel 

class DesignModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('design_obj')

    def define(self):
        design_obj = self.parameters['design_obj']

        design_geometry_obj = design_obj.design_geometry
       
        # Needs to unpack and add DesignGeometryModel and Features
        # feature_dict = design_obj.feature_dict

        design_geometry_model = DesignGeometryModel(design_geometry_obj=design_geometry_obj)
        self.add(design_geometry_model, name='design_geometry_model', promotes=[])
