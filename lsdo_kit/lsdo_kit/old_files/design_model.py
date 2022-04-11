from lsdo_kit.design.design_geometry.design_geometry_model import DesignGeometryModel 

class DesignModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('features')

    def define(self):
        features = self.parameters('features')

        design_geometry_model = DesignGeometryModel(ffd_blocks=ffd_blocks, geometry=geometry)
        self.add(ffd_model, name='ffd_model', promotes=[])
