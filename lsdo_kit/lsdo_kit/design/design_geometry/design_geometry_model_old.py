from lsdo_kit.design.design_geometry.ffd_model.inner_optimization_model import InnerOptimizationModel 
from lsdo_kit.design.design_geometry.geometric_outputs_model import GeometricOutputsModel 
import csdl
import numpy as np

class DesignGeometryModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('geometry')
    
    def define(self):
        geometry = self.parameters('geometry')

        # Since geometry contains all the components, and the components contain the respective FFD blocks, couldnt we just feed in geometry and extract the ffd block info
        # I added a dictionary in geometry that appends every time a component is added with a key which is the component name, and the value is the respective ffd block

        inner_optimization_model = InnerOptimizationModel(ffd_blocks=ffd_blocks, geometry=geometry)
        self.add(inner_optimization_model, name='inner_optimization_model', promotes=[])

        geometric_outputs_model = GeometricOutputsModel()