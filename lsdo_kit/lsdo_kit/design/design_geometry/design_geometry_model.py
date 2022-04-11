import csdl
import numpy as np

from lsdo_kit.design.design_geometry.geometric_outputs_model import GeometricOutputsModel 
from lsdo_kit.design.design_geometry.ffd_model.ffd_model import FFDModel
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets


class DesignGeometryModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('design_geometry_obj')

    def define(self):
        design_geometry_obj = self.parameters['design_geometry_obj']

        # Add the Forward Evaluation Model so that we can use for Twist Optimization

        ffd_model = FFDModel(design_geometry_obj=design_geometry_obj)
        self.add(ffd_model, name='ffd_model', promotes=[])

        eval_pointsets = EvaluatePointsets(design_geometry_obj=design_geometry_obj)
        self.add(eval_pointsets, name='eval_pointsets', promotes=[])
        self.connect(f'ffd_model.application.total_geometry_control_points', f'eval_pointsets.control_points')

        # geometric_outputs_model = GeometricOutputsModel(geometry=geometry, geometric_outputs=geometric_outputs)
        # self.add(geometric_outputs_model, name='geo_outputs', promotes=[])

        # Add the GeometricOutputsModel and Geometric Constraints 

        if design_geometry_obj.geometric_outputs is not None:
            geometric_outputs_model = GeometricOutputsModel(design_geometry_obj=design_geometry_obj)
            self.add(geometric_outputs_model, name='geometric_outputs', promotes=[])
            self.connect('eval_pointsets.points', 'geometric_outputs.points')   


