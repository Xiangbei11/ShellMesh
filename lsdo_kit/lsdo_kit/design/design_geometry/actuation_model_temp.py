
from locale import normalize
import csdl
from csdl_om import Simulator
import numpy as np

from lsdo_kit.design.design_geometry.core.actuation import Actuation
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets
from lsdo_kit.quat_rotate_vector_custom import QuatRotateVectorCustom

class ActuationModelTemp(csdl.Model):
    def initialize(self):
      self.parameters.declare('design')
      # self.parameters.declare('meshes')
      self.parameters.declare('actuation_dict')
      self.parameters.declare('nt')


    def define(self):
        design = self.parameters['design']
        nt = self.parameters['nt']
        actuation_dict = self.parameters['actuation_dict']
        geo = design.design_geometry

        pre_actuated_control_points = self.declare_variable('total_geometry_control_points', val=geo.total_cntrl_pts_vector)
        # print('SHAPE: ', pre_actuated_control_points.shape)
        pre_actuated_control_points_nt = csdl.expand(pre_actuated_control_points, (nt,) + pre_actuated_control_points.shape, 'ij->tij')
        self.register_output('actuated_control_points', pre_actuated_control_points_nt)