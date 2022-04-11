from lsdo_kit.design.design_geometry.core.component import Component
import numpy as np

class Nacelle(Component):
    def __init__(self, name: str, stp_entity_names: list, ffd_obj=None, nxp=10, nyp=10, nzp=10):
        super().__init__(name, stp_entity_names, ffd_obj, nxp, nyp, nzp)

    def _compute_thrust_vector(self):
        _nacelle_length = self.x_max - self.x_min
        _quarter_nacelle_length = (1/4) * _nacelle_length
        _tip_indices = np.argwhere(self.embedded_entities_control_points == self.x_min)

        x_thrust_origin = self.x_min + _quarter_nacelle_length
        # y_thrust_origin = 
