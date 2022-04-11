from lsdo_kit.design.design_geometry.core.component import Component
import numpy as np 

class LiftingSurface(Component):
    instance_index = 0
    def __init__(self, stp_entity_names: list, name=None, ffd_control_points=None):
        self.name = name
        if self.name is None:
            self.name = f'LiftingSurface_{LiftingSurface.instance_index}'
        super().__init__(self.name, stp_entity_names, ffd_control_points)
        LiftingSurface.instance_index += 1

    # TODO: Modify Section Properties so that it can take a string value and issue a connection at FFD model level
    def add_constant_twist(self, twist):
        if type(twist) is str:
            self.ffd.add_shape_parameter(property_name = 'rot_x', parameter_name='constant', order=1, num_cp=2, dv=True, val=twist)

    def add_linear_twist(self, design_var=False, val=np.ones(2)):
        self.ffd.add_shape_parameter(property_name = 'rot_x', parameter_name='linear', order=2, num_cp=3, dv=design_var, val=val)

    def add_nth_order_twist(self, order, num_cp, design_var, val):
        self.ffd.add_shape_parameter(property_name ='rot_x', parameter_name=f'{order}_order', order=order, num_cp=num_cp, dv=design_var, val=val)

if __name__ == "__main__":
    from lsdo_kit.design.design_geometry import DesignGeometry 
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' Using rect_wing as an example for a Component '''
    path_name = '../../../lsdo_geo/lsdo_geo/examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    LiftingSurface1 = LiftingSurface(name='UniqueWing', stp_entity_names=['RectWing'])
    LiftingSurface2 = LiftingSurface(stp_entity_names=['RectWing'])

    geo.add_component(LiftingSurface1)
    geo.add_component(LiftingSurface2)
    
    print(LiftingSurface1.name)
    print(LiftingSurface2.name)

    print(geo.components)