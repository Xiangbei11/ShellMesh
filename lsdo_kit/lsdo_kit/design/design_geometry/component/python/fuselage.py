from lsdo_kit.components.component import Component

class Fuselage(Component):
    instance_index = 0

    def __init__(self, stp_entity_names: list, name=None, ffd_control_points=None):
        self.name = name
        if self.name is None:
            self.name = f'Fuselage_{Fuselage.instance_index}'
        super().__init__(self.name, stp_entity_names, ffd_control_points)
        Fuselage.instance_index += 1




if __name__ == "__main__":
    from lsdo_kit.design_geometry import DesignGeometry 
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' Using rect_wing as an example for a Component '''
    path_name = '../../../lsdo_geo/lsdo_geo/examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    fus1 = Fuselage(name='UniqueFuselage', stp_entity_names=['RectWing'])
    fus2 = Fuselage(stp_entity_names=['RectWing'])

    geo.add_component(fus1)
    geo.add_component(fus2)
    
    print(fus1.name)
    print(fus2.name)

    print(geo.components)