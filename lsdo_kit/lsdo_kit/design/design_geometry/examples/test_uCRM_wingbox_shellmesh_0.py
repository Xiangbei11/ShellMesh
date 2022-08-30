from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

path_name = 'CAD/'
file_name = 'uCRM-9_wingbox.stp' #_wing_with_tip #_wing_structure_PENGoLINS
geo = DesignGeometry(path_name + file_name, plot = True)

#top_projected_point1 = geo.extract_pointset(geo.input_bspline_entity_dict[''], np.array([0]), np.array([1]))


