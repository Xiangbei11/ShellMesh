import lsdo
import csdl


lsdo_problem = lsdo.LSDOProblem()

lsdo_problem.design = design = lsdo.Design()

design.design_geometry = design_geometry = lsdo.DesignGeometry('geo.stp')

'''
Example 1:
'''

wing = lsdo.LiftingSurfaceComponent('main_wing', stp_entity_names=['left_wing', 'right_wing'])
wing.add_shape_parameter('aspect_ratio')
wing.add_shape_parameter('area')
design_geometry.add_component(wing)


'''
Example 2:
'''

wing = lsdo.LiftingSurfaceComponent('bot_wing', stp_entity_names=['left_wing', 'right_wing'])
wing.add_shape_parameter('twist', num=10, order=4)
wing.add_shape_parameter('shape', num_spanwise=6, num_chordwise=10, order_spanwise=4, order_chordwise=4)
design_geometry.add_component(wing)

wing = lsdo.LiftingSurfaceComponent('top_wing', stp_entity_names=['left_wing', 'right_wing'])
wing.add_twist_parameter('twist', num=10, order=4)
wing.add_shape_parameter('shape', num_spanwise=6, num_chordwise=10, order_spanwise=4, order_chordwise=4)
design_geometry.add_component(wing)
...

sim = csdl.Simulator()
# Insert code to define my preprocessing model
sim.model.create_input('twist')
sim.model.add_design_variable('twist')
sim.model.add(MyPreprocessingModel())
lsdo_model = lsdo_problem.setup()  # We need to implement this - it instantiates the LSDOProblemModel CSDL model 
sim.model.add(lsdo_model)


'''
Example 3
'''

wing = lsdo.LiftingSurfaceComponent('wing', stp_entity_names=['left_wing', 'right_wing'])
wing.add_span_parameter()
wing.add_chord_parameter()
design_geometry.add_component(wing)

tail = lsdo.LiftingSurfaceComponent('tail', stp_entity_names=['left_tail', 'right_tail'])
tail.add_span_parameter()
tail.add_chord_parameter()
design_geometry.add_component(tail)

fuselage = lsdo.Fuselage('fuselage', stp_entity_names=['fuselage'])
wing.add_length_parameter()
wing.add_width_parameter()
design_geometry.add_component(fuselage)

rotor = lsdo.Rotor('rotor', stp_entity_names=['rotor'])
rotor.add_radius_parameter()
design_geometry.add_component(rotor)

tail_fuselage_joint = GeometricConstraint('rigid', tail_pointset, fuse_pointset)
design_geometry.add_geometric_constraint(tail_fuselage_joint)

wing_span_model = csdl.Model()
rotor_radius = wing_span_model.declare_variable('rotor_radius')
wing_span = wing_span_model.declare_variable('wing_span')
constraint = 8 * rotor_radius + wing_span
wing_span_model.register_output('wing_span_constraint', constraint)
wing_span_constraint = GeometricConstraint('wing_span_constraint', wing_span_model)

# rotor_point_center = pointset_center
# rotor_point_tip = pointset_tip

# # geo_outputs = GeometricOuputs()
# # geo_outputs('rotor_radius', rotor_point_center, rotor_point_tip)


# # wing_point_left = pointset_left
# # wing_point_right = pointset_right
# # geo_outputs('wing_span', wing_point_left, wing_point_right)

rotor_wing_relation = DesignConstraint('GeoConstraint.rotor_wing_relation', )
