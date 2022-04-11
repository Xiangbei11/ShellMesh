# from lsdo import LSDOProblem, Design, Simulation
import lsdo
import csdl

from lsdo_kit.design.design_geometry.core.pointset import PointSet


lsdo_problem = lsdo.LSDOProblem()

lsdo_problem.design = design = lsdo.Design()

design.design_geometry = design_geometry = lsdo.DesignGeometry('geo.stp')

'''
Design and DesignGeometry Discussion
-------------------------------------
DESIGN PYTHON: 
. Add both components and features to the deisgn class. The Design class then feeds the component information down to DeisgnGeometry
    . Adding components/features to design just appends a list or dictionary. 
        . Make dictionary key = to the component or feature name 

. Geometric outputs also get added to the design class. 

. Treat this class as an information parser

DESIGN CSDL:
. Calls other CSDL models 
. Feeds the parameters down the line as necessary 

DesignGeomtry PYTHON:
. Is just the geometry class itself

DesignGeometry CSDL:
. Feeds parameters into the respective model
. Runs InnerOptModel, GeoOutputsModel, GeoConstraintsModel

'''

'''
Feature Discussion
-------------------
PYTHON:

. User feeds in a pointset, along that may contain addition physical information to describe that pointset. 

. Feautures are a heirarchy. The feature object needs an "add_children" method to provide ability to create this heirarchy.  


'''


wing = lsdo.LiftingSurfaceComponent('main_wing', stp_entity_names=['left_wing', 'right_wing'])
wing.add_shape_parameter('aspect_ratio')
# twist will be left as a hanging input, and it will be up to the user to create it as a dv or input upstream.
# user can also specify val=  if they just want to specify it as a constant input.
wing.add_shape_parameter('twist', num=10, order=4)
wing.add_shape_parameter('shape', num_spanwise=6, num_chordwise=10, order_spanwise=4, order_chordwise=4)
design_geometry.add_component(wing)

rotor1 = lsdo.RotorComponent('rotor1', stp_entity_names=['rotor1'])
rotor1.add_shape_parameter('radius')
design_geometry.add_component(rotor1)

# Insert code here to compute rib pointsets
ribs = lsdo.Feature('ribs')
for rib_pointset in rib_pointsets:
    ribs.add_pointset(rib_pointset)
design.add_feature(ribs)

# TransitionSimulation.tilt_actuation will be left as a hanging input, and it will be up to the user to create it as a dv or input upstream.
# user can also specify val=  if they just want to specify it as a constant input.
tilt_actuation = lsdo.Actuation(name='tilt_actuation', pointset1=wing_root_quarter_chord, pointset2=wing_tip_quarter_chord, components=[wing, rotor1], child_actuations=[rotor_actuation1, rotor_actuation...])
design.add_actuation(tilt_actuation) # in csdl, we connect pointset1 and poinset2 (and normalize) as the axis vector for the quaternion


cruise = CruiseSimulation(altitude_km=1., cruise_speed=65)
lsdo_problem.add_simulation(cruise)
hover = HoverSimulation(tilt_actuation, angle=90)
nt = 50
transition = TransitionSimulation(tilt_actuation, angle=np.linspace(90,0,nt))

vlm_solver = VLMSolver()    # get inputs from solver developer (Jiayao)

# Insert code here to compute vlm mesh pointsets
vlm_solver.mesh = vlm_mesh = VLMMesh()
for vlm_pointset in vlm_pointsets:
    vlm_mesh.add_pointset(vlm_pointset)

cruise.add_solver(vlm_solver)   
transition.add_solver(vlm_solver)  # this can work in python, but in CSDL, we'll make 2 vlm solver instances


sim = csdl.Simulator()

# Insert code to define my preprocessing model
sim.model.add(MyPreprocessingModel())

lsdo_model = lsdo_problem.setup()  # We need to implement this - it instantiates the LSDOProblemModel CSDL model 
sim.model.add(lsdo_model)

# Insert code to define my postprocessing model
sim.model.add(MyPostprocessingModel())

sim.run()


##########################################################################################
import numpy as np
import lsdo_vehicle.api as lv
ac.set_geo(geo=geo)
# short form

''' 
The number can be a part of the stp_name argument. It is the way the stp file organizes 
the different surfaces. 

alias : a way to rename the stp_name in order to allow for easier naming conventions of hte different 
surfaces

--------

In the different specific components, have a clever underlying algorithm that can create the 
FFD block for the user. 

In the general component class, the user will have to provide more information on how that 
FFD block is created. 

'''
# boom_right = ac.create_component(stp_entity_names=['InnerLiftNacelles'], 0, alias='boom_right')

# create_liftingsurface 
# create_fuselage 
# create_rotor
# create_motor etc... 


# form form
boom_right = lg.Component('InnerLiftNacelles', 0)
boom_right.add_ffd(LiftingSurfaceFFD(), name='default') #could consider set_ffd instead
ac.add_component(boom_right)

fuse = lg.Component('Fuselage', 0)
ac.add_component(fuse)

path_name = '../../examples/CAD/'
file_name = 'rect_wing.stp'
geo = lg.DesignGeometry(path_name + file_name)
ac.set_geometry(geo)

ac.components['boom_right']
ac.components['InnerLightNacelles', 0]



# parameter_name, property_name, order, num
boom_right.ffd.add_shape_parameter('my_span_var', 'trans_y', order=2, num=2)

kin_con = Rigid(
    boom_right, boom_right_pt_set,
    fuse, fuse_pt_set,
)
ac.add_kinematic_constraint(kin_con)

cr = Cruise(
    speed=150,
    altitude=3000.,
)
ac.add_operating_condition(cr)

uvlm = UVLM(
    wake_type='free',
    mesh=uvlm_mesh,
)
ac.add_solver(uvlm)

model = ac.assemble()

s = Simulator(model)
s.run()


# everything above could be
ac = lg.Vehicle()
geo = get_ecrm2_geometry() # this function imports the stp files, adds all components
structural_geometry = get_ecrm2_structural_geometry_1() # structures mesh
cabin_geometry = get_ecrm2_cabin_geometry_layout_2()    # payload/dyanmics mesh potentially
operating_conditions = get_ecrm2_75nm_mission()         # Operating conditions

vlm = get_uvlm()
ac.add_solver(vlm)

model = ac.assemble()
s = Simulator(model)
s.run()
