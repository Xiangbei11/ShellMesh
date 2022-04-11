import csdl
import numpy as np

from parametrization.section_properties_model import SectionPropertiesModel 
from lsdo_kit.geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

from csdl_om import Simulator

nxp = 5
nyp = 5
nzp = 5

point000 = np.array([170. ,0. ,100.])
point010 = np.array([130., 230., 100.])
point001 = np.array([170., 0., 170.])
point011 = np.array([130., 230., 170.])

point100 = np.array([240. ,0. ,100.])
point101 = np.array([240. ,0. ,170.])
point110 = np.array([200. ,230. ,100.])
point111 = np.array([200. ,230. ,170.])

control_points = np.zeros((2,2,2,3))

control_points[0,0,0,:] = point000
control_points[0,0,1,:] = point001

control_points[0,1,0,:] = point010
control_points[0,1,1,:] = point011

control_points[1,0,0,:] = point100
control_points[1,0,1,:] = point101

control_points[1,1,0,:] = point110
control_points[1,1,1,:] = point111

ffd_control_points1 = create_ffd(control_points, nxp, nyp, nzp)

point000 = np.array([90. ,0. ,100.])
point010 = np.array([180., 230., 100.])
point001 = np.array([90., 0., 170.])
point011 = np.array([180., 230., 170.])

point100 = np.array([130. ,0. ,100.])
point101 = np.array([130. ,0. ,170.])
point110 = np.array([70. ,230. ,100.])
point111 = np.array([70. ,230. ,170.])

control_points = np.zeros((2,2,2,3))

control_points[0,0,0,:] = point000
control_points[0,0,1,:] = point001

control_points[0,1,0,:] = point010
control_points[0,1,1,:] = point011

control_points[1,0,0,:] = point100
control_points[1,0,1,:] = point101

control_points[1,1,0,:] = point110
control_points[1,1,1,:] = point111

ffd_control_points2 = create_ffd(control_points, nxp, nyp, nzp)

ffd_block_1 = FFD('test_ffd1', ffd_control_points1)
ffd_block_2 = FFD('test_ffd2', ffd_control_points2)



'''
Go through the FFD classes and add an 'add' method 

ADD METHOD options:
name 
dv
parameter_name 
order 
num_cp
val : optional, if val is none then declare variable

'''

# Create ffd_blocks list, and the csdl variables needed
# ffd_block_1.add_shape_parameter(property_name = 'rot_x', parameter_names=['linear'], order=[2], num_cp=[3], dv=[False])
# TODO: Make each function call add only one parameter.
ffd_block_1.add_shape_parameter(property_name = 'rot_y', parameter_names=['linear', 'quadratic'], order=[2, 3], num_cp=[3, 4], dv=[True, False], val=[np.arange(3), None])

# ffd_block_2.add_shape_parameter(property_name ='rot_x', parameter_names=['linear', 'quadratic'], order=[2,3], num_cp=[3,5], dv=[False, False], val=[np.arange(3), np.arange(5)])

ffd_blocks = [ffd_block_1]

sim = Simulator(SectionPropertiesModel(ffd_blocks=ffd_blocks))
sim.run()
sim.prob.model.list_inputs(prom_name=True)
sim.prob.model.list_outputs(prom_name=True)

sim.visualize_model()

