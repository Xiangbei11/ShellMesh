import vedo
import numpy as np
import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit/lsdo_kit/design/design_geometry/examples')
import os
cwd = os.getcwd() 
CAD_file_path = cwd + '/lsdo_kit/lsdo_kit/design/design_geometry/examples'
os.chdir(CAD_file_path)
from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import geo, members_ctrl_pointset_list

os.chdir(cwd)
from shellmesh import ShellMesh

for pointset in members_ctrl_pointset_list:
    print(pointset.pointset_id,pointset.name)
print()

shell_mesh = ShellMesh('shell_mesh', members_ctrl_pointset_list)
bspline_surface_list = list(geo.input_bspline_entity_dict.values())
OML_pointset_list = shell_mesh.extract_pointset_list_from_bspline_surface(geo, bspline_surface_list)

for i, pointset in enumerate(OML_pointset_list):
    print(pointset.pointset_id,pointset.name)

merged_OML_list = [[142, 144, 'OML_upper_wing'],
                [141, 143, 'OML_lower_wing'],] 
OML_ctrl_pointset_list = shell_mesh.merge_OML(geo, [merged_OML_list[0]]) 
OML_ctrl_pointset_list.append(shell_mesh.merge_OML(geo, [merged_OML_list[1]])[0]) 

vd_points2 = []
print()
for pointset in OML_ctrl_pointset_list:
    print(pointset.pointset_id,pointset.name)
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    vd_points2.append(vedo.Points(points, r=20, c='blue',alpha=0.3))
print()

num_u_wing = 76
num_v_wing = 40
num_u_spar = 56
num_v_spar = 4
num_u_rib = 15
num_v_rib = 4
intersection_list_upper_wing = []
intersection_list_lower_wing = []
intersection_list_primary_spar = []
intersection_list_rear_spar = []

for pointset in members_ctrl_pointset_list:
    if 'spar' in pointset.name:
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_wing, num_v_wing, num_u_spar, num_v_spar #65, 40, 56, 4
        intersection_list_upper_wing.append([147,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u0','T'])
        intersection_list_lower_wing.append([150,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T'])
    elif 'rib' in pointset.name:
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_wing, num_v_wing, num_u_rib, num_v_rib #65, 40, 15, 4
        intersection_list_upper_wing.append([147,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u0','T']) 
        intersection_list_lower_wing.append([150,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T']) 
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_spar, num_v_spar, num_u_rib, num_v_rib     
        intersection_list_primary_spar.append([127,pointset.pointset_id,num_points_u0, num_points_v0, num_points_u1, num_points_v1,'v0','T'])
        intersection_list_rear_spar.append([128,pointset.pointset_id,num_points_u0, num_points_v0, num_points_u1, num_points_v1,'v1','T'])
num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_wing, num_v_wing, num_u_wing, num_v_wing
intersection_list_wing = [[147,150,num_points_u0,num_points_v0,num_points_u1,num_points_v1,['u0','u1'],'-'], [147,150,num_points_u0,num_points_v0,num_points_u1,num_points_v1,['u1','u0'],'-']]


shell_mesh.identify_intersection_list(geo, intersection_list_upper_wing)
shell_mesh.identify_intersection_list(geo, intersection_list_primary_spar)
shell_mesh.identify_intersection_list(geo, intersection_list_rear_spar)
shell_mesh.identify_intersection_list(geo, intersection_list_lower_wing)
shell_mesh.identify_intersection_list(geo,intersection_list_wing)
for memb in shell_mesh.members_dict.values():
    if memb.options['id'] == 1 or memb.options['id'] == 2: 
        print('test_eVTOL_shellmesh')
        num_u_spar = 56
        num_v_spar = 4
        indices_v0 = np.arange(num_v_spar)
        indices_v1 = np.arange(num_v_spar)+ (num_u_spar*num_v_spar) - num_v_spar  
        indices_list = [indices_v0, indices_v1]
        for indices in indices_list:
            memb.options['constrained_node_indices'] += list(indices)
            memb.options['constrained_boundary_node_indices'] += list(indices)   
            memb.options['constrained_boundary_node_indices'] = list(set(memb.options['constrained_boundary_node_indices'])) 
shell_mesh.construct_whole_structure_mesh(plot = True)#plot = True
#shell_mesh.save_tri_vtk('CAD_test_eVTOL_shellmesh_1', shell_mesh.total_points, shell_mesh.tri_connectivity)

'''Just test the path of gurobi licence file'''
import gurobipy as gp
test = gp.Model('test') # export GRB_LICENSE_FILE=/Users/Sansara/Public/Code/A/gurobi.lic

shell_mesh.optimize_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_test_eVTOL_shellmesh')#itr=[1,2,3]  #w4 = 0.7 26486

exit()
#itr=[1,2,3] #w4 = 0.9 26734 
#itr=[1,2,3] #w4 = 0.8 26570 
#itr=[1,2,2,2,2,3] #w4 = 1  27240 
#itr=[1,2,3]  #w4 = 0.7 26486
# if 'upper_wing' in memb.options['name']:
#     m.vertexCoords[:,2] = m.vertexCoords[:,2]+25
# elif 'lower_wing' in memb.options['name']:
#     m.vertexCoords[:,2] = m.vertexCoords[:,2]-25