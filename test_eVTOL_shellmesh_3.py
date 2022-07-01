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

'''Just test the path of gurobi licence file'''
import gurobipy as gp
test = gp.Model('test') # export GRB_LICENSE_FILE=/Users/Sansara/Public/Code/A/gurobi.lic

from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh_1 import geo, members_ctrl_pointset_list

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
OML_ctrl_pointset_list  = shell_mesh.merge_OML(geo, merged_OML_list) 
vd_points2 = []
print()
for pointset in OML_ctrl_pointset_list:
    print(pointset.pointset_id,pointset.name)
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    vd_points2.append(vedo.Points(points, r=20, c='blue',alpha=0.3))
print()

num_u_wing = 50
num_v_wing = 30
num_u_spar = 42
num_v_spar = 4
num_u_rib = 12
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

'''            
            if pointset1.name != 'ctrl_pts_OML_lower_wing':
                print(pointset1.name)
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [435,405,406], [404,405,435]], dtype = np.int32), axis = 0)#
                #plot = 1
            else:
                print(pointset1.name)
                #plot = 1
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [404,405,406], [404,435,406]], dtype = np.int32), axis = 0)  #'''
shell_mesh.identify_intersection_list(geo, intersection_list_upper_wing)

shell_mesh.identify_intersection_list(geo, intersection_list_primary_spar)
shell_mesh.identify_intersection_list(geo, intersection_list_rear_spar)

shell_mesh.identify_intersection_list(geo, intersection_list_lower_wing)
shell_mesh.identify_intersection_list(geo,intersection_list_wing)
'''
            if memb.options['id'] == 5:
                print(inde[0])
                print(cons)
                print(inde[cons])
                inde[0]= 404
                memb.options['node_indices'] = list(inde)
            if memb.options['id'] == 11:
                print(inde[0])
                print(cons)
                print(inde[cons])                
                inde[0]= 1095
                memb.options['node_indices'] = list(inde)'''
shell_mesh.construct_whole_structure_mesh(plot = True)#

#shell_mesh.save_tri_vtk('CAD_test_eVTOL_shellmesh_1', shell_mesh.total_points, shell_mesh.tri_connectivity)

shell_mesh.optimizie_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_test_eVTOL_shellmesh_3')

exit()

'''
        connectivity_check_list = list(np.copy(connectivity_check).flatten()) 
        if np.amax(connectivity_check) >= len(self.members_dict[pointset0.name].options['u_v_vec']):
            print('test_eVTOL_shellmesh_3')
            m = max(connectivity_check_list)
            max_flatten_indices_list = [i for i, j in enumerate(connectivity_check_list) if j == m]
            max_indices_list = []
            for i in max_flatten_indices_list:
                print(i, (i)//3, connectivity_check[(i)//3,:]) 
                max_indices_list.append((i)//3)
                #print()
            connectivity_check = np.delete(connectivity_check, max_indices_list,axis = 0)
            print(len(connectivity_check))
            if pointset1.name != 'ctrl_pts_OML_lower_wing':
                print(pointset1.name)
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [435,405,406], [404,405,435]], dtype = np.int32), axis = 0)#
            else:
                print(pointset1.name)
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [404,405,406], [404,435,406]], dtype = np.int32), axis = 0)  #
            self.members_dict[pointset_ini.name].options['tri_connectivity'] = connectivity_check'''