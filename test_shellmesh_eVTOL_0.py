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
from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import *

from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import geo, structures_mesh, members_ctrl_pointset_list

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

merged_OML_list = [#[pointset_wing_tip0, pointset_wing_tip1, 'wing_tip'],
                    [142, 144, 'OML_upper_wing'],
                    [141, 143, 'OML_lower_wing'],                   
                    ] 
OML_ctrl_pointset_list  = shell_mesh.merge_OML(geo, merged_OML_list) 
vd_points2 = []
print()
for pointset in OML_ctrl_pointset_list:
    print(pointset.pointset_id,pointset.name)
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    vd_points2.append(vedo.Points(points, r=20, c='blue',alpha=0.3))
print()

intersection_list_upper_wing = []
intersection_list_lower_wing = []
intersection_list_primary_spar = []
intersection_list_rear_spar = []
for pointset in members_ctrl_pointset_list:
    if 'spar' in pointset.name:
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = 67, 40, 56, 4 #65, 40, 56, 4
        intersection_list_upper_wing.append([147,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u0','T'])
        intersection_list_lower_wing.append([150,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T'])
    elif 'rib' in pointset.name:
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = 67, 40, 15, 4 #65, 40, 15, 4
        intersection_list_upper_wing.append([147,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u0','T']) 
        intersection_list_lower_wing.append([150,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T']) 
        num_points_u0, num_points_v0, num_points_u1, num_points_v1 = 56, 4, 15, 4      
        intersection_list_primary_spar.append([127,pointset.pointset_id,num_points_u0, num_points_v0, num_points_u1, num_points_v1,'v0','T'])
        intersection_list_rear_spar.append([128,pointset.pointset_id,num_points_u0, num_points_v0, num_points_u1, num_points_v1,'v1','T'])


# intersection_list_upper_wing_test = [intersection_list_upper_wing[0],intersection_list_upper_wing[1],intersection_list_upper_wing[3]]
shell_mesh.identify_intersection_list(geo, intersection_list_upper_wing)
#shell_mesh.identify_intersection_list(geo, intersection_list_lower_wing)
#shell_mesh.identify_intersection_list(geo, [intersection_list_primary_spar[0], intersection_list_primary_spar[1], intersection_list_primary_spar[2], intersection_list_primary_spar[3], intersection_list_primary_spar[4], intersection_list_primary_spar[5], intersection_list_primary_spar[6]], plot = True)#
shell_mesh.identify_intersection_list(geo, intersection_list_primary_spar)
shell_mesh.identify_intersection_list(geo, intersection_list_rear_spar)
#shell_mesh.identify_intersection_list(geo, intersection_list_lower_wing,plot = True)
shell_mesh.identify_intersection_list(geo, [intersection_list_lower_wing[0]])
#shell_mesh.identify_intersection_list(geo, [intersection_list_rear_spar[0], intersection_list_rear_spar[1], intersection_list_rear_spar[2], intersection_list_rear_spar[3], intersection_list_rear_spar[4], intersection_list_rear_spar[5]], plot = True)
shell_mesh.construct_whole_structure_mesh(plot = True)
exit()



vd_points3 = vedo.Points(points, r=10, c='black',alpha=1.0)
vd_points4 = vedo.Points(pts, r=15, c='red',alpha=0.5)
print(len(points))
length = len(B['triangles'].tolist())
mesh = vedo.Mesh([points, np.array(B['triangles'].tolist()).reshape((length,3))])
mesh.backColor().lineColor('green').lineWidth(3)
vd_test = vedo.Plotter(axes=1)#, mesh
vd_test.show(vd_points2, 'Test', viewup="z", interactive=True) #vd_points0, vd_points1, #vd_points2, vd_points3, vd5, mesh




temp = vedo.Plotter()
temp.show(interactive=True)

shell_mesh.transform_to_quad()
shell_mesh.write_vtk('eVTOL.vtk')




#TODO change name 
shell_mesh.add_design_geometry(geo)

top_wing_surface_names = [ 
    'Surface, 1', 
    'Surface, 3', 
    ]
bot_wing_surface_names = [
    'Surface, 0',
    'Surface, 2', 
    ]

temp = vedo.Plotter()
temp.show(interactive=True)

shell_mesh.identify_intersection(intersection_list)
shell_mesh.create_triangulation()
shell_mesh.improve_mesh_quality(weight_parameters)
shell_mesh.transform_to_quad()
shell_mesh.write_vtk('eVTOL.vtk')
shell_mesh.improve_mesh_quality(weight_parameters)


for i, pointset in enumerate(OML_pointset_list):
    print(pointset.pointset_id,pointset.name)
    if pointset.name == 'OML: Surface, {}'.format(i):
        pass

    # if pointset.name == 'OML: Surface, 1':
    #     pointset_upper_wing0 = pointset
    # if pointset.name == 'OML: Surface, 3':
    #     pointset_upper_wing1 = pointset         
    # if pointset.name == 'OML: Surface, 0':
    #     pointset_lower_wing0 = pointset
    # if pointset.name == 'OML: Surface, 2':
    #     pointset_lower_wing1 = pointset      
    # if pointset.name == 'OML: Surface, 4':
    #     pointset_wing_tip0 = pointset
    #     geo.assemble(pointset = pointset)
    #     points = geo.evaluate(pointset = pointset)
    #     vd_points0 = vedo.Points(points, r=10, c='red',alpha=0.8)          
    # if pointset.name == 'OML: Surface, 5':
    #     pointset_wing_tip1 = pointset
    #     geo.assemble(pointset = pointset)
    #     points = geo.evaluate(pointset = pointset)
    #     vd_points1 = vedo.Points(points, r=15, c='green',alpha=0.5) 