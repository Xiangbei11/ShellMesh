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

shell_mesh = ShellMesh('shell_mesh', members_ctrl_pointset_list)
bspline_surface_list = list(geo.input_bspline_entity_dict.values())
OML_pointset_list = shell_mesh.extract_pointset_list_from_bspline_surface(geo, bspline_surface_list)

for pointset in OML_pointset_list:
    print(pointset.pointset_id,pointset.name)
    if pointset.name == 'OML: Surface, 1':
        pointset_upper_wing0 = pointset
      
    if pointset.name == 'OML: Surface, 3':
        pointset_upper_wing1 = pointset
         
    if pointset.name == 'OML: Surface, 0':
        pointset_lower_wing0 = pointset
    if pointset.name == 'OML: Surface, 2':
        pointset_lower_wing1 = pointset
    if pointset.name == 'OML: Surface, 4':
        pointset_wing_tip0 = pointset
        geo.assemble(pointset = pointset)
        points = geo.evaluate(pointset = pointset)
        vd_points0 = vedo.Points(points, r=10, c='red',alpha=0.8)          
    if pointset.name == 'OML: Surface, 5':
        pointset_wing_tip1 = pointset
        geo.assemble(pointset = pointset)
        points = geo.evaluate(pointset = pointset)
        vd_points1 = vedo.Points(points, r=15, c='green',alpha=0.5) 

merged_OML_list = [[pointset_wing_tip0, pointset_wing_tip1, 'wing_tip'],
                    [pointset_upper_wing0, pointset_upper_wing1, 'upper_wing'],
                    [pointset_lower_wing0, pointset_lower_wing1, 'lower_wing'],                   
                    ] 
OML_ctrl_pointset_list  = shell_mesh.merge_OML(geo, merged_OML_list) 
vd_points2 = []
for pointset in OML_ctrl_pointset_list:
    print(pointset.pointset_id,pointset.name, pointset.shape)
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    vd_points2.append(vedo.Points(points, r=20, c='blue',alpha=0.3))
vd_test = vedo.Plotter(axes=1)#, mesh
vd_test.show(vd_points2, 'Test', viewup="z", interactive=True) #vd_points0, vd_points1 #vd_points2, vd_points3, vd5, mesh
exit()


# for pointset in shell_mesh.pointest_dict.values():
#     print(pointset.pointset_id,pointset.name)
intersection_list =[[77, 64, 'T'],\
                    [77, 67, 'T'],\
                    #[77, 69, 'T'],\
                    # ['ctrl_pts_upper_wing', 'ctrl_pts_lower_wing', '-'],\
                        ] #'', '', '+' #'', '', '-'

points, pts, B, vd5 = shell_mesh.identify_intersection(geo, intersection_list)

vd_points3 = vedo.Points(points, r=10, c='black',alpha=1.0)
vd_points4 = vedo.Points(pts, r=15, c='red',alpha=0.5)
print(len(points))
length = len(B['triangles'].tolist())
mesh = vedo.Mesh([points, np.array(B['triangles'].tolist()).reshape((length,3))])
mesh.backColor().lineColor('green').lineWidth(3)





pointset_list = oml_pointset_list + structures_mesh.pointset_list

for pointset in geo.pointsets_dict.values():
    print(pointset.pointset_id,pointset.name, type(pointset).__name__)
print()
for pointset in pointset_list:
    print(pointset.pointset_id,pointset.name)
print(geo.pointsets_dict[1])
print(geo.pointsets_dict[1].physical_coordinates)
geo.assemble(pointset = geo.pointsets_dict[1])
points = geo.evaluate(pointset = geo.pointsets_dict[1])
print(points)



intersection_list =[['OML: Surface, 0', 'OML: Surface, 2', 'T'],\
                    ['', '', 'T'],\
                    ['', '', 'T']] #'', '', '+' #'', '', '-'
shell_mesh.identify_intersection(intersection_list)
temp = vedo.Plotter()
temp.show(interactive=True)

shell_mesh.create_triangulation()
shell_mesh.transform_to_quad()
shell_mesh.write_vtk('eVTOL.vtk')



# for pointset in structures_mesh.pointset_list:
#     print(pointset.pointset_id,pointset.name)

#whole_structure_pointset_list = pointsets_dict



shell_mesh = ShellMesh('shell_mesh', structures_mesh.pointset_list)

''' shell_mesh.extract_pointsets_from_OML()'''

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



