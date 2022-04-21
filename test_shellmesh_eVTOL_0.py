import vedo

import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit/lsdo_kit/design/design_geometry/examples')
import os
cwd = os.getcwd() 
CAD_file_path = cwd + '/lsdo_kit/lsdo_kit/design/design_geometry/examples'
os.chdir(CAD_file_path)
from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import *
from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import geo, structures_mesh, members_ctrl_pointsets 

os.chdir(cwd)
from shellmesh import ShellMesh
shell_mesh = ShellMesh('shell_mesh', members_ctrl_pointsets)
oml_pointset_list = geo.extract_pointset_list_from_OML()
for pointset in oml_pointset_list:
    print(pointset.pointset_id,pointset.name, type(pointset).__name__)
    if pointset.pointset_id == 72:
        geo.assemble(pointset = pointset)
        points = geo.evaluate(pointset = pointset)
        vd_points0 = vedo.Points(points, r=10, c='red',alpha=0.8)
        pointset_72 = pointset

    if pointset.pointset_id == 74:
        geo.assemble(pointset = pointset)
        points = geo.evaluate(pointset = pointset)
        vd_points1 = vedo.Points(points, r=15, c='green',alpha=0.5)
        pointset_74 = pointset

# intersection_list =[[pointset_72, pointset_74, '-'],\
#                     ['', '', 'T'],\
#                     ['', '', 'T']] 
# pointset_72_74,intersection = shell_mesh.identify_intersection(geo, intersection_list)

pointset_72_74,intersection = shell_mesh.merge_OML(geo, [[pointset_72, pointset_74, '-']])
geo.assemble(pointset = pointset_72_74)
points = geo.evaluate(pointset = pointset_72_74)
vd_points2 = vedo.Points(points, r=20, c='blue',alpha=0.3) 
vd_points3 = vedo.Points(intersection, r=20, c='black',alpha=0.9) 
print('pointset_72_74', pointset_72_74.shape)
vd_test = vedo.Plotter(axes=1)
vd_test.show(vd_points0, vd_points1, vd_points2,  'Test', viewup="z", interactive=True)
exit()


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



