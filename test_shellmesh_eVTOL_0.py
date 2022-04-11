import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit/lsdo_kit/design/design_geometry/examples')
import os
cwd = os.getcwd() 
CAD_file_path = cwd + '/lsdo_kit/lsdo_kit/design/design_geometry/examples'
os.chdir(CAD_file_path)
from lsdo_kit.design.design_geometry.examples.test_eVTOL_wing_structure_shellmesh import *

os.chdir(cwd)
print(os.getcwd())
from shellmesh import ShellMesh

import vedo

# for pointset in structures_mesh.pointset_list:
#     print(pointset.pointset_id,pointset.name)

#whole_structure_pointset_list = pointsets_dict

geo.extract_pointsets_from_OML()

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
exit()



