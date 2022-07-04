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

from lsdo_kit.design.design_geometry.examples.test_uCRM_wingbox_shellmesh_0 import geo#, members_ctrl_pointset_list

os.chdir(cwd)
from shellmesh import ShellMesh

shell_mesh = ShellMesh('shell_mesh')
bspline_surface_list = list(geo.input_bspline_entity_dict.values())[:120]
OML_pointset_list = shell_mesh.extract_pointset_list_from_bspline_surface(geo, bspline_surface_list)
for i, pointset in enumerate(OML_pointset_list):
    print(pointset.pointset_id,pointset.name)
print(geo.current_id)
merged_OML_list = [[15, 16, 'OML_wingbox_side44']]
OML_ctrl_pointset_list  = shell_mesh.merge_OML(geo, merged_OML_list, plot = True)

num_points_u0, num_points_v0, num_points_u1, num_points_v1 = 10, 4, 10, 4
print(OML_pointset_list[14].shape,OML_pointset_list[111].shape)
OML_pointset_list[14].shape = np.array([11, 11, 3])
OML_pointset_list[111].shape = np.array([11, 11, 3])
intersection_list = [[OML_pointset_list[14].pointset_id,OML_pointset_list[111].pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'v1','T']]
shell_mesh.identify_intersection_list(geo, intersection_list, plot =True)
vd = vedo.Plotter()
vd.show(interactive = True) 
exit()
