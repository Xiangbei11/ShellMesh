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
bspline_surface_list = list(geo.input_bspline_entity_dict.values())
OML_pointset_list = shell_mesh.extract_pointset_list_from_bspline_surface(geo, bspline_surface_list)
for i, pointset in enumerate(OML_pointset_list):
    print(pointset.pointset_id,pointset.name)
print(geo.current_id)

merged_OML_list_wingbox_side_00 = []
for i in range(24):
    merged_OML_list_wingbox_side_00.append([2*i+11, 2*i+12, 'OML_wingbox_side'+'%s'%i])
'''247 ctrl_pts_OML_wingbox_side0
249 ctrl_pts_OML_wingbox_side1'''
OML_ctrl_pointset_list_wingbox_side_00  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_00) 

merged_OML_list_wingbox_side_01 = []
for i in np.arange(2,24,2):
    merged_OML_list_wingbox_side_01.append([2*i+247, int(2*(i-1)+251), 'OML_wingbox_side'+'%s'%int(((i/2)+23))])
'''315 ctrl_pts_OML_wingbox_side34'''
OML_ctrl_pointset_list_wingbox_side_01  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_01) 

merged_OML_list_wingbox_side_02 = []
for i in np.arange(0,10,2):
    merged_OML_list_wingbox_side_02.append([2*i+295, int(2*(i-1)+299), 'OML_wingbox_side'+'%s'%int(((i/2)+35))])
OML_ctrl_pointset_list_wingbox_side_02  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_02)

merged_OML_list_wingbox_side_03 = [[317, 319, 'OML_wingbox_side40'], [321, 323, 'OML_wingbox_side41'],[325, 315, 'OML_wingbox_side42']]
OML_ctrl_pointset_list_wingbox_side_03  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_03)

merged_OML_list_wingbox_side_04 = [[327, 329, 'OML_wingbox_side43']]
OML_ctrl_pointset_list_wingbox_side_04  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_04,plot = True)
vd_points2 = [] 
for i, pointset in enumerate(OML_ctrl_pointset_list_wingbox_side_04):
    print(pointset.pointset_id,pointset.name)
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    color = list(vedo.colors.colors.values())[i]
    if i%2==0:
        vd_points2.append(vedo.Points(points, r=20, c=color, alpha=0.5)) 
    else:
        vd_points2.append(vedo.Points(points, r=15, c=color, alpha=0.8))

vd_test = vedo.Plotter(axes=1)
vd_test.show(vd_points2, 'Test', viewup="z", interactive=True)

# vd = vedo.Plotter()
# vd.show(interactive = True) 
exit()