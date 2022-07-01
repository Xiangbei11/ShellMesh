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
# print()
# for i, pointset in enumerate(OML_ctrl_pointset_list_wingbox_side_03):
#     print(pointset.pointset_id,pointset.name)
# print()
merged_OML_list_wingbox_side_04 = [[327, 329, 'OML_wingbox_side43']]
OML_ctrl_pointset_list_wingbox_side_04  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_04)
print(OML_ctrl_pointset_list_wingbox_side_04[0].pointset_id)
pointset_331 = OML_ctrl_pointset_list_wingbox_side_03[2]
pointset_331.shape= np.array([121, 11, 3])
ctrl_pointset_331 = geo.fit_bspline_ctrl_pointsets([pointset_331])
print(ctrl_pointset_331[0].pointset_id,ctrl_pointset_331[0].name)
merged_OML_list_wingbox_side_05 = [[334, 335, 'OML_wingbox_side44']]
OML_ctrl_pointset_list_wingbox_side_05  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_05)
for pointset in OML_ctrl_pointset_list_wingbox_side_05:
    print(pointset.pointset_id, pointset.name)



print()
print()
print('Another side')
merged_OML_list_wingbox_side_10 = []
for i in range(22):
    merged_OML_list_wingbox_side_10.append([2*i+59, 2*i+60, 'OML_wingbox_side'+'%s'%(i+45)])
OML_ctrl_pointset_list_wingbox_side_10  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_10) 
# for pointset in OML_ctrl_pointset_list_wingbox_side_10:
#     print(pointset.pointset_id, pointset.name)
merged_OML_list_wingbox_side_11 = []
for i in np.arange(0,22,2):
    merged_OML_list_wingbox_side_11.append([2*i+340, int(2*(i-1)+344), 'OML_wingbox_side'+'%s'%int(((i/2)+67))])
OML_ctrl_pointset_list_wingbox_side_11  = shell_mesh.merge_OML(geo, [merged_OML_list_wingbox_side_11[0]], plot =True) 
print('START')
for pointset in OML_ctrl_pointset_list_wingbox_side_11:
    print(pointset.pointset_id, pointset.name)

merged_OML_list_wingbox_side_12 = []
for i in np.arange(0,8,2):
    merged_OML_list_wingbox_side_12.append([2*i+387, int(2*(i-1)+391), 'OML_wingbox_side'+'%s'%int(((i/2)+78))])
merged_OML_list_wingbox_side_12.append([403, 406, 'OML_wingbox_side'+'%s'%82])
print()
for i in merged_OML_list_wingbox_side_12:
    print(i)
print()

OML_ctrl_pointset_list_wingbox_side_12  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_12[:2])

vd_points2 = [] 
print()
for i, pointset in enumerate(OML_ctrl_pointset_list_wingbox_side_05):
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
# exit()