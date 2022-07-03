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
print()

merged_OML_list_wingbox_side_00 = []
for i in range(24):
    merged_OML_list_wingbox_side_00.append([2*i+11, 2*i+12, 'OML_wingbox_side'+'%s'%i])
OML_ctrl_pointset_list_wingbox_side_00  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_00, num=11) 
merged_OML_list_wingbox_side_01 = []
for i in np.arange(2,24,2):
    merged_OML_list_wingbox_side_01.append([2*i+247, int(2*(i-1)+251), 'OML_wingbox_side'+'%s'%int(((i/2)+23))])
'''315 ctrl_pts_OML_wingbox_side34'''
OML_ctrl_pointset_list_wingbox_side_01  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_01, num=11) 
merged_OML_list_wingbox_side_02 = []
for i in np.arange(0,10,2):
    merged_OML_list_wingbox_side_02.append([2*i+295, int(2*(i-1)+299), 'OML_wingbox_side'+'%s'%int(((i/2)+35))])
OML_ctrl_pointset_list_wingbox_side_02  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_02, num=11)
merged_OML_list_wingbox_side_03 = [[317, 319, 'OML_wingbox_side40'], [321, 323, 'OML_wingbox_side41'],[325, 315, 'OML_wingbox_side42']]
OML_ctrl_pointset_list_wingbox_side_03  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_03, num=11)
merged_OML_list_wingbox_side_04 = [[327, 329, 'OML_wingbox_side43']]
OML_ctrl_pointset_list_wingbox_side_04  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_04, num=11)
# print(OML_ctrl_pointset_list_wingbox_side_04[0].pointset_id)
pointset_331 = OML_ctrl_pointset_list_wingbox_side_03[2]
pointset_331.shape= np.array([121, 11, 3])
ctrl_pointset_331 = geo.fit_bspline_ctrl_pointsets([pointset_331])
merged_OML_list_wingbox_side_05 = [[334, 335, 'OML_wingbox_side44']]
OML_ctrl_pointset_list_wingbox_side_05  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_05, num=11)


print()
print()
print('Another side')
merged_OML_list_wingbox_side_10 = []
for i in range(22):
    merged_OML_list_wingbox_side_10.append([2*i+59, 2*i+60, 'OML_wingbox_side'+'%s'%(i+45)])
merged_OML_list_wingbox_side_10.append([382, 103, 'OML_wingbox_side'+'%s'%67])
OML_ctrl_pointset_list_wingbox_side_10  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_10, num=11) 

merged_OML_list_wingbox_side_11 = []
for i in np.arange(0,20,2):
    merged_OML_list_wingbox_side_11.append([2*i+340, int(2*(i-1)+344), 'OML_wingbox_side'+'%s'%int(((i/2)+68))])
merged_OML_list_wingbox_side_11.append([380, 384, 'OML_wingbox_side'+'%s'%78])
OML_ctrl_pointset_list_wingbox_side_11  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_11, num=11) 
'''384 ctrl_pts_OML_wingbox_side67'''
merged_OML_list_wingbox_side_12 = []
for i in np.arange(0,10,2):
    merged_OML_list_wingbox_side_12.append([2*i+388, int(2*(i-1)+392), 'OML_wingbox_side'+'%s'%int(((i/2)+79))])
for i in merged_OML_list_wingbox_side_12:
    print(i)
OML_ctrl_pointset_list_wingbox_side_12  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_12, num=11)
merged_OML_list_wingbox_side_13 = [[408, 410, 'OML_wingbox_side84'], [412, 414, 'OML_wingbox_side85']]
OML_ctrl_pointset_list_wingbox_side_13  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_13, num=11)
merged_OML_list_wingbox_side_14 = [[418, 420, 'OML_wingbox_side86']]
OML_ctrl_pointset_list_wingbox_side_14  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_14, num=11)
pointset_416 = OML_ctrl_pointset_list_wingbox_side_12[4]
pointset_416.shape= np.array([91, 11, 3])
ctrl_pointset_416 = geo.fit_bspline_ctrl_pointsets([pointset_416])
merged_OML_list_wingbox_side_15 = [[423, 424, 'OML_wingbox_side87']]
OML_ctrl_pointset_list_wingbox_side_15 = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_15, num=11)

side_spar_0 = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(num_points_u):
    if i== 0 or i==34 or i==47 or i==2:
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve0',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve0 = geo.extract_pointset(side_spar_0, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve0)
points_side_curve0 = geo.evaluate(pointset = side_curve0)

side_spar_1 = OML_ctrl_pointset_list_wingbox_side_15[0]
num_points_u = 45
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(num_points_u):
    if i==43:
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve1',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve1 = geo.extract_pointset(side_spar_1, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve1)
points_side_curve1 = geo.evaluate(pointset = side_curve1)


side_curve = OML_pointset_list[157]
geo.assemble(pointset = side_curve)
points_side_curve = geo.evaluate(pointset = side_curve)
indices = np.arange(-10,-1)
print('side_curve2',len(indices))
side_curve2 = geo.extract_pointset(side_curve, indices, len(indices))
geo.assemble(pointset = side_curve2)
side_curve2_points = geo.evaluate(pointset = side_curve2)

side_curve = OML_pointset_list[195]
geo.assemble(pointset = side_curve)
points_side_curve = geo.evaluate(pointset = side_curve)
indices = np.arange(1,10)+22
print('side_curve3',len(indices))
side_curve3 = geo.extract_pointset(side_curve, indices, len(indices))
geo.assemble(pointset = side_curve3)
side_curve3_points = geo.evaluate(pointset = side_curve3)

top_point1 = geo.extract_pointset(side_curve0, np.array([-1]), np.array([1]))
top_point2 = geo.extract_pointset(side_curve0, np.array([-22]), np.array([1]))
bot_point1 = geo.extract_pointset(side_curve1, np.array([-1]), np.array([1]))
bot_point2 = geo.extract_pointset(side_curve1, np.array([-22]), np.array([1]))
geo.assemble(pointset = top_point1)
top_point1_vd = geo.evaluate(pointset = top_point1)
geo.assemble(pointset = top_point2)
top_point2_vd = geo.evaluate(pointset = top_point2)
geo.assemble(pointset = bot_point1)
bot_point1_vd = geo.evaluate(pointset = bot_point1)
geo.assemble(pointset = bot_point2)
bot_point2_vd = geo.evaluate(pointset = bot_point2)
# side_curve2 = geo.add_pointsets(top_point1, side_curve2)
# side_curve3 = geo.add_pointsets(top_point2, side_curve3)
# side_curve2 = geo.add_pointsets(side_curve2, bot_point1)
# side_curve3 = geo.add_pointsets(side_curve3, bot_point1)
# geo.assemble(pointset = side_curve2)
# side_curve2_points = geo.evaluate(pointset = side_curve2)
# geo.assemble(pointset = side_curve3)
# side_curve3_points = geo.evaluate(pointset = side_curve3)
# num_side_pts_spars = 10
# rear_side_root_curve = geo.perform_linear_interpolation(top_point1, bot_point1, np.array([num_side_pts_spars]))
# rear_side_tip_curve = geo.perform_linear_interpolation(top_point2, bot_point2, np.array([num_side_pts_spars]))
# surface0 = geo.perform_2d_transfinite_interpolation(side_curve0, side_curve1, 
#     rear_side_root_curve, rear_side_tip_curve)


# ctrl_surface0 = geo.fit_bspline_ctrl_pointsets([surface0])[0]
# print('ctrl_surface0', ctrl_surface0.pointset_id,ctrl_surface0.name,ctrl_surface0.shape)
# geo.assemble(pointset = ctrl_surface0)

vd_points2 = [] 
print()
print('Plotting')
for i, pointset in enumerate(OML_ctrl_pointset_list_wingbox_side_15):   #
    geo.assemble(pointset = pointset)
    points = geo.evaluate(pointset = pointset)
    print(pointset.pointset_id,pointset.name, points.shape)
    unique_points = np.unique(points, axis=0)
    print('unique',len(points),len(unique_points))    
    color = list(vedo.colors.colors.values())[i]
    if i%2==0:
        vd_points2.append(vedo.Points(points, r=20, c=color, alpha=0.5)) 
    else:
        vd_points2.append(vedo.Points(points, r=15, c=color, alpha=0.8))
vd_points3 = vedo.Points(points_side_curve0, r=25, c='red', alpha=0.3)
vd_points4 = vedo.Points(points_side_curve1, r=25, c='red', alpha=0.3)
#vd_points5 = vedo.Points(rib0_points, r=10, c='black', alpha=0.3)#surface_points0
vd_points5 = vedo.Points(side_curve2_points, r=30, c='black', alpha=0.3)
vd_points6 = vedo.Points(side_curve3_points, r=30, c='black', alpha=0.3)
vd_points7 = vedo.Points(top_point1_vd, r=30, c='green', alpha=0.3)
vd_points8 = vedo.Points(top_point2_vd, r=30, c='green', alpha=0.3)
vd_points9 = vedo.Points(bot_point1_vd, r=30, c='blue', alpha=0.3)
vd_points10 = vedo.Points(bot_point2_vd, r=30, c='blue', alpha=0.3)
#vd_points5 = vedo.Points(surface_points0, r=10, c='black', alpha=1.0)#surface_points0
vd_test = vedo.Plotter(axes=1)
vd_test.show(vd_points3,vd_points4,vd_points5,vd_points6,vd_points7,vd_points8,vd_points9,vd_points10, 'Test15', viewup="z", interactive=True)
exit()
vd = vedo.Plotter()
vd.show(interactive = True) 

merged_OML_list_wingbox_side_21 = [[158,162, 'OML_wingbox_side85'], [165,169, 'OML_wingbox_side86']]
OML_ctrl_pointset_list_wingbox_side_21  = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_21, num=11, plot = True)
exit()
rib_list = OML_pointset_list[120:157]
for i, rib in enumerate(rib_list):
    print('rib', rib.pointset_id,rib.name,rib.shape)
    rib.shape = np.array([11,11,3])
geo.assemble(pointset = rib)
surface_points0 = geo.evaluate(pointset = rib)

num_u_upper = 45
num_v_lower = 6
num_u_rib = 6
num_v_rib = 2
intersection_list_upper =[]
for pointset in rib_list:
    num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_upper, num_v_lower, num_u_rib, num_v_rib #65, 40, 56, 4
    intersection_list_upper.append([436,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T'])
shell_mesh.identify_intersection_list(geo, intersection_list_upper, plot =True)

side_spar_0 = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(num_points_u):
    if i== 0 or i==34 or i==47 or i==2:
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve0',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve0 = geo.extract_pointset(side_spar_0, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve0)
points_side_curve0 = geo.evaluate(pointset = side_curve0)