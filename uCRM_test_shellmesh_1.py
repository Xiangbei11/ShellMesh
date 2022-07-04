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

'''Upper'''
vd_points7= []
surface00 = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(21,48):
    if i==34 or i==47 or i==22:#i== 0 or  or i==2 or i ==25
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve0',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve00 = geo.extract_pointset(surface00, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve00)
points_side_curve00 = geo.evaluate(pointset = side_curve00)

surface01 = OML_ctrl_pointset_list_wingbox_side_15[0]
num_points_u = 45
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(20,45):
    if i==43:
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve1',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve01 = geo.extract_pointset(surface01, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve01)
points_side_curve01 = geo.evaluate(pointset = side_curve01)

k = 157
print('k',k)
surface02 = OML_pointset_list[k]
geo.assemble(pointset = surface02)
side_curve_points02 = geo.evaluate(pointset = surface02)
indices = np.arange(-10,-1)
print('side_curve2',len(indices))
side_curve02 = geo.extract_pointset(surface02, indices, len(indices))
geo.assemble(pointset = side_curve02)
side_curve_points02 = geo.evaluate(pointset = side_curve02)

k = 207
print('k',k)
surface03 = OML_pointset_list[k]
geo.assemble(pointset = surface03)
side_curve_points03 = geo.evaluate(pointset = surface03)
indices = np.array([-11,-10,-9,-8,-6,-4,-3,-2,-1])
print('side_curve3',len(indices))
side_curve03 = geo.extract_pointset(surface03, indices, len(indices))
geo.assemble(pointset = side_curve03)
side_curve_points03 = geo.evaluate(pointset = side_curve03)

surface000 = geo.perform_2d_transfinite_interpolation(side_curve00, side_curve01, 
    side_curve02, side_curve03)
geo.assemble(pointset = surface000)
surface00_points = geo.evaluate(pointset = surface000)

vd_points7.append(vedo.Points(surface00_points, r=30, c='green', alpha=0.3))

surface10 = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(0,21):
    if i==0 or i==2 or i==18 or i==1:# 
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve0',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve10 = geo.extract_pointset(surface10, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve10)
points_side_curve0 = geo.evaluate(pointset = side_curve10)

surface11 = OML_ctrl_pointset_list_wingbox_side_15[0]
num_points_u = 45
num_points_v = 11
indices_u1_reduced = []
indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1#indices_u0 = num_points_v * np.arange(num_points_u)
for i in range(0,20):
    if i==0 or i==2 or i==1:
        pass
    else:
        indices_u1_reduced.append(indices_u1[i])
print('side_curve1',len(indices_u1_reduced))
indices_u1_reduced = np.array(indices_u1_reduced)
side_curve11 = geo.extract_pointset(surface11, indices_u1_reduced, len(indices_u1_reduced))
geo.assemble(pointset = side_curve11)
points_side_curve1 = geo.evaluate(pointset = side_curve11)

k = 207
print('k',k)
surface12 = OML_pointset_list[k]
geo.assemble(pointset = surface12)
side_curve_points2 = geo.evaluate(pointset = surface12)
indices = np.array([-11,-10,-9,-8,-6,-4,-3,-2,-1])#
print('side_curve2',len(indices))
side_curve12 = geo.extract_pointset(surface12, indices, len(indices))
geo.assemble(pointset = side_curve12)
side_curve_points2 = geo.evaluate(pointset = side_curve12)

k = 231
print('k',k)
surface13 = OML_pointset_list[k]
geo.assemble(pointset = surface13)
side_curve_points3 = geo.evaluate(pointset = surface13)
indices = np.array([0,1,2,3,4,6,8,9,10])
print('side_curve3',len(indices))
side_curve13 = geo.extract_pointset(surface13, indices, len(indices))
geo.assemble(pointset = side_curve13)
side_curve_points3 = geo.evaluate(pointset = side_curve13)

surface111 = geo.perform_2d_transfinite_interpolation(side_curve10, side_curve11, 
    side_curve12, side_curve13)
geo.assemble(pointset = surface111)
surface00_points = geo.evaluate(pointset = surface111)
vd_points7.append(vedo.Points(surface00_points, r=30, c='green', alpha=0.3))

'''lower'''
surface00b = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u0_reduced = []
indices_u0 = num_points_v * np.arange(num_points_u) 
for i in range(21,48):
    if i==47 or i==22 or i==34:# 
        pass
    else:
        indices_u0_reduced.append(indices_u0[i])
print('side_curve0',len(indices_u0_reduced))
indices_u0_reduced = np.array(indices_u0_reduced)
side_curve00b = geo.extract_pointset(surface00b, indices_u0_reduced, len(indices_u0_reduced))
geo.assemble(pointset = side_curve00b)
points_side_curve0 = geo.evaluate(pointset = side_curve00b)

surface01b = OML_ctrl_pointset_list_wingbox_side_15[0]
num_points_u = 45
num_points_v = 11
indices_u0_reduced = []
indices_u0 = num_points_v * np.arange(num_points_u) 
for i in range(20,45):
    if i==43:#
        pass
    else:
        indices_u0_reduced.append(indices_u0[i])
print('side_curve1',len(indices_u0_reduced))
indices_u0_reduced = np.array(indices_u0_reduced)
side_curve01b = geo.extract_pointset(surface01b, indices_u0_reduced, len(indices_u0_reduced))
geo.assemble(pointset = side_curve01b)
points_side_curve1 = geo.evaluate(pointset = side_curve01b)

k = 158
print('k',k)
surface02b = OML_pointset_list[k]
geo.assemble(pointset = surface02b)
side_curve_points2 = geo.evaluate(pointset = surface02b)
indices = np.array([-10,-9,-8,-7,-6,-5,-4,-3,-2])
print('side_curve3',len(indices))
side_curve02b = geo.extract_pointset(surface02b, indices, len(indices))
geo.assemble(pointset = side_curve02b)
side_curve_points2 = geo.evaluate(pointset = side_curve02b)

k = 208
print('k',k)
surface03b = OML_pointset_list[k]
geo.assemble(pointset = surface03b)
side_curve_points3 = geo.evaluate(pointset = surface03b)
indices = np.array([-11,-10,-9,-8,-7,-6,-5,-3,-1])-11
print('side_curve2',len(indices))
side_curve03b = geo.extract_pointset(surface03b, indices, len(indices))
geo.assemble(pointset = side_curve03b)
side_curve_points3 = geo.evaluate(pointset = side_curve03b)

surface000b = geo.perform_2d_transfinite_interpolation(side_curve00b, side_curve01b, 
    side_curve02b, side_curve03b)
geo.assemble(pointset = surface000b)
surface00_points = geo.evaluate(pointset = surface000b)
vd_points7.append(vedo.Points(surface00_points, r=30, c='green', alpha=0.3))
# print(surface111.pointset_id)

surface10b = OML_ctrl_pointset_list_wingbox_side_05[0]
num_points_u = 48
num_points_v = 11
indices_u0_reduced = []
indices_u0 = num_points_v * np.arange(num_points_u) 
for i in range(0,21):
    if i==18 or i==0 or i==1 or i ==3:# 
        pass
    else:
        indices_u0_reduced.append(indices_u0[i])
print('side_curve0',len(indices_u0_reduced))
indices_u0_reduced = np.array(indices_u0_reduced)
side_curve10b = geo.extract_pointset(surface10b, indices_u0_reduced, len(indices_u0_reduced))
geo.assemble(pointset = side_curve10b)
points_side_curve0 = geo.evaluate(pointset = side_curve10b)

surface11b = OML_ctrl_pointset_list_wingbox_side_15[0]
num_points_u = 45
num_points_v = 11
indices_u0_reduced = []
indices_u0 = num_points_v * np.arange(num_points_u) 
for i in range(0,20):
    if i==0 or i==1 or i==2 :#
        pass
    else:
        indices_u0_reduced.append(indices_u0[i])
print('side_curve1',len(indices_u0_reduced))
indices_u0_reduced = np.array(indices_u0_reduced)
side_curve11b = geo.extract_pointset(surface11b, indices_u0_reduced, len(indices_u0_reduced))
geo.assemble(pointset = side_curve11b)
points_side_curve1 = geo.evaluate(pointset = side_curve11b)

k = 208
print('k',k)
surface12b = OML_pointset_list[k]
geo.assemble(pointset = surface12b)
side_curve_points2 = geo.evaluate(pointset = surface12b)
indices = np.array([-11,-10,-9,-8,-7,-6,-5,-3,-1])-11
print('side_curve2',len(indices))
side_curve12b = geo.extract_pointset(surface12b, indices, len(indices))
geo.assemble(pointset = side_curve12b)
side_curve_points2 = geo.evaluate(pointset = side_curve12b)

k = 232
print('k',k)
surface13b = OML_pointset_list[k]
geo.assemble(pointset = surface13b)
side_curve_points3 = geo.evaluate(pointset = surface13b)
indices = np.array([0,1,2,3,4,5,6,8,10])
print('side_curve3',len(indices))
side_curve13b = geo.extract_pointset(surface13b, indices, len(indices))
geo.assemble(pointset = side_curve13b)
side_curve_points3 = geo.evaluate(pointset = side_curve13b)

surface111b = geo.perform_2d_transfinite_interpolation(side_curve10b, side_curve11b, 
    side_curve12b, side_curve13b)
geo.assemble(pointset = surface111b)
surface00_points = geo.evaluate(pointset = surface111b)
vd_points7.append(vedo.Points(surface00_points, r=30, c='green', alpha=0.3))
print(surface111.pointset_id)

# merged_OML_list_wingbox_side_20 = [[surface000.pointset_id, surface111.pointset_id, 'OML_wingbox_side88'],[surface000b.pointset_id, surface111b.pointset_id, 'OML_wingbox_side89']]
# OML_ctrl_pointset_list_wingbox_side_20 = shell_mesh.merge_OML(geo, merged_OML_list_wingbox_side_20, num=9)
# shell_mesh.sequence_points(geo, OML_ctrl_pointset_list_wingbox_side_20[0], [0,1], 40, 9, plot = True)
# OML_ctrl_pointset_list_wingbox_side_20[0].shape = np.array([40,9,3])

geo.assemble(pointset = OML_ctrl_pointset_list_wingbox_side_05[0])
points = geo.evaluate(pointset = OML_ctrl_pointset_list_wingbox_side_05[0])
vd_points4 = vedo.Points(points, r=25, c='red', alpha=0.8)
geo.assemble(pointset = OML_ctrl_pointset_list_wingbox_side_15[0])
points = geo.evaluate(pointset = OML_ctrl_pointset_list_wingbox_side_15[0])
vd_points5 = vedo.Points(points, r=25, c='blue', alpha=0.8)

# rib_list = OML_pointset_list[134:157]
# for i, rib in enumerate(rib_list):
#     rib.shape = np.array([11,11,3])
# num_u_upper = 24
# num_v_lower = 5
# num_u_rib = 5
# num_v_rib = 2
# intersection_list_upper =[]
# surface000.shape = np.array([24,9,3])
# for pointset in rib_list:
#     num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_upper, num_v_lower, num_u_rib, num_v_rib 
#     intersection_list_upper.append([surface000.pointset_id,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T'])
# shell_mesh.identify_intersection_list(geo, intersection_list_upper, plot =True)

rib_list = OML_pointset_list[113:125]+OML_pointset_list[126:130]#OML_pointset_list[112:125]+  OML_pointset_list[126:132]
print('rib_list',len(rib_list))
for i, rib in enumerate(rib_list):
    rib.shape = np.array([11,11,3])
num_u_upper = 20
num_v_lower = 5
num_u_rib = 5
num_v_rib = 2
intersection_list_upper =[]
surface111.shape = np.array([17,9,3])
for pointset in rib_list:
    num_points_u0, num_points_v0, num_points_u1, num_points_v1 = num_u_upper, num_v_lower, num_u_rib, num_v_rib 
    intersection_list_upper.append([surface111.pointset_id,pointset.pointset_id,num_points_u0,num_points_v0,num_points_u1,num_points_v1,'u1','T'])
shell_mesh.identify_intersection_list(geo, intersection_list_upper, plot =True)

# print()
# print('Plotting')
# vd_points2 = [] 
# for i, pointset in enumerate(OML_ctrl_pointset_list_wingbox_side_20):   #
#     geo.assemble(pointset = pointset)
#     points = geo.evaluate(pointset = pointset)
#     print(pointset.pointset_id,pointset.name, points.shape)
#     unique_points = np.unique(points, axis=0)
#     #print('unique',len(points),len(unique_points))    
#     color = list(vedo.colors.colors.values())[i]
#     if i%2==0:
#         vd_points2.append(vedo.Points(points, r=20, c=color, alpha=1.0)) 
#     else:
#         vd_points2.append(vedo.Points(points, r=15, c=color, alpha=1.0))
# vd_points3 = vedo.Points(points_side_curve0, r=25, c='red', alpha=0.3)
# vd_points4 = vedo.Points(points_side_curve1, r=25, c='red', alpha=0.3)
# #vd_points5 = vedo.Points(rib0_points, r=10, c='black', alpha=0.3)#surface_points0
# vd_points5 = vedo.Points(side_curve_points2, r=30, c='black', alpha=0.3)
# vd_points6 = vedo.Points(side_curve_points3, r=30, c='black', alpha=0.3)

# vd_points8 = vedo.Points(top_point2_vd, r=30, c='green', alpha=0.3)
# vd_points9 = vedo.Points(bot_point1_vd, r=30, c='blue', alpha=0.3)
# vd_points10 = vedo.Points(bot_point2_vd, r=30, c='blue', alpha=0.3),vd_points7 vd_points3,vd_points4,vd_points5,vd_points6,
#vd_points5 = vedo.Points(surface_points0, r=10, c='black', alpha=1.0)#surface_points0 ,vd_points7,vd_points8,vd_points9,vd_points10
vd_test = vedo.Plotter(axes=1)
vd_test.show(vd_points4,vd_points5, vd_points7, 'Test15', viewup="z", interactive=True)
exit()




print(surface111.pointset_id)



















