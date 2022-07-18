import numpy as np
import vedo
import triangle as tr
from meshopt import meshopt

points_temp0 = np.empty((0,3))
points_temp1 = np.empty((0,3))
for j in range(6):
    for i in range(10):
        points_temp0 = np.append(points_temp0, np.array([[i,10,j]]).reshape(1,3), axis = 0)
for j in range(6):
    for i in range(12):
        points_temp1 = np.append(points_temp1, np.array([[3.4,i-1,j]]).reshape(1,3), axis = 0)

vd_points_temp0 = vedo.Points(points_temp0, r=15, c='blue',alpha = 1)
A = dict(vertices=points_temp0[:,[0,2]])
B = tr.triangulate(A,'pc')   
connectivity0 = np.copy(B['triangles'])

mesh00 = vedo.Mesh([points_temp0,connectivity0], alpha=0.7,c='lightblue')
mesh00.backColor().lineColor('grey').lineWidth(0)  

vd_points_temp1 = vedo.Points(points_temp1, r=15, c='yellow',alpha = 1)
A = dict(vertices=points_temp1[:,[1,2]])
B = tr.triangulate(A,'pc')   
connectivity1 = np.copy(B['triangles'])
mesh10 = vedo.Mesh([points_temp1,connectivity1], alpha=1, c='lightyellow')
mesh10.backColor().lineColor('grey').lineWidth(0) 

vd_plotter = vedo.Plotter(shape=(2,2))
vd_plotter.show(vd_points_temp0,mesh00,vd_points_temp1,mesh10,'a',at=0, axes=0, viewup="z", interactive = False)      

indices1 = range(11,72,12)
vd_points_temp11 = vedo.Points(points_temp1[indices1,:], r=15, c='red',alpha = 1)

vd_plotter.show(vd_points_temp0,mesh00,vd_points_temp1,mesh10,vd_points_temp11,'b',at=1, axes=0, viewup="z", interactive = False)

indices_reduced = np.arange(0,60,10)+3
indices = []
for i in range(60):
    if i not in indices_reduced:
        indices.append(i)
vd_points_temp0 = vedo.Points(points_temp0[indices,:], r=15, c='blue',alpha = 1)
vd_plotter.show(vd_points_temp0,mesh00,vd_points_temp1,mesh10,vd_points_temp11,'c',at=2, axes=0, viewup="z", interactive = False)

# print('points_temp1[indices1,:]',points_temp1[indices1,:])
# print('points_temp0[indices_reduced,:]',points_temp0[indices_reduced,:])
points_temp0[indices_reduced,:] = points_temp1[indices1,:]
# print()
# print('points_temp1[indices1,:]',points_temp1[indices1,:])
# print('points_temp0[indices_reduced,:]',points_temp0[indices_reduced,:])
mesh01 = vedo.Mesh([points_temp0,connectivity0], alpha=0.9,c='lightblue')
mesh01.backColor().lineColor('black').lineWidth(5)
mesh11 = vedo.Mesh([points_temp1,connectivity1], alpha=1, c='lightyellow')
mesh11.backColor().lineColor('black').lineWidth(5)
vd_plotter.show(vd_points_temp0,mesh01,vd_points_temp1,mesh11,vd_points_temp11,'d',at=3, axes=0, viewup="z", interactive = True)