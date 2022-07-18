import numpy as np
import vedo
import triangle as tr
from meshopt import meshopt

points_temp0 = np.empty((0,3))
points_temp1 = np.empty((0,3))
for j in range(2):
    for i in range(3):
        points_temp0 = np.append(points_temp0, np.array([[i,10,j]]).reshape(1,3), axis = 0)

A = dict(vertices=points_temp0[:,[0,2]])
B = tr.triangulate(A,'pc')   
connectivity0 = np.copy(B['triangles'])

quadlist = np.empty((0,4),dtype=np.int32)
vertexCoords = points_temp0.astype('float32')
trilist = np.copy(connectivity0).astype('int32')
fixedvert = np.array([0,1,2], dtype='int32')
print('trilist',trilist.shape)
print('vertexCoords',vertexCoords.shape)
itr = [2]
m = meshopt(vertexCoords,trilist,quadlist, fixedvert=fixedvert ,itr=itr, w1=1.,w2=1.,w3=1., plot = 0)#,w4 =1.
m.optimization()
