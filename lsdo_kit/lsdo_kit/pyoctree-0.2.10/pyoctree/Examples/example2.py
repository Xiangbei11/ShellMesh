# Imports
from __future__ import print_function
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import sys, vtk
sys.path.append('../')
import pyoctree
from pyoctree import pyoctree as ot
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
import scipy as sp

print('pyoctree version = ', pyoctree.__version__)
print('vtk version = ', vtk.vtkVersion.GetVTKVersion())

# Read in stl file using vtk
reader = vtk.vtkSTLReader()
# reader.SetFileName("knot.stl")
# reader.SetFileName("aircraft_114k.stl")
reader.SetFileName("456.STL")
# reader.SetFileName("aircraft_7k.stl")
# reader.SetFileName("aircraft_7k.stl")
# reader.SetFileName("aircraft_7k.stl")
reader.MergingOn()
reader.Update()
stl = reader.GetOutput()
print("Number of points    = %d" % stl.GetNumberOfPoints())
print("Number of triangles = %d" % stl.GetNumberOfPolys())



decimate = vtk.vtkDecimatePro()
decimate.SetInputData(stl)
decimate.SetTargetReduction(.5)
decimate.Update()
decimate.SetMaximumError(0.0001)
decimatedPoly = vtk.vtkPolyData()
decimatedPoly.ShallowCopy(decimate.GetOutput())

print("After decimation \n"
        "-----------------\n")
print("Number of points    = %d" % decimatedPoly.GetNumberOfPoints())
print("Number of triangles = %d" % decimatedPoly.GetNumberOfPolys())






a=2


# Extract polygon info from stl
if a==1:
    stl=stl
else:
    stl=decimatedPoly
# 1. Get array of point coordinates
numPoints   = stl.GetNumberOfPoints()
pointCoords = np.zeros((numPoints,3),dtype=float)
for i in range(numPoints):
    pointCoords[i,:] = stl.GetPoint(i)
    
# 2. Get polygon connectivity
numPolys     = stl.GetNumberOfCells()
connectivity = np.zeros((numPolys,3),dtype=np.int32)
for i in range(numPolys):
    atri = stl.GetCell(i)
    ids = atri.GetPointIds()
    for j in range(3):
        connectivity[i,j] = ids.GetId(j)

        # Show format of pointCoords
pointCoords

# # Show format of connectivity
# connectivity

# # Create octree structure containing stl poly mesh
# tree = ot.PyOctree(pointCoords,connectivity)

# # Print out basic Octree data
# print("Size of Octree               = %.3fmm" % tree.root.size)
# print("Number of Octnodes in Octree = %d" % tree.getNumberOfNodes())
# print("Number of polys in Octree    = %d" % tree.numPolys)


# # Get the root node
# print(tree.root)

# # Get the branches of the root node, OctNode 0
# tree.root.branches

# # Get OctNode 0-0 (first branch of the root node)
# print(tree.root.branches[0])


# # Get the branches of OctNode 0-0 (first branch of the root node) using getNodeFromId function
# tree.getNodeFromId('0-0').branches

# # An octnode that has no branches is a leaf
# print('OctNode 0-0-1 is a leaf = ', tree.getNodeFromId('0-3-1').isLeaf)
# print(tree.getNodeFromId('0-3-1'))

# # Get list of first 10 polygons stored in tree
# tree.polyList[0:10]

# # Get details of the first tri in the tree
# tri = tree.polyList[0]
# s = []
# for i in range(3):
#     s.append('[%.3f, %.3f, %.3f]' % tuple(tri.vertices[i][:3]))
# s = ','.join(s)

# print('Tri label = %d' % tri.label)
# print('Vertex coordinates = %s' % s)
# print('Face normal direction = [%.2f, %.2f, %.2f]' % tuple(tri.N))
# print('Perp. distance from origin = %.3f' % tri.D)

# # Find the OctNode that the tri with label=0 lies within
# tree.getNodesFromLabel(0)
# # Test if OctNode contains a particular tri label
# tree.getNodeFromId('0-3-1').hasPolyLabel(0)

# # Create a list containing a single ray
# xs,xe,ys,ye,zs,ze = stl.GetBounds()
# x = 0.5*np.mean([xs,xe])
# y = np.mean([ys,ye])
# rayPointList = np.array([[[x,y,zs],[x,y,ze]],[[10,0,-5],[x,y,ze]]],dtype=np.float32)
# # # Find if an intersection occurred
# # for i in tree.rayIntersections(rayPointList):
# #     print(1)    
# #     print('Intersection coords = [%.2f, %.2f, %.2f]' % tuple(i.p), ',  Parametric. dist. along ray = %.2f' % i.s)

# #     # Get intersection points for a single ray




# print(111)
# ray = rayPointList[0]
# for i in tree.rayIntersection(ray):
#     print('Intersected tri = %d,' % i.triLabel, 'Intersection coords = [%.2f, %.2f, %.2f]' % tuple(i.p), ',  Parametric. dist. along ray = %.2f' % i.s)
#     # Create a vtk representation of Octree
# # projpoint = np.array([10,10,10],dtype=np.float32)
# # a = tree.findmindist(projpoint)
# # tree.getOctreeRep()
# # print(a.position,a.distance)

# a=i

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# ax.plot_trisurf(pointCoords[:,0],pointCoords[:,1],pointCoords[:,2], linewidth=0.2, antialiased=True)

ax.scatter(pointCoords[:,0],pointCoords[:,1],pointCoords[:,2],color='blue',s=0.1)
# ax.scatter(a.p[0],a.p[1],a.p[2],color='black',s=5)
# ax.scatter(x,y,zs,color='red',s=2)
# ax.scatter(p_proj[0],p_proj[1],p_proj[2],color='green',s=2)
X = pointCoords[:,0]
Y = pointCoords[:,1]
Z = pointCoords[:,2]
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')
plt.show()



# ax = a3.Axes3D(pl.figure())
# for i in range(len(connectivity)):
#     vtx = np.vstack((pointCoords[connectivity[i,0],:],pointCoords[connectivity[i,1],:],pointCoords[connectivity[i,2],:]))
#     # vtx=np.array([[1,1,1],[2,1,1],[1,0,3]])
#     tri = a3.art3d.Poly3DCollection([vtx])
#     # tri.set_color(colors.rgb2hex(sp.rand(3)))
#     tri.set_edgecolor('k')
#     ax.add_collection3d(tri)
# ax.set_xlim3d(-10, 40)
# ax.set_ylim3d(-20, 20)
# ax.set_zlim3d(-20, 20)
# pl.show()

