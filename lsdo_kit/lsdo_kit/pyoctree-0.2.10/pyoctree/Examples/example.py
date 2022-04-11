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


print('pyoctree version = ', pyoctree.__version__)
print('vtk version = ', vtk.vtkVersion.GetVTKVersion())

# Read in stl file using vtk
reader = vtk.vtkSTLReader()
# reader.SetFileName("knot.stl")
# reader.SetFileName("aircraft_114k.stl")
reader.SetFileName("aircraft_7k.stl")
# reader.SetFileName("aircraft_7k.stl")
# reader.SetFileName("aircraft_7k.stl")
# reader.SetFileName("aircraft_7k.stl")
reader.MergingOn()
reader.Update()
stl = reader.GetOutput()
stl = reader.GetOutput()
print("Number of points    = %d" % stl.GetNumberOfPoints())
print("Number of triangles = %d" % stl.GetNumberOfCells())

# Extract polygon info from stl

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

# Show format of connectivity
connectivity

# Create octree structure containing stl poly mesh
tree = ot.PyOctree(pointCoords,connectivity)

# # Print out basic Octree data
# print("Size of Octree               = %.3fmm" % tree.root.size)
# print("Number of Octnodes in Octree = %d" % tree.getNumberOfNodes())
# print("Number of polys in Octree    = %d" % tree.numPolys)

tree.calculateavgedgesize()
# Get the root node
# print(tree.root)

# Get the branches of the root node, OctNode 0
tree.root.branches

# Get OctNode 0-0 (first branch of the root node)
# print(tree.root.branches[0])


# Get the branches of OctNode 0-0 (first branch of the root node) using getNodeFromId function
tree.getNodeFromId('0-0').branches

# An octnode that has no branches is a leaf
# print('OctNode 0-0-1 is a leaf = ', tree.getNodeFromId('0-3-1').isLeaf)
# print(tree.getNodeFromId('0-3-1'))

# Get list of first 10 polygons stored in tree
tree.polyList[0:10]

# Get details of the first tri in the tree
tri = tree.polyList[0]
s = []
for i in range(3):
    s.append('[%.3f, %.3f, %.3f]' % tuple(tri.vertices[i][:3]))
s = ','.join(s)

print('Tri label = %d' % tri.label)
print('Vertex coordinates = %s' % s)
print('Face normal direction = [%.2f, %.2f, %.2f]' % tuple(tri.N))
print('Perp. distance from origin = %.3f' % tri.D)

# Find the OctNode that the tri with label=0 lies within
tree.getNodesFromLabel(0)
# Test if OctNode contains a particular tri label
tree.getNodeFromId('0-3-1').hasPolyLabel(0)

# Create a list containing a single ray
xs,xe,ys,ye,zs,ze = stl.GetBounds()
x = 0.5*np.mean([xs,xe])
y = np.mean([ys,ye])
# ---------------------------------------------------------------------------------------------------
p_starting = np.array([18,10,1.5],dtype=np.float32)# x and y coordinates for projection
p_ending = np.array([20,10,1.5],dtype=np.float32)
p_proj = np.array([20,4,-2])
p_vec1 = p_ending-p_starting
p_vec2 = p_ending-p_proj
p_normal = np.cross(p_vec1,p_vec2)/np.linalg.norm(np.cross(p_vec1,p_vec2))

# num_int = int(round(np.linalg.norm(p_starting-p_ending)/e_size))
num_int=100
Pt = np.zeros([num_int,3])
Pt[:,0] = np.linspace(p_starting[0],p_ending[0],num=num_int,endpoint=True,dtype=np.float32)#.reshape(num_int,1)
Pt[:,1] = np.linspace(p_starting[1],p_ending[1],num=num_int,endpoint=True,dtype=np.float32)#.reshape(num_sint,1)
Pt[:,2] = np.linspace(p_starting[2],p_ending[2],num=num_int,endpoint=True,dtype=np.float32)#.reshape(num_int,1)
Pt=np.array(Pt,dtype=np.float32)
Pproj=np.zeros((len(Pt),3),dtype=np.float32)
Pproj[:,0]=Pt[:,0]
Pproj[:,1]=Pt[:,1]-2
Pproj[:,2]=Pt[:,2]-4
rayList=np.zeros((len(Pt),2,3))
for i in range(len(Pt)):
    rayList[i,:,:]=np.array([Pt[i,:],Pproj[i,:]])
    # rayList[2*i+1,:]=Pproj[i,:]
rayList=np.array(rayList,dtype=np.float32)
# ---------------------------------------------------------------------------------------------------
rayPointList = np.array([[[x,y,zs],[x,y,ze]],[[10,0,5],[13,0,4]]],dtype=np.float32)
rayPointList = np.array([[[30,5,0],[30,5,-5]]],dtype=np.float32)
# Find if an intersection occurred
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

print('------------------------------------------------------------------------------------------')
'''
j=0
for i in tree.rayIntersections(rayList):
    print('Intersection coords = [%.2f, %.2f, %.2f]' % tuple(i.p), ',  Parametric. dist. along ray = %.2f' % i.s)
    ax.scatter(i.p[0],i.p[1],i.p[2],color='black',s=5)
    ax.scatter(rayList[j,0,0],rayList[j,0,1],rayList[j,0,2],color='red',s=2)
    ax.quiver(rayList[j,0,0],rayList[j,0,1],rayList[j,0,2],rayList[j,1,0]-rayList[j,0,0],rayList[j,1,1]-rayList[j,0,1],rayList[j,1,2]-rayList[j,0,2],length=0.6, normalize=True)
    # ax.scatter(rayList[j,0,0],rayList[j,0,1],rayList[j,0,2],color='red',s=2)

    j=j+1
    # Get intersection points for a single ray

# ax.scatter(a.position[0],a.position[1],a.position[2],color='black',s=5)



ax.scatter(pointCoords[:,0],pointCoords[:,1],pointCoords[:,2],color='blue',s=0.1)
# ax.scatter(r[0,0,0],rayPointList[0,0,1],rayPointList[0,0,2],color='red',s=2)
# ax.scatter(18,7,10,color='red',s=2)
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

'''
j=0
a=tree.findmindists(Pt)
for i in tree.findmindists(Pt):
    ax.scatter(i.position[0],i.position[1],i.position[2],color='black',s=5)
    # ax.scatter(rayList[j,0,0],rayList[j,0,1],rayList[j,0,2],color='red',s=2)
    # ax.quiver(rayList[j,0,0],rayList[j,0,1],rayList[j,0,2],rayList[j,1,0]-rayList[j,0,0],rayList[j,1,1]-rayList[j,0,1],rayList[j,1,2]-rayList[j,0,2],length=0.6, normalize=True)
    ax.scatter(Pt[j,0],Pt[j,1],Pt[j,2],color='red',s=2)
    j=j+1




ax.scatter(pointCoords[:,0],pointCoords[:,1],pointCoords[:,2],color='blue',s=0.1)
# ax.scatter(r[0,0,0],rayPointList[0,0,1],rayPointList[0,0,2],color='red',s=2)
# ax.scatter(18,7,10,color='red',s=2)
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














# ray = rayPointList[0]
# for i in tree.rayIntersection(ray):
#     print('Intersected tri = %d,' % i.triLabel, 'Intersection coords = [%.2f, %.2f, %.2f]' % tuple(i.p), ',  Parametric. dist. along ray = %.2f' % i.s)
    # Create a vtk representation of Octree
# projpoint = np.array([18,7,10],dtype=np.float32)
# a = tree.findmindist(projpoint)
# tree.getOctreeRep()
# print(a.position,a.distance)

# a=i
# # ax.scatter(a.position[0],a.position[1],a.position[2],color='black',s=5)
# ax.scatter(a.p[0],a.p[1],a.p[2],color='black',s=5)



# ax.scatter(pointCoords[:,0],pointCoords[:,1],pointCoords[:,2],color='blue',s=0.1)
# ax.scatter(rayPointList[0,0,0],rayPointList[0,0,1],rayPointList[0,0,2],color='red',s=2)
# # ax.scatter(18,7,10,color='red',s=2)
# # ax.scatter(p_proj[0],p_proj[1],p_proj[2],color='green',s=2)
# X = pointCoords[:,0]
# Y = pointCoords[:,1]
# Z = pointCoords[:,2]
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# # Create cubic bounding box to simulate equal aspect ratio
# max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
# Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
# Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
# Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
# # Comment or uncomment following both lines to test the fake bounding box:
# for xb, yb, zb in zip(Xb, Yb, Zb):
#    ax.plot([xb], [yb], [zb], 'w')
# # plt.show()