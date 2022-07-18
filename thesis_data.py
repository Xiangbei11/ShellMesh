import vtk
import numpy as np
import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/pymeshopt')
import pymeshopt
import matplotlib.pyplot as plt

path = 'uCRM_ICEM_mesh/uCRM_ICEM_coarse.vtk'
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(path)
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()
data = reader.GetOutput()
num_cells = data.GetNumberOfCells()
print('num_cells', num_cells)
num_points = data.GetNumberOfPoints()
print('num_points', num_points)
poly = 4
connectivity = np.zeros((num_cells,poly),dtype=np.int32)
for i in range(num_cells):
    atri = data.GetCell(i)
    ids = atri.GetPointIds()
    for j in range(ids.GetNumberOfIds()):
        connectivity[i,j] = ids.GetId(j)

pointCoords = np.zeros((num_points,3),dtype=float)
for i in range(num_points):
    pointCoords[i,:] = data.GetPoint(i)

print('connectivity',len(connectivity))
print('pointCoords',len(pointCoords))

trilist = np.empty((0,3),dtype=np.int32)
meshopt_test = pymeshopt.Pymeshopt(pointCoords.astype(np.float32),trilist.astype(np.int32),connectivity.astype(np.int32),1.,1.,1.) 
aspect_ratio = meshopt_test.calculate_aspect_ratio()
internal_angles = meshopt_test.calculate_internal_angles()
# np.savetxt('data/aspect_ratio.txt', aspect_ratio)
# np.savetxt('data/internal_angles.txt', internal_angles)
print('aspect_ratio', aspect_ratio.shape,'internal_angles', internal_angles.shape)
fig, axs = plt.subplots(2,1, figsize= (14,8))
axs[0].hist(aspect_ratio, color = 'blue', alpha = 0.5, bins = 50, label = 'ICEM')#, density = True)
axs[0].title.set_text('Aspect Ratio')
axs[1].hist(internal_angles.flatten(), color = 'red', alpha = 0.5, bins =20, label = 'ICEM')#, density = True)
axs[1].title.set_text('Internal Angles')
axs[1].set_xticks([0, 45, 90, 135, 180])
aspect_ratio = np.loadtxt('/Users/Sansara/Public/Code/Geomesh/ShellMesh/data/aspect_ratio.txt')#eVTOL_26486
internal_angles = np.loadtxt('/Users/Sansara/Public/Code/Geomesh/ShellMesh/data/internal_angles.txt')
print('aspect_ratio', aspect_ratio.shape,'internal_angles', internal_angles.shape)
axs[0].hist(aspect_ratio, color = 'blue', alpha = 0.8, bins = 50, label = 'ShellMesh')#, density = True)
axs[1].hist(internal_angles.flatten(), color = 'red', alpha = 0.8, bins =20, label = 'ShellMesh')#, density = True)
#axs[1].set_xlim(left=0, right=180)
fig.tight_layout()
axs[0].legend(loc = 'upper right')
axs[1].legend(loc = 'upper right')
# plt.savefig('data/histogram')            
plt.show()
