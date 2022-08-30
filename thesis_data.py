import vtk
import numpy as np
import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/pymeshopt')
import pymeshopt
import matplotlib.pyplot as plt
import vedo
file_path = '/Users/Sansara/Public/Code/Geomesh/ShellMesh/data/'
aspect_ratio = np.loadtxt(file_path + 'aspect_ratio.txt')#eVTOL_26486
internal_angles = np.loadtxt(file_path + 'internal_angles.txt')
print('aspect_ratio', aspect_ratio.shape,'internal_angles', internal_angles.shape)
# aspect_ratio = np.delete(aspect_ratio, np.argmax(aspect_ratio))
# aspect_ratio = np.delete(aspect_ratio, np.argmax(aspect_ratio))
print(len(aspect_ratio[aspect_ratio<1.5]))
print(len(aspect_ratio[aspect_ratio>=1.5]))
inter = internal_angles.flatten()
print(len(inter[inter>=110]))
num_bins = 15
bin_a = list(np.arange(1,3.8,0.35))
bin_i = list(np.arange(0,180,10))
fig, axs = plt.subplots(1,2, figsize= (13,4))
axs[0].hist(aspect_ratio, color = 'blue', alpha = 0.6, bins = bin_a, label = 'Our mesh', density = True)#+' - ' '%s'%len(aspect_ratio)+' elements'
axs[1].hist(internal_angles.flatten(), color = 'blue', alpha = 0.6, bins = bin_i, label = 'Our mesh', density = True)#+' - ' '%s'%len(aspect_ratio)+' elements'
vertexCoords = np.loadtxt(file_path + 'vertlist.txt')
quadlist = np.loadtxt(file_path + 'quadlist.txt').astype('int')
print('np.argmax(internal_angles)',np.argmax(internal_angles),np.max(internal_angles),int(np.argmax(internal_angles)/4))
print('vertexCoords',vertexCoords.shape,'quadlist',quadlist.shape)
#print(quadlist[np.argmax(internal_angles),:])
print(internal_angles[internal_angles>110])
vd_points = vedo.Points(vertexCoords[quadlist[int(np.argmax(internal_angles)/4),:],:], r=20, c='black')
mesh_quad = vedo.Mesh([vertexCoords, quadlist], alpha=0.9)
mesh_quad.backColor().lineColor('red').lineWidth(6) 
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,mesh_quad,'222', axes=1, viewup="z", interactive = False)#,vd_points1

path = 'ICEM_uCRM_mesh/uCRM_ICEM_medium.vtk'
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

axs[0].hist(aspect_ratio, color = 'red', alpha = 0.3, bins = bin_a, label = 'ICEMCFD mesh', density = True) #+' - ' '%s'%len(aspect_ratio)+' elements'
axs[1].hist(internal_angles.flatten(), color = 'red', alpha = 0.3, bins = bin_i, label = 'ICEMCFD mesh', density = True)#+' - ' '%s'%len(aspect_ratio)+' elements'

axs[0].set_xlabel('Aspect ratio')
axs[1].set_xlabel('Internal angles')
axs[0].set_ylabel('Density')
axs[1].set_ylabel('Density')
axs[1].set_xticks([0, 45, 90, 135, 180])

#axs[1].set_xlim(left=0, right=180)
fig.tight_layout()
axs[0].legend(loc = 'upper right')
axs[1].legend(loc = 'upper right')
plt.savefig(file_path + 'comparison')            
plt.show()
