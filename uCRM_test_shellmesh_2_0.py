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
import triangle as tr
import matplotlib.pyplot as plt
from member import Member

shell_mesh = ShellMesh('shell_mesh')
bspline_surface_list = list(geo.input_bspline_entity_dict.values())[:160]
OML_pointset_list = shell_mesh.extract_pointset_list_from_bspline_surface(geo, bspline_surface_list)
for i, pointset in enumerate(OML_pointset_list):
    print(pointset.pointset_id,pointset.name)
print()
scale = 100
vd0 = []
surface_upper0 = np.empty((0,3))
surface_lower0 = np.empty((0,3))
k = 0
for i in range(111,157):
    if i == 125 or i ==123:
        pass
    else:
        pointset = OML_pointset_list[i]
        print('rib', i, pointset.pointset_id,pointset.name)
        pointset.shape= np.array([11, 11, 3])
        num_points_u0, num_points_v0  = 5, 2       
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
        basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
        geo.assemble(pointset = pointset)
        geo.evaluate(pointset = pointset)    
        test_points = basis_matrix.dot(pointset.physical_coordinates)
        color = list(vedo.colors.colors.values())[i]
        vd0.append(vedo.Points(test_points, r=8, c = color)) 

        indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
        surface_upper0 = np.append(surface_upper0, test_points[indices_u1,:],axis=0)
        indices_u0 = num_points_v0 * np.arange(num_points_u0)
        surface_lower0 = np.append(surface_lower0, test_points[indices_u0,:],axis=0) 

        if i == 111:
            indices = [0,1,num_points_u0*num_points_v0-2,num_points_u0*num_points_v0-1]
            print(indices)
            A = dict(vertices=u0_v0_vec[indices,:])#, segments=constrained_edges
            B = tr.triangulate(A,'pc') 
            # tr.compare(plt, A, B)         
            # plt.show()#block = False 
            connectivity = np.copy(B['triangles'])
            shell_mesh.members_dict['rib'+'%s'%i] = Member( 
                id = k,
                name = 'rib'+'%s'%k,   
                coordinates = test_points[indices,:]/scale,
                tri_connectivity = connectivity,
                constrained_node_indices = [],
                )    
            k+=1            
        else:
            A = dict(vertices=u0_v0_vec)#, segments=constrained_edges
            B = tr.triangulate(A,'pc') 
            connectivity = np.copy(B['triangles'])
            shell_mesh.members_dict['rib'+'%s'%i] = Member( 
                id = k,   
                name = 'rib'+'%s'%k,
                coordinates = test_points/scale,
                tri_connectivity = connectivity,
                constrained_node_indices = [],
                )      
            k+=1
vd1 = vedo.Points(surface_upper0, r=15, c = 'black', alpha = 0.5)

vd2 = []
surface_upper1 = np.empty((0,3))
surface_lower1 = np.empty((0,3))
for i in range(103,108):
    pointset = OML_pointset_list[i]
    print('rib', i, pointset.pointset_id,pointset.name)
    pointset.shape= np.array([11, 11, 3])
    if i == 107:
        num_points_u0, num_points_v0  = 2, 2 
    else:
        num_points_u0, num_points_v0  = 5, 2      
    u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
    v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
    u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
    basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
    geo.assemble(pointset = pointset)
    geo.evaluate(pointset = pointset)    
    test_points = basis_matrix.dot(pointset.physical_coordinates)
    color = list(vedo.colors.colors.values())[i]
    vd2.append(vedo.Points(test_points, r=10, c = color)) 

    indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
    surface_upper1 = np.append(surface_upper1, test_points[indices_u1,:],axis=0)
    indices_u0 = num_points_v0 * np.arange(num_points_u0)
    surface_lower1 = np.append(surface_lower1, test_points[indices_u0,:],axis=0)
    if i == 107:
        pass
    else:
        A = dict(vertices=u0_v0_vec)#, segments=constrained_edges
        B = tr.triangulate(A,'pc') 
        connectivity = np.copy(B['triangles'])
        shell_mesh.members_dict['rib'+'%s'%i] = Member(  
            id = k,  
            name = 'rib'+'%s'%k,
            coordinates = test_points/scale,
            tri_connectivity = connectivity,
            constrained_node_indices = [],
            )  
        k+=1
vd3 = vedo.Points(surface_upper1, r=20, c = 'black', alpha = 0.5)#

total_surface_upper = np.append(surface_upper0[0,:].reshape((1,3)),surface_upper0[5:,:],axis=0)
total_surface_upper = np.append(surface_upper1, total_surface_upper,axis=0) 
print('total_surface_upper',total_surface_upper.shape)
vd4 = vedo.Points(total_surface_upper, r=15, c = 'blue', alpha = 0.5)

total_surface_lower = np.append(surface_lower0[0,:].reshape((1,3)),surface_lower0[5:,:],axis=0)
total_surface_lower = np.append(surface_lower1, total_surface_lower,axis=0) 
print('total_surface_lower',total_surface_lower.shape)
vd5 = vedo.Points(total_surface_lower, r=15, c = 'red', alpha = 0.5)

constrained_edges = np.empty((0,2),dtype = int)
constrained_edges = np.append(constrained_edges, \
    np.array([[0,1],[1,2],[2,3],[3,4],[3,4],[4,9],[9,14],[14,19],[19,37]\
        ]), axis = 0)
num_points_uu0 = 45
num_points_vv0 = 5
in_start = 37
for i in np.arange(in_start,in_start+(num_points_uu0-5)*num_points_vv0,num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i+5]).reshape(1,2), axis = 0)  
for i in np.arange(in_start+(num_points_uu0-5)*num_points_vv0,in_start+(num_points_uu0-5)*num_points_vv0-num_points_vv0+1,-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
for i in np.arange(in_start+(num_points_uu0-5)*num_points_vv0-num_points_vv0+1,23,-5):
    constrained_edges = np.append(constrained_edges, np.array([i, i-5]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[23,22],[22,20],[20,15],[15,10],[10,5],[5,0]]), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[20,21],[21,27],[27,32]]), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[23,24],[24,25],[25,26],[26,27]]), axis = 0)
#print('constrained_edges',constrained_edges)
A = dict(vertices=total_surface_upper[:,:2], segments=constrained_edges)#
B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show(block = False)
connectivity = np.copy(B['triangles'])
print('total_surface_upper',len(total_surface_upper),'connectivity',len(connectivity))
shell_mesh.members_dict['upper_surface'] = Member( 
    id = k,
    name = 'upper_surface',
    coordinates = total_surface_upper/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
k+=1
shell_mesh.members_dict['lower_surface'] = Member( 
    id = k,   
    name = 'lower_surface',
    coordinates = total_surface_lower/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
k+=1
indices = [20,21,27,32,37] 
points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
ids = [0,2]
A = dict(vertices=points[:,ids])
B = tr.triangulate(A,'pc') 
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['rib'+'%s'%i] = Member(   
    id = k, 
    name = 'rib'+'%s'%(k-2),
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
print('k-2',k-2)
k+=1
indices = [4,9,14,19] + list(np.arange(in_start,in_start+(num_points_uu0-5)*num_points_vv0+5,num_points_vv0))
#print('indices',indices)
points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
#np.append(constrained_edges, np.array([[0,1],[1,2],[2,3],[3,4]]).reshape(4,2), axis = 0)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(2*len(indices)-1,len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
#print(constrained_edges)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
connectivity = np.copy(B['triangles'])
# tr.compare(plt, A, B)         
# plt.show(block=False)#block=False
shell_mesh.members_dict['side0'] = Member(  
    id = k,  
    name = 'side0',
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
k+=1
indices = [0,5,10,15,20,22] + list(np.arange(23,23+num_points_uu0*num_points_vv0-2*num_points_vv0,num_points_vv0))
#print('indices',indices)
vd6 = vedo.Points(total_surface_upper[indices,:], r=25, c = 'black', alpha = 0.5)
vd7 = vedo.Points(total_surface_lower[indices,:], r=25, c = 'black', alpha = 0.5)
points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
#np.append(constrained_edges, np.array([[0,1],[1,2],[2,3],[3,4]]).reshape(4,2), axis = 0)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(2*len(indices)-1,len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
#print(constrained_edges)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
connectivity = np.copy(B['triangles'])
# tr.compare(plt, A, B)         
# plt.show()#block=False
shell_mesh.members_dict['side1'] = Member(
    id = k,     
    name = 'side1',
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  


shell_mesh.optimizie_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_uCRM_shellmesh_2')

# mesh0 = vedo.Mesh([total_surface_upper, connectivity], alpha=0.5)
# mesh0.backColor().lineColor('green').lineWidth(3) 
# mesh1 = vedo.Mesh([total_surface_lower, connectivity], alpha=0.5)
# mesh1.backColor().lineColor('green').lineWidth(3) 
# vd_plotter = vedo.Plotter()
# vd_plotter.show(mesh0,vd4,mesh1,vd5,'111', axes=1, viewup="z", interactive = False)

vd_mesh = []
vd_points = []
count = 0
for memb in shell_mesh.members_dict.values():
    print(memb.options['name'])
    conn = np.copy(memb.options['tri_connectivity']) 
    count += len(conn)
    pts = np.copy(memb.options['coordinates'])
    mesh = vedo.Mesh([pts,conn], alpha=0.5)
    mesh.backColor().lineColor('green').lineWidth(3)
    vd_mesh.append(mesh)
    if memb.options['name'] == 'upper_surface':
        vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    else:
        vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
print('count',count)
vd_plotter1 = vedo.Plotter()#vd6,vd7,
vd_plotter1.show(vd_points,vd_mesh,'222', axes=1, viewup="z", interactive = True)
exit()

