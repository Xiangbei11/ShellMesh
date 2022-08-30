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
scale = 1000
k = 0
vd0 = []
surface_upper0 = np.empty((0,3))
surface_lower0 = np.empty((0,3))
side0_indices = []
side1_indices = []
num_points = 0
constrained_node_indices = []
constrained_node_indices_lower = []
num_u0_root= 7
for i in range(103,108):
    pointset = OML_pointset_list[i]
    print('rib', i, pointset.pointset_id,pointset.name)
    pointset.shape= np.array([11, 11, 3])
    if i == 107:
        num_points_u0, num_points_v0  = 2, 3 
    # elif i == 106:
    #     num_points_u0, num_points_v0  = 6, 2    
    # elif i == 103:
    #     num_points_u0, num_points_v0  = 7, 2         
    else:
        num_points_u0, num_points_v0  = num_u0_root, 3   

    side0_indices.append(num_points)   
    side1_indices.append(num_points+num_points_u0-1) 
    constrained_node_indices += list(np.arange(num_points,num_points+num_points_u0))
    num_points += num_points_u0   

    u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
    v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
    u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
    basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
    geo.assemble(pointset = pointset)
    geo.evaluate(pointset = pointset)    
    test_points = basis_matrix.dot(pointset.physical_coordinates)
    color = list(vedo.colors.colors.values())[i]
    vd0.append(vedo.Points(test_points/scale, r=10, c = color)) 

    indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
    surface_upper0 = np.append(surface_upper0, test_points[indices_u1,:],axis=0)
    indices_u0 = num_points_v0 * np.arange(num_points_u0)
    surface_lower0 = np.append(surface_lower0, test_points[indices_u0,:],axis=0)
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
            constrained_node_indices = list(indices_u1) + list(indices_u0),
            )  
        k+=1
vd1 = vedo.Points(surface_upper0/scale, r=15, c = 'black', alpha = 0.5)#
side1_indices = side1_indices
# print('constrained_node_indices',constrained_node_indices)
# print('num_points',num_points)

total_surface_upper = np.append(surface_upper0[:-1,:], np.array([32360.04710327, 3100., 4388.06828032]).reshape((1,3)),axis=0)
total_surface_lower = np.append(surface_lower0[:-1,:], np.array([32360.05802888, 3100., 3155.82527481]).reshape((1,3)),axis=0)

# indices = side1_indices
# points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
# vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 1.0)
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd0,vd1,vd_points1,'000', axes=1, viewup="z", interactive = True)
# exit()
rib_48_index = num_points-1 
vd2 = []
surface_upper1 = np.empty((0,3))
surface_lower1 = np.empty((0,3))
rib48_indices = [num_points-2]
rib49_indices = []
num_v = 3
for i in range(111,157):
    if i == 125 or i ==123:
        pass
    else:
        pointset = OML_pointset_list[i]
        print('rib', i, pointset.pointset_id,pointset.name)
        pointset.shape= np.array([11, 11, 3])
        if i == 111:
            num_points_u0, num_points_v0  = 2, num_v
        elif i==112:             
            num_points_u0, num_points_v0  = 3, num_v
        elif i == 113:
            num_points_u0, num_points_v0  = 5, num_v
        elif i == 114:
            num_points_u0, num_points_v0  = 6, num_v
        elif i >= 154:
            num_points_u0, num_points_v0  = 2, num_v#3, 2
        elif i >= 150 and i <= 156:
            num_points_u0, num_points_v0  = 3, num_v
        elif i >= 140 and i <= 149:
            num_points_u0, num_points_v0  = 4, num_v
        elif i >= 132 and i <= 139:
            num_points_u0, num_points_v0  = 5, num_v   
        elif i >= 123 and i <= 131:
            num_points_u0, num_points_v0  = 6, num_v 
        elif i >= 118 and i <= 122:
            num_points_u0, num_points_v0  = 7, num_v                          
        else:        
            num_points_u0, num_points_v0  = 8, num_v

        side0_indices.append(num_points)   
        if i >= 114:
            side1_indices.append(num_points+num_points_u0-1) 
        else:
            rib48_indices.append(num_points+num_points_u0-1)
            print('rib49_indices',)
        if i == 114:
            rib48_indices.append(rib_48_index)
        constrained_node_indices += list(np.arange(num_points,num_points+num_points_u0))
        num_points += num_points_u0
        
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
        basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
        geo.assemble(pointset = pointset)
        geo.evaluate(pointset = pointset)    
        test_points = basis_matrix.dot(pointset.physical_coordinates)
        if i<=114:
            test_points[-1,1] = 3.1*scale
            test_points[-2,1] = 3.1*scale
            test_points[-3,1] = 3.1*scale
        if i==124:
            # print('i==124',k)
            # print(len(test_points),test_points)
            rib49_indices.append(num_points-num_points_u0+1)
            test_points[3,:] = np.array([32980,10458,4055])
            test_points[4,:] = np.array([32980,10458,4432.5])
            test_points[5,:] = np.array([32980,10458,4810])         
        if i==126:
            # print()
            # print('i==126',k)
            # print(len(test_points),test_points)
            rib49_indices.append(num_points-num_points_u0+1)
            test_points[3,:] = np.array([32981,11303,4144])
            test_points[4,:] = np.array([32981,11303,4465])
            test_points[5,:] = np.array([32981,11303,4786]) 
        if i == 114:
            print('i==115')
            print(len(test_points),test_points)  
            test_points[15,:] = np.array([32360.0390864,3185.8868092,3155.9209041])
            test_points[16,:] = np.array([32360.0390864,3185.8868092,3771.94291074])            
            test_points[17,:] = np.array([32360.0390864,3185.8868092,4387.96491738])                 
        
        color = list(vedo.colors.colors.values())[i]
        vd2.append(vedo.Points(test_points/scale, r=8, c = color)) 

        indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
        surface_upper1 = np.append(surface_upper1, test_points[indices_u1,:],axis=0)
        indices_u0 = num_points_v0 * np.arange(num_points_u0)
        surface_lower1 = np.append(surface_lower1, test_points[indices_u0,:],axis=0) 

        A = dict(vertices=u0_v0_vec)#, segments=constrained_edges
        B = tr.triangulate(A,'pc') 
        connectivity = np.copy(B['triangles'])
        shell_mesh.members_dict['rib'+'%s'%i] = Member( 
            id = k,   
            name = 'rib'+'%s'%k,
            coordinates = test_points/scale,
            tri_connectivity = connectivity,
            constrained_node_indices = list(indices_u1) + list(indices_u0),
            )      
        k+=1
vd3 = vedo.Points(surface_upper1/scale, r=15, c = 'black', alpha = 0.5)
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd0,vd2,vd3,'000', axes=1, viewup="z", interactive = True)
# exit()

total_surface_upper = np.append(total_surface_upper, surface_upper1[:-2,:],axis=0)
total_surface_lower = np.append(total_surface_lower, surface_lower1[:-2,:],axis=0)
# vd4 = vedo.Points(total_surface_upper/scale, r=15, c = 'blue', alpha = 0.7)
# vd5 = vedo.Points(total_surface_lower/scale, r=15, c = 'red', alpha = 0.7)
# vd_plotter0 = vedo.Plotter()
# vd_plotter0.show(vd4,vd5,'000', axes=1, viewup="z", interactive = True)#
# exit()
print('total_surface_upper',total_surface_upper.shape, 'num_points',num_points)
print('side0',side0_indices)
print('side1',side1_indices)
total_surface_upper = np.append(total_surface_upper, ((surface_upper1[-4,:]+surface_upper1[-2,:])/2).reshape((1,3)),axis=0)
total_surface_upper = np.append(total_surface_upper, ((surface_upper1[-3,:]+surface_upper1[-1,:])/2).reshape((1,3)),axis=0)
total_surface_upper = np.append(total_surface_upper, surface_upper1[-2:,:],axis=0)
total_surface_upper[rib48_indices,1] = 3.1*scale
total_surface_upper = np.append(total_surface_upper, ((total_surface_upper[rib48_indices[1],:]+total_surface_upper[rib48_indices[2],:])/2).reshape((1,3)),axis=0)
total_surface_upper = np.append(total_surface_upper, ((total_surface_upper[rib48_indices[2],:]+total_surface_upper[rib48_indices[3],:])/2).reshape((1,3)),axis=0)
# total_surface_upper = np.append(total_surface_upper, ((total_surface_upper[rib48_indices[0],:]+total_surface_upper[rib48_indices[1],:])/2).reshape((1,3)),axis=0)##1795
total_surface_lower = np.append(total_surface_lower, ((surface_lower1[-4,:]+surface_lower1[-2,:])/2).reshape((1,3)),axis=0)
total_surface_lower = np.append(total_surface_lower, ((surface_lower1[-3,:]+surface_lower1[-1,:])/2).reshape((1,3)),axis=0)
total_surface_lower = np.append(total_surface_lower, surface_lower1[-2:,:],axis=0)
total_surface_lower[rib48_indices,1] = 3.1*scale
# print(side0_indices)
# print()
# print(side1_indices)
# vd4 = vedo.Points(total_surface_upper/scale, r=15, c = 'blue', alpha = 0.7)
# vd5 = vedo.Points(total_surface_lower/scale, r=15, c = 'red', alpha = 0.7)
# vd_plotter0 = vedo.Plotter()
# vd_plotter0.show(vd4,vd5,'000', axes=1, viewup="z", interactive = True)#
# exit()
total_surface_lower = np.append(total_surface_lower, ((total_surface_lower[rib48_indices[1],:]+total_surface_lower[rib48_indices[2],:])/2).reshape((1,3)),axis=0)
total_surface_lower = np.append(total_surface_lower, ((total_surface_lower[rib48_indices[2],:]+total_surface_lower[rib48_indices[3],:])/2).reshape((1,3)),axis=0)
print('total_surface_upper',total_surface_upper.shape,'k',k)
rib48_indices_updated = rib48_indices[0:2] + [num_points+num_points_u0] + [rib48_indices[2]] + [num_points+num_points_u0+1] + rib48_indices[3:]
# rib48_indices_updated = [rib48_indices[0]]+ [num_points+num_points_u0+2] + [rib48_indices[1]] + [num_points+num_points_u0] + [rib48_indices[2]] + [num_points+num_points_u0+1] + rib48_indices[3:]##
print('rib48_indices',rib48_indices)
print('rib48_indices_updated',rib48_indices_updated)
side0_indices.append(num_points)   
side1_indices.append(num_points+num_points_u0-1) 
constrained_node_indices += list(np.arange(num_points,num_points+num_points_u0))
# print('side0',side0_indices)
# print('side1',side1_indices)

# vd4 = vedo.Points(total_surface_upper/scale, r=15, c = 'blue', alpha = 0.7)
# vd5 = vedo.Points(total_surface_lower/scale, r=15, c = 'red', alpha = 0.7)
# vd00 = vedo.Points(total_surface_upper[side0_indices,:]/scale, r=20, c = 'green', alpha = 0.5)
# vd01 = vedo.Points(total_surface_upper[side1_indices,:]/scale, r=20, c = 'black', alpha = 0.5)
# vd48 = vedo.Points((total_surface_upper[rib48_indices_updated,:]/scale), r=30, c = 'yellow', alpha = 0.3)
# vd49 = vedo.Points((total_surface_upper[rib49_indices,:]/scale), r=30, c = 'brown', alpha = 0.3)
# vd_plotter0 = vedo.Plotter()
# vd_plotter0.show(vd4,vd5,vd00,vd01,vd48,vd49,'000', axes=1, viewup="z", interactive = True)#
# exit()

constrained_edges = np.empty((0,2),dtype = int)
for i in np.arange(0,num_u0_root-1,1):
    constrained_edges = np.append(constrained_edges, np.array([i,i+1]).reshape(1,2), axis = 0)
for i in range(len(side1_indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([side1_indices[i], side1_indices[i+1]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[side1_indices[i+1],side0_indices[-1]]]).reshape(1,2), axis = 0) 
for i in np.arange(len(side0_indices)-1,0,-1):
    constrained_edges = np.append(constrained_edges, np.array([side0_indices[i], side0_indices[i-1]]).reshape(1,2), axis = 0)
for i in range(len(rib48_indices_updated)-1):
    constrained_edges = np.append(constrained_edges, np.array([rib48_indices_updated[i], rib48_indices_updated[i+1]]).reshape(1,2), axis = 0)
indices = [32,33,34]
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([indices[i], indices[i+1]]).reshape(1,2), axis = 0)
A = dict(vertices=total_surface_upper[:,:2], segments=constrained_edges)
B = tr.triangulate(A,'p') 
tr.compare(plt, A, B) 
# plt.show()
# exit()
plt.show(block = False)

connectivity = np.copy(B['triangles'])
mesh0 = vedo.Mesh([total_surface_upper/scale, connectivity], alpha=0.5)
mesh0.backColor().lineColor('green').lineWidth(3) 

print('total_surface_upper',len(total_surface_upper),'connectivity',len(connectivity))

constrained_node_indices += side0_indices
constrained_node_indices += side1_indices
constrained_node_indices += rib48_indices_updated #28, 31, 247, 34, 248, 39, 29
constrained_node_indices += rib49_indices
constrained_node_indices += [33,247]
constrained_node_indices += [50,43]
shell_mesh.members_dict['upper_surface'] = Member( 
    id = k,
    name = 'upper_surface',
    coordinates = total_surface_upper/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = constrained_node_indices,
    )  
k+=1

shell_mesh.members_dict['lower_surface'] = Member( 
    id = k,   
    name = 'lower_surface',
    coordinates = total_surface_lower/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = constrained_node_indices,
    )  
k+=1

indices = rib48_indices_updated
points = np.append(total_surface_upper[indices,:],(total_surface_upper[indices,:]+total_surface_lower[indices,:])/2, axis = 0)
points = np.append(points,total_surface_lower[indices,:], axis = 0)
ids = [0,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices)-1,3*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(3*len(indices)-1,2*len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices),len(indices)]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show()
# exit()
connectivity = np.copy(B['triangles'])
print('TTTT',len(connectivity), connectivity)
shell_mesh.members_dict['rib'+'%s'%i] = Member(   
    id = k, 
    name = 'rib'+'%s'%(k-2),
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = list(np.arange(len(rib48_indices_updated))) + list(np.arange(len(rib48_indices_updated),2*len(rib48_indices_updated))),
    )  
print('rib',k-2)
k+=1
i+=1

indices = rib49_indices
points = np.append(total_surface_upper[indices,:],(total_surface_upper[indices,:]+total_surface_lower[indices,:])/2, axis = 0)
points = np.append(points,total_surface_lower[indices,:], axis = 0)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices)-1,3*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(3*len(indices)-1,2*len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices),len(indices)]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['rib'+'%s'%i] = Member(   
    id = k, 
    name = 'rib'+'%s'%(k-2),
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [0,1,2,3],
    )  
print('rib',k-2)
k+=1

# vd_plotter = vedo.Plotter()
# vd_plotter.show(vd4,vd5,mesh0,'111', axes=1, viewup="z", interactive = False)

indices = side0_indices
points = np.append(total_surface_upper[indices,:],(total_surface_upper[indices,:]+total_surface_lower[indices,:])/2, axis = 0)
points = np.append(points,total_surface_lower[indices,:], axis = 0)
vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 1.0)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices)-1,3*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(3*len(indices)-1,2*len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices),len(indices)]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['side0'] = Member(  
    id = k,  
    name = 'side0',
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = list(np.arange(len(side0_indices))) + list(np.arange(len(side0_indices),2*len(side0_indices))),
    )  
k+=1

indices = side1_indices
points = np.append(total_surface_upper[indices,:],(total_surface_upper[indices,:]+total_surface_lower[indices,:])/2, axis = 0)
points = np.append(points,total_surface_lower[indices,:], axis = 0)
vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 1.0)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices)-1,3*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(3*len(indices)-1,2*len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[2*len(indices),len(indices)]]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
print('constrained_edges',constrained_edges)
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['side1'] = Member(
    id = k,     
    name = 'side1',
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = list(np.arange(len(side1_indices))) + list(np.arange(len(side1_indices),2*len(side1_indices))),
    )  

shell_mesh.optimize_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_uCRM_shellmesh_7_1')
vd_mesh = []
vd_points = []
count = 0
for memb in shell_mesh.members_dict.values():
    #print(memb.options['name'])
    conn = np.copy(memb.options['tri_connectivity']) 
    count += len(conn)
    pts = np.copy(memb.options['coordinates'])
    mesh = vedo.Mesh([pts,conn], alpha=0.5)
    mesh.backColor().lineColor('green').lineWidth(3)
    vd_mesh.append(mesh)
    if memb.options['name'] == 'upper_surface':
        pass
        # vd_points.append(vedo.Points(pts[[29,30],:], r=40, c='black'))
        # vd_points.append(vedo.Points(pts[[31,32,33],:], r=40, c='green'))
        # vd_points.append(vedo.Points(pts[[34,35,36,37,38],:], r=40, c='blue'))
        # vd_points.append(vedo.Points(pts[[39,40,41,42,43,44,45],:], r=40, c='pink'))
        #vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    elif memb.options['name'] == 'side1':
        #print(len(pts),pts)
        vd_points.append(vedo.Points(pts, r=20, c='black'))
    else:
        pass
        #vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
print('count',count)
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,vd_mesh,'222', axes=1, viewup="z", interactive = True)#,vd_points1
exit()
vd_mesh = []
vd_points = []
count = 0
for memb in shell_mesh.members_dict.values():
    #print(memb.options['name'])
    conn = np.copy(memb.options['tri_connectivity']) 
    count += len(conn)
    pts = np.copy(memb.options['coordinates'])
    mesh = vedo.Mesh([pts,conn], alpha=0.5)
    mesh.backColor().lineColor('green').lineWidth(3)
    vd_mesh.append(mesh)
    if memb.options['name'] == 'upper_surface':
        pass
        #print(pts[[40,41],:])
        #vd_points.append(vedo.Points(pts[[121,129],:], r=30, c='black'))
        #vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    elif memb.options['name'] == 'rib48':
        #print(memb.options['name'], len(pts),pts)
        vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
    elif memb.options['name'] == 'side0':
        vd_points.append(vedo.Points(pts, r=30, c='green'))
        #vd_points.append(vedo.Points(pts[421,:].reshape(1,3), r=30, c='black'))
    elif memb.options['name'] == 'side1':
        vd_points.append(vedo.Points(pts, r=30, c='black'))
    else:
        pass
        # color = list(vedo.colors.colors.values())[memb.options['id']]
        # vd_points.append(vedo.Points(pts, r=20, c=color,alpha = 0.7))
print('count',count,'num_u0_root',num_u0_root)
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,'222', axes=1, viewup="z", interactive = True)
exit()
