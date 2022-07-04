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
        if i == 111:
            num_points_u0, num_points_v0  = 4, 2 
        elif i==112:
            num_points_u0, num_points_v0  = 6, 2
        else:        
            num_points_u0, num_points_v0  = 10, 2       
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
        basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
        geo.assemble(pointset = pointset)
        geo.evaluate(pointset = pointset)    
        test_points = basis_matrix.dot(pointset.physical_coordinates)
        color = list(vedo.colors.colors.values())[i]
        vd0.append(vedo.Points(test_points/scale, r=8, c = color)) 

        indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
        indices_u0 = num_points_v0 * np.arange(num_points_u0)
        if i == 111 or i==112:
            pass
        elif i == 113:
            points_before_upper = test_points[indices_u1,:]
            points_before_lower = test_points[indices_u0,:]
        else:
            surface_upper0 = np.append(surface_upper0, np.add(points_before_upper,test_points[indices_u1,:])/2,axis=0)
            surface_lower0 = np.append(surface_lower0, np.add(points_before_lower,test_points[indices_u0,:])/2,axis=0)
            points_before_upper = test_points[indices_u1,:]
            points_before_lower = test_points[indices_u0,:] 

        surface_upper0 = np.append(surface_upper0, test_points[indices_u1,:],axis=0)
        surface_lower0 = np.append(surface_lower0, test_points[indices_u0,:],axis=0) 
           


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
vd1 = vedo.Points(surface_upper0/scale, r=15, c = 'black', alpha = 0.5)

vd2 = []
surface_upper1 = np.empty((0,3))
surface_lower1 = np.empty((0,3))
for i in range(103,108):
    pointset = OML_pointset_list[i]
    print('rib', i, pointset.pointset_id,pointset.name)
    pointset.shape= np.array([11, 11, 3])
    if i == 107:
        num_points_u0, num_points_v0  = 3, 2 
    else:
        num_points_u0, num_points_v0  = 10, 2      
    u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
    v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
    u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
    basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
    geo.assemble(pointset = pointset)
    geo.evaluate(pointset = pointset)    
    test_points = basis_matrix.dot(pointset.physical_coordinates)
    color = list(vedo.colors.colors.values())[i]
    vd2.append(vedo.Points(test_points/scale, r=10, c = color)) 
    indices_u0 = num_points_v0 * np.arange(num_points_u0)
    indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
    if i == 103:
        points_before_upper = test_points[indices_u1,:]
        points_before_lower = test_points[indices_u0,:] 
    elif i == 107:
        pass
    else:
        surface_upper1 = np.append(surface_upper1,np.add(points_before_upper,test_points[indices_u1,:])/2,axis=0)
        surface_lower1 = np.append(surface_lower1,np.add(points_before_lower,test_points[indices_u0,:])/2,axis=0)
        points_before_upper = test_points[indices_u1,:]
        points_before_lower = test_points[indices_u0,:] 

    surface_upper1 = np.append(surface_upper1, test_points[indices_u1,:],axis=0)
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


vd3 = vedo.Points(surface_upper1/scale, r=20, c = 'black', alpha = 0.5)#

total_surface_upper = np.append(surface_upper0[[0,1,2],:],surface_upper0[4:,:],axis=0)
total_surface_upper = np.append(surface_upper1, total_surface_upper,axis=0)
print('total_surface_upper',total_surface_upper.shape)
vd4 = vedo.Points(total_surface_upper/scale, r=15, c = 'blue', alpha = 0.5)

total_surface_lower = np.append(surface_lower0[[0,1,2],:],surface_lower0[4:,:],axis=0)
total_surface_lower = np.append(surface_lower1, total_surface_lower,axis=0)
print('total_surface_lower',total_surface_lower.shape)
vd5 = vedo.Points(total_surface_lower/scale, r=15, c = 'red', alpha = 0.5)


constrained_edges = np.empty((0,2),dtype = int)
for i in np.arange(0,9,1):
    constrained_edges = np.append(constrained_edges, np.array([i,i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[9,19],[19,29],[29,39],[39,49],[49,59],[59,69],[69,111]\
        ]), axis = 0)
# constrained_edges = np.append(constrained_edges, \
#     np.array([[101,91]\
#         ]), axis = 0)

num_points_uu0 = 80
num_points_vv0 = 10
in_start = 111
for i in np.arange(in_start,in_start+num_points_uu0*num_points_vv0,num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i+num_points_vv0]).reshape(1,2), axis = 0)  
for i in np.arange(in_start+num_points_uu0*num_points_vv0,in_start+num_points_uu0*num_points_vv0-num_points_vv0+1,-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
for i in np.arange(in_start+num_points_uu0*num_points_vv0-num_points_vv0+1,82,-num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i-num_points_vv0]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[82,76],[76,73],[73,70],[70,60],[60,50],[50,40],[40,30],[30,20],[20,10],[10,0]]), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[70,71],[71,72],[72,81],[81,91],[91,101]]), axis = 0)#
A = dict(vertices=total_surface_upper[:,:2], segments=constrained_edges)#
B = tr.triangulate(A,'p') 
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
indices = [70,71,72,81,91,101,111] 
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


indices0 = [9,19,29,39,49,59,69,111]  + list(np.arange(in_start,in_start+num_points_uu0*num_points_vv0+10,num_points_vv0))
print('len(indices0)',len(indices0))
points0 = np.append(total_surface_upper[indices0,:],total_surface_lower[indices0,:], axis = 0)
vd_points1 = vedo.Points(points0/scale, r=30, c='green',alpha = 0.8)
points_temp = np.empty((0,2),dtype = int)
for i in range(89):
    points_temp = np.append(points_temp, np.array([[i,10]]).reshape(1,2), axis = 0)
for i in range(89):
    points_temp = np.append(points_temp, np.array([[i,0]]).reshape(1,2), axis = 0)
A = dict(vertices=points_temp)
B = tr.triangulate(A,'pc') 
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['side0'] = Member(  
    id = k,  
    name = 'side0',
    coordinates = points0/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
k+=1

indices = [0,10,20,30,40,50,60,70,73,76,82] + list(np.arange(92,92+num_points_uu0*num_points_vv0+2*num_points_vv0,num_points_vv0))
points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 0.8)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices)-1]]).reshape(1,2), axis = 0)
for i in np.arange(2*len(indices)-1,len(indices),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['side1'] = Member(
    id = k,     
    name = 'side1',
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  

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
        #print(pts[[40,41],:])
        vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    elif memb.options['name'] == 'rib48':
        #print(len(pts),pts)
        vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
    else:
        vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
print('count',count)
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,vd_points1,vd_mesh,'222', axes=1, viewup="z", interactive = False)
shell_mesh.optimizie_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_uCRM_shellmesh_2')




