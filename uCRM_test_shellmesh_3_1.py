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
            num_points_u0, num_points_v0  = 8, 2
        else:        
            num_points_u0, num_points_v0  = 15, 2       
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T    
        basis_matrix = shell_mesh.discritize_ctrl_pointset(pointset, uv_vec = u0_v0_vec)
        geo.assemble(pointset = pointset)
        geo.evaluate(pointset = pointset)    
        test_points = basis_matrix.dot(pointset.physical_coordinates)
        if k==12:
            print('k==12',i)
            print(len(test_points),test_points)
            test_points[4,:] = np.array([32980,10458,4055])
            test_points[5,:] = np.array([32980,10458,4810])        
        if k==13:
            print()
            print('k==13',i)
            print(len(test_points),test_points)
            test_points[2,:] = np.array([32981,11303,4144])
            test_points[3,:] = np.array([32981,11303,4786])
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
            # 
            #     points_before_upper[2,:] = points_before_upper[3,:]
            surface_upper_in = np.add(points_before_upper,test_points[indices_u1,:])/2
            surface_lower_in = np.add(points_before_lower,test_points[indices_u0,:])/2
            if k ==13:
                surface_upper_in[2,:] = np.array([32980.5,10880.5,4796])
                surface_lower_in[2,:] = np.array([32980.5,10880.5,4099.5])
            surface_upper0 = np.append(surface_upper0, surface_upper_in,axis=0)
            surface_lower0 = np.append(surface_lower0, surface_lower_in,axis=0)
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
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd0,vd1,'000', axes=1, viewup="z", interactive = True)
# exit()
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
        num_points_u0, num_points_v0  = 15, 2      
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
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd0,vd2,vd3,'000', axes=1, viewup="z", interactive = True)
# exit()
total_surface_upper = np.append(surface_upper0[[0,1,2],:],surface_upper0[4:,:],axis=0)
total_surface_upper = np.append(surface_upper1, total_surface_upper,axis=0)
print('total_surface_upper',total_surface_upper.shape)
vd4 = vedo.Points(total_surface_upper/scale, r=15, c = 'blue', alpha = 0.5)
total_surface_lower = np.append(surface_lower0[[0,1,2],:],surface_lower0[4:,:],axis=0)
total_surface_lower = np.append(surface_lower1, total_surface_lower,axis=0)
print('total_surface_lower',total_surface_lower.shape)
vd5 = vedo.Points(total_surface_lower/scale, r=15, c = 'red', alpha = 0.5)
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd0,vd2,vd4,vd5,'111', axes=1, viewup="z", interactive = False)
# exit()

constrained_edges = np.empty((0,2),dtype = int)
for i in np.arange(0,14,1):
    constrained_edges = np.append(constrained_edges, np.array([i,i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[14,29],[29,44],[44,59],[59,74],[74,89],[89,104],[104,163]\
        ]), axis = 0)

num_points_uu0 = 80
num_points_vv0 = 15
in_start = 163
for i in np.arange(in_start,in_start+num_points_uu0*num_points_vv0,num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i+num_points_vv0]).reshape(1,2), axis = 0)  
for i in np.arange(in_start+num_points_uu0*num_points_vv0,in_start+num_points_uu0*num_points_vv0-num_points_vv0+1,-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
for i in np.arange(in_start+num_points_uu0*num_points_vv0-num_points_vv0+1,119,-num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i-num_points_vv0]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[119,111],[111,108],[108,105]]), axis = 0)
for i in np.arange(105,0,-num_points_vv0):
    constrained_edges = np.append(constrained_edges, np.array([i, i-num_points_vv0]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, \
    np.array([[105,106],[106,107],[107,118],[118,133]]), axis = 0)#
# print('constrained_edges',constrained_edges)
# A = dict(vertices=total_surface_upper[:,:2], segments=constrained_edges)
# B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
A = dict(vertices=total_surface_upper[:,:2], segments=constrained_edges)#
B = tr.triangulate(A,'p') 
connectivity = np.copy(B['triangles'])

print('total_surface_upper',len(total_surface_upper),'connectivity',len(connectivity))
indices = [105,106,107,118,133,148,163] 
total_surface_upper[indices,1] = 3.1*scale
total_surface_lower[indices,1] = 3.1*scale
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

points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
ids = [0,2]
A = dict(vertices=points[:,ids])
B = tr.triangulate(A,'pc') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['rib'+'%s'%i] = Member(   
    id = k, 
    name = 'rib'+'%s'%(k-2),
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
print('rib',k-2)
k+=1
i+=1

indices =[421,436,450]
points = np.append(total_surface_upper[indices,:],total_surface_lower[indices,:], axis = 0)
ids = [1,2]
A = dict(vertices=points[:,ids])
B = tr.triangulate(A,'pc') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['rib'+'%s'%i] = Member(   
    id = k, 
    name = 'rib'+'%s'%(k-2),
    coordinates = points/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
print('rib',k-2)
k+=1

indices0 = [14,29,44,59,74,89,104,163]  + list(np.arange(in_start,in_start+num_points_uu0*num_points_vv0+num_points_vv0,num_points_vv0))
print('len(indices0)',len(indices0))
points0 = np.append(total_surface_upper[indices0,:],total_surface_lower[indices0,:], axis = 0)
vd_points1 = vedo.Points(points0/scale, r=30, c='green',alpha = 0.8)
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd1,vd3,vd_points1,'111', axes=1, viewup="z", interactive = True)
# exit()
points_temp = np.empty((0,2),dtype = int)
for i in range(89):
    points_temp = np.append(points_temp, np.array([[i,10]]).reshape(1,2), axis = 0)
for i in range(89):
    points_temp = np.append(points_temp, np.array([[i,0]]).reshape(1,2), axis = 0)
A = dict(vertices=points_temp)
B = tr.triangulate(A,'pc') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
connectivity = np.copy(B['triangles'])
shell_mesh.members_dict['side0'] = Member(  
    id = k,  
    name = 'side0',
    coordinates = points0/scale,
    tri_connectivity = connectivity,
    constrained_node_indices = [],
    )  
k+=1

indices1 = list(np.arange(0,105,num_points_vv0)) +[105,108,111] + list(np.arange(119,119+num_points_uu0*num_points_vv0+3*num_points_vv0,num_points_vv0))#
print('len(indices1)',len(indices1))
points = np.append(total_surface_upper[indices1,:],total_surface_lower[indices1,:], axis = 0)
vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 0.8)
# vd_plotter1 = vedo.Plotter()
# vd_plotter1.show(vd1,vd3,vd_points1,'111', axes=1, viewup="z", interactive = True)
# exit()
points = np.append(total_surface_upper[indices1,:],total_surface_lower[indices1,:], axis = 0)
vd_points1 = vedo.Points(points/scale, r=30, c='green',alpha = 0.8)
ids = [1,2]
constrained_edges = np.empty((0,2),dtype = int)
for i in range(len(indices1)-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i+1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[i+1,2*len(indices1)-1]]).reshape(1,2), axis = 0)
for i in np.arange(2*len(indices1)-1,len(indices1),-1):
    constrained_edges = np.append(constrained_edges, np.array([i, i-1]).reshape(1,2), axis = 0)
constrained_edges = np.append(constrained_edges, np.array([[len(indices1),0]]).reshape(1,2), axis = 0)
A = dict(vertices=points[:,ids], segments=constrained_edges)
B = tr.triangulate(A,'p') 
# tr.compare(plt, A, B)         
# plt.show()#block=False
# exit()
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
        #vd_points.append(vedo.Points(pts[293,:].reshape(1,3), r=40, c='blue'))
        vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    elif memb.options['name'] == 'rib49':
        #print(memb.options['name'], len(pts),pts)
        vd_points.append(vedo.Points(pts, r=30, c='red',alpha = 0.7))
    elif memb.options['name'] == 'lower_surface':
        vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
        #vd_points.append(vedo.Points(pts[303,:].reshape(1,3), r=40, c='blue'))
    else:
        pass
        #vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
print('count',count)

shell_mesh.optimizie_mesh()
shell_mesh.construct_whole_structure_optmesh('CAD_uCRM_shellmesh_3')
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,vd_points1,vd_mesh,'222', axes=1, viewup="z", interactive = True)



exit()
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
        vd_points.append(vedo.Points(pts[[105,106,107,118,133,148,163],:], r=30, c='black'))
        vd_points.append(vedo.Points(pts, r=25, c='blue',alpha = 0.3))
    elif memb.options['name'] == 'rib49':
        #print(memb.options['name'], len(pts),pts)
        vd_points.append(vedo.Points(pts, r=30, c='red',alpha = 0.7))
    elif memb.options['name'] == 'lower_surface':
        vd_points.append(vedo.Points(pts, r=25, c='green',alpha = 0.3))
        vd_points.append(vedo.Points(pts[421,:].reshape(1,3), r=30, c='black'))
    else:
        pass
        #vd_points.append(vedo.Points(pts, r=20, c='red',alpha = 0.7))
print('count',count)
vd_plotter1 = vedo.Plotter()
vd_plotter1.show(vd_points,vd_mesh,'222', axes=1, viewup="z", interactive = True)
exit()