import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/pymeshopt')
import pymeshopt

import numpy as np
import scipy
from scipy.sparse import csc_matrix, coo_matrix, vstack, hstack
from scipy.optimize import linprog


import gurobipy as gp
from gurobipy import GRB
from pyoctree import pyoctree as ot

from stl import mesh
import meshio
import vedo

import matplotlib.pyplot as plt

class meshopt(object):
    def __init__(self,vertexCoords,trilist,quadlist,w1,w2,w3,itr,plot,fixedvert,dupvert=np.empty([0],dtype=np.int32)):
        '''
        vertexCoords : np.ndarray (float32) of 3 coordinates of all the vertices (n,3)
        trilist : np.ndarray (int32) of 3 indices of the vertices in all the trianlges (n,3)
        quadlist : np.ndarray (int32) of 4 indices of the vertices in all the trianlges (n,4)
        w1, w2, and w3 are weights for calculating the energy, w1 for regularity wrt angles, w2 for regularity wrt edge lengths, w3 for dihedral angles (quads)
        itr : list (int) of designed methods, 0 => splitting, 1 => merging, 2 => smoothing, 3 => transformation to fully quad mesh
        plot : 1 = plot, 0 = not plot
        fixedvert : np.ndarray (int32) of the indices of fixed vertices (n)
        refv, reft, refq : duplicates of initial vertexCoords, trilist, and quadlist
        dupvert : np.ndarray (int32) of the indices of the duplicated vertices (n) 
        '''
        self.vertexCoords = vertexCoords
        self.trilist = trilist
        self.quadlist = quadlist
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.itr = itr
        self.plot = plot
        self.fixedvert = fixedvert
        self.refv = vertexCoords.astype(np.float64)
        self.reft = trilist.astype(np.int32)
        self.refq = quadlist.astype(np.int32)
        self.dupvert = dupvert
    
    # arrange the optimzation 
    def optimization(self): 
        k = 0
        plot = False
        for i in self.itr:
            if i == 0:
                self.vertexCoords,self.trilist,self.quadlist = self.splitting()
            elif i == 1:
                self.trilist,self.quadlist = self.merging()
            elif i == 2:
                if k > 0:
                    if len(self.vertexCoords) == 39:
                        self.fixedvert = np.array([1,16,3,21,5,24,7,28,9,33,11,37,13,0,18,2,22,4,25,6,29,8,34,10,38,12,1,17,0,13,36,12]).astype(np.int32)            
                        if abs(self.vertexCoords[0,1]-3.1) < 0.001:                      
                            self.fixedvert = np.array([0, 18, 1, 21, 2, 26, 3, 30, 4, 33, 5, 37, 6, 7, 16, 8, 22, 9, 24, 10, 28, 11, 34, 12, 38, 13, 0,15,7,6,36,13 ]).astype(np.int32) #   12,36,13
                    elif len(self.vertexCoords) == 45:
                        self.fixedvert = np.array([1, 18, 3, 23, 5, 27, 7, 32, 9, 35, 11, 39, 13, 43, 15, 0, 20, 2, 24, 4, 28, 6, 30, 8, 36, 10, 40, 12, 44, 14, 1,19,0,15,42,14]).astype(np.int32)
                    elif len(self.vertexCoords) == 33:
                        self.fixedvert = np.array([1, 14, 3, 19, 5, 24, 7, 27, 9, 31, 11, 0, 16, 2, 20, 4, 22, 6, 28, 8, 32, 10, 1,15,0,11,30,10]).astype(np.int32)
                    elif len(self.vertexCoords) == 27:
                        self.fixedvert = np.array([1, 12, 3, 18, 5, 20, 7, 25, 9, 0, 14, 2, 16, 4, 21, 6, 26, 8, 1,13,0,9,24,8]).astype(np.int32)#
                    elif len(self.vertexCoords) == 21:
                        self.fixedvert = np.array([1, 10, 3, 16, 5, 19, 7, 0, 12, 2, 14, 4, 20, 6, 1,11,0,7,18,6]).astype(np.int32)#
                    elif len(self.vertexCoords) == 15:
                        self.fixedvert = np.array([1, 9, 3, 12, 5, 0, 7, 2, 13, 4, 1,10,0,5,14,4]).astype(np.int32)#
                    elif len(self.vertexCoords) == 9:
                        self.fixedvert = np.array([1, 6, 3, 0, 8, 2, 1, 7, 0, 3, 5, 2]).astype(np.int32)
                    elif len(self.vertexCoords) == 297: # or len(self.vertexCoords) == 273
                        print('self.vertexCoords',self.vertexCoords)
                        print(type(self.fixedvert),self.fixedvert)
                        fixedvert_update = []
                        fixedvert = np.array(np.where((self.vertexCoords[:,2] > 4.5 ) & (self.vertexCoords[:,1]<22))).flatten()
                        print(len(fixedvert),list(fixedvert))
                        num_points_u0 = int(len(fixedvert)/2)
                        for i in range(num_points_u0):
                            fixedvert_update.append(fixedvert[i])
                            fixedvert_update.append(fixedvert[i+num_points_u0])
                        print(len(fixedvert_update),list(fixedvert_update))

                        fixedvert = np.array(np.where((self.vertexCoords[:,2] > 4.42) & (self.vertexCoords[:,1]>22) & (self.vertexCoords[:,1]<27))).flatten()
                        print(len(fixedvert),list(fixedvert))
                        num_points_u0 = int(len(fixedvert)/2)
                        for i in range(num_points_u0):
                            fixedvert_update.append(fixedvert[i])
                            fixedvert_update.append(fixedvert[i+num_points_u0])
                        print(len(fixedvert_update),list(fixedvert_update))

                        fixedvert = np.array(np.where((self.vertexCoords[:,2] > 4.263) & (self.vertexCoords[:,1]>26.8))).flatten()
                        print(len(fixedvert),list(fixedvert))
                        num_points_u0 = int(len(fixedvert)/2)
                        for i in range(num_points_u0):
                            fixedvert_update.append(fixedvert[i])
                            fixedvert_update.append(fixedvert[i+num_points_u0])
                        print(len(fixedvert_update),list(fixedvert_update))

                        # fixedvert = np.array(np.where((self.vertexCoords[:,2] <3.7) & (self.vertexCoords[:,1]<3))).flatten()
                        # print(len(fixedvert),list(fixedvert))
                        # num_points_u0 = int(len(fixedvert)/2)
                        # for i in range(num_points_u0):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # print(len(fixedvert_update),list(fixedvert_update))

                        # fixedvert = np.array(np.where((self.vertexCoords[:,2] <4.3) & (self.vertexCoords[:,1]>3) & (self.vertexCoords[:,1]<25))).flatten()
                        # print(len(fixedvert),list(fixedvert))
                        # num_points_u0 = int(len(fixedvert)/2)
                        # for i in range(num_points_u0):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # print(len(fixedvert_update),list(fixedvert_update))

                        # fixedvert = np.array(np.where((self.vertexCoords[:,2] <4.2) & (self.vertexCoords[:,1]>25) & (self.vertexCoords[:,1]<26.8))).flatten()
                        # print(len(fixedvert),list(fixedvert))
                        # num_points_u0 = int(len(fixedvert)/2)
                        # for i in range(num_points_u0):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # print(len(fixedvert_update),list(fixedvert_update))

                        # fixedvert = np.array(np.where((self.vertexCoords[:,2] <4.1) & (self.vertexCoords[:,1]>28))).flatten()
                        # print(len(fixedvert),list(fixedvert))
                        # num_points_u0 = int(len(fixedvert)/2)
                        # for i in range(num_points_u0):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # print(len(fixedvert_update),list(fixedvert_update))

                        # fixedvert = np.array(np.where((self.vertexCoords[:,2] <4.1) & (self.vertexCoords[:,1]>29.3))).flatten()
                        # print(len(fixedvert),list(fixedvert))
                        # num_points_u0 = int(len(fixedvert)/2)
                        # for i in range(num_points_u0):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # print(len(fixedvert_update),list(fixedvert_update))

                        self.fixedvert = np.array(fixedvert_update).astype(np.int32)
                        # fixedvert_update.append(fixedvert[i+1])
                        # fixedvert = np.array(np.where(self.vertexCoords[:,2] < 3.9)).flatten()
                        # print(len(fixedvert),fixedvert)
                        # for i in range(num_points_u0-1):
                        #     fixedvert_update.append(fixedvert[i])
                        #     fixedvert_update.append(fixedvert[i+num_points_u0])
                        # fixedvert_update.append(fixedvert[i+1])                        
                        # print(len(fixedvert_update),list(fixedvert_update))
                        plot = True
                    else:
                        print('self.vertexCoords',len(self.vertexCoords))
                        plot = True
                        
                        
                self.vertexCoords = self.qpmanual()
            elif i == 3:
                k += 1
                self.fullyquad(plot)

    # splitting optimization
    def splitting(self, plot = False):
        print('Splitting setup...')
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3) 
        matrix = meshopt.create_splitting_matrix()# coefficient matrix
        num_polygons = matrix.shape[0]
        Cs = matrix.flatten()
        
        data = np.ones((11*num_polygons,1)).flatten()
        row = np.einsum('i,j->ij',np.arange(num_polygons),np.ones(11)).flatten()
        col = np.arange(len(data))
        A1 = coo_matrix((data,(row,col)), shape=[num_polygons,num_polygons*11])
        B1 = np.ones((num_polygons,1)).flatten()
        equcon=meshopt.equcon()
        data=equcon.data
        row=equcon.row
        col=equcon.col
        A2 = coo_matrix((data,(row,col)), shape=[equcon.sharededges*3,num_polygons*11])
        B2=equcon.b
        A = scipy.sparse.vstack([A1,A2])
        B = np.concatenate((B1,B2)).flatten()
        bound = meshopt.createbounds(self.fixedvert)

        opt = gp.Model("split")
        x = opt.addMVar(shape=len(Cs),lb=bound[0,:], ub=bound[1,:], vtype=GRB.BINARY, name="x")
        opt.setObjective(Cs @ x, GRB.MINIMIZE)
        opt.addMConstrs(A,x,'=',B)
        a = np.eye(len(Cs))
        print('splitting optimizing start...')
        opt.optimize()
        print('splitting optimize finish')
        a=(np.argwhere(x.X > 0.5).flatten()-np.linspace(0,11*(num_polygons-1),num_polygons)).astype(int)
        # print('a', a.shape, a)
        # print('x.X', np.shape(x.X), x.X)
        # print('np.argwhere(x.X > 0.5)', np.shape(np.argwhere(x.X > 0.5)), np.argwhere(x.X > 0.5))
        # print('np.argwhere(x.X > 0.5).flatten()', np.shape(np.argwhere(x.X > 0.5).flatten()),np.linspace(0,11*(num_polygons-1),num_polygons))
        # print('np.linspace(0,11*(num_polygons-1),num_polygons)', np.shape(np.linspace(0,11*(num_polygons-1),num_polygons)))
        # exit()
        output = meshopt.splitupdatemesh(a.astype(np.int32),self.fixedvert)
        vertexCoords = output.vertlist
        trilist = output.trilist
        quadlist = output.quadlist
        self.fixedvert = output.fixedvert.astype(np.int32)
        print('splitting computenormal...')
        vertexCoords=vertexCoords.astype(np.float32)
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)
        if 0:
            mesh_tri = vedo.Mesh([vertexCoords, trilist], alpha=0.5)
            mesh_tri.backColor().lineColor('blue').lineWidth(6) 
            mesh_quad = vedo.Mesh([vertexCoords, quadlist], alpha=0.9)
            mesh_quad.backColor().lineColor('red').lineWidth(6) 
            mesh_ini = vedo.Mesh([self.vertexCoords, self.trilist], alpha=0.2)
            mesh_ini.backColor().lineColor('green').lineWidth(3) 
            vd_ini = vedo.Plotter(axes=1)
            vd_ini.show(mesh_ini,'Ini_mesh', viewup="z", interactive=False) 
            vd_test = vedo.Plotter(axes=1)
            print('plot finished', len(trilist),len(quadlist))
            vd_test.show(mesh_tri, mesh_quad,mesh_ini,'optimizie_mesh', viewup="z",interactive=True)    
            exit()         
        return vertexCoords, trilist, quadlist

    # merging optimization
    def merging(self): #export GRB_LICENSE_FILE=/Users/Sansara/Public/Code/A/gurobi.lic     
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3)
        C = meshopt.create_merging_vector()
        A_ub = meshopt.ubcon()
        b_ub = np.ones(len(self.trilist))
        bound = meshopt.create_merging_bounds(self.fixedvert)
        opt = gp.Model("merge")
        x = opt.addMVar(shape=len(C), lb=bound[0,:], ub=bound[1,:], vtype=GRB.BINARY, name="x")
        opt.setObjective(C @ x, GRB.MINIMIZE)
        opt.addMConstrs(A_ub,x,'<',b_ub)
        a = np.eye(len(C))
        opt.optimize()
        res = x.X
        a = res+0.5
        output = meshopt.mergeupdatemesh(a.astype(np.int32))
        trilist = output.trilist
        quadlist = output.quadlist
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)
        
        if 0:
            connectivity = trilist
            for i in range(len(connectivity)):
                index = 29
                if index in list(connectivity[i,:]): 
                    print(i, index, connectivity[i,:]) 
                # index = 50
                # if index in list(connectivity[i,:]): 
                #     print(i, index, connectivity[i,:])    
            print()    
            connectivity = quadlist
            for i in range(len(connectivity)):
                index = 29
                if index in list(connectivity[i,:]): 
                    print(i, index, connectivity[i,:]) 
                # index = 50
                # if index in list(connectivity[i,:]): 
                #     print(i, index, connectivity[i,:])  
         
            mesh_tri = vedo.Mesh([self.vertexCoords, trilist], alpha=0.5)
            mesh_tri.backColor().lineColor('blue').lineWidth(6) 
            mesh_quad = vedo.Mesh([self.vertexCoords, quadlist], alpha=0.9)
            mesh_quad.backColor().lineColor('red').lineWidth(6) 
            mesh_ini = vedo.Mesh([self.vertexCoords, self.trilist], alpha=0.2)
            mesh_ini.backColor().lineColor('green').lineWidth(3) 
            vd_ini = vedo.Plotter(axes=1)
            index = [44,29]
            print('index',index)
            points = vedo.Points(self.vertexCoords[index,:], r=25, c='blue',alpha = 0.3)
            vd_ini.show(mesh_ini,points,'Ini_mesh', viewup="z", interactive=False) 
            
            print()
            print('self.fixedvert',len(self.fixedvert),self.fixedvert)
            print()
            print('plot finished', len(self.trilist),len(trilist),len(quadlist),len(self.vertexCoords))
            vd_test = vedo.Plotter(axes=1)
            vd_test.show(mesh_tri, mesh_quad,'optimizie_mesh', viewup="z",interactive=True)  
            exit()  
        return trilist, quadlist

    # apply the KKT conditions and get the following work (faster and is used)
    def qpmanual(self,num_itr=50, tol=1e-4, w=20, limit=0.7):#
        '''
        num_itr : number of times of doing the smoothing and projection methods
        tol : if the norm of the displacement vector is smaller than tol, for loop stops
        w : weight coefficient for the regularization term
        limit : limit for the measure of the local curvature, if weight is smaller than limit, the corresponding vertex is fixed
        '''
        dupvert = self.dupvert
        print('dupvert',dupvert)
        vertlist = self.vertexCoords
        vertlist_flatten = self.vertexCoords.flatten()

        meshopt = pymeshopt.Pymeshopt(vertlist,self.trilist,self.quadlist,self.w1,self.w2,self.w3)
        optprep = meshopt.computenormal()
        n = optprep.normal.flatten()
        #print('n',n.shape)
        weight = optprep.weight
        weight[weight>=limit] = 1/weight[weight>=limit]
        #print('weight',weight.shape,weight)
        extraweight = np.zeros(weight.shape[0])
        extraweight[dupvert] = 1e8
        weight = weight+extraweight
        weight[weight<limit] = 1e20
        weight = np.outer(weight,np.ones(3)).flatten()
        #print('weight',weight.shape,weight)
        for i in range(num_itr):    
            num_v = len(self.vertexCoords)
            edges = meshopt.uniedges().astype(np.int32)
            num_e = len(edges)
            # print('edges',edges)
            # print('num_v',num_v,'num_e',num_e)
            row = np.outer(np.arange(num_e*3),np.ones(2)).flatten()
            col = np.linspace(edges*3,edges*3+2,num=3,axis=1).flatten()
            #print('col',np.linspace(edges*3,edges*3+2,num=3,axis=1))
            data = np.repeat(np.array([[1,-1]]),num_e*3,axis=0).flatten()
            dPdd = csc_matrix((data,(row,col)),shape=(num_e*3,num_v*3))
            #print('dPdd',np.shape(dPdd.toarray()),dPdd.toarray())
            v0 = vertlist[edges[:,0]].flatten()
            v1 = vertlist[edges[:,1]].flatten()
            v = csc_matrix(v0-v1)
            #print('v', np.shape(v.toarray()))
            b = dPdd.T.dot(v.T)
            b = vstack((b,csc_matrix((num_v,1))))
            b = b.toarray().reshape(num_v*4)*(-0.5)

            row = col = np.arange(num_v*3)
            data = np.ones(num_v*3)*w
            index = np.linspace(self.fixedvert*3,self.fixedvert*3+2,num=3).flatten().astype(np.int32)
            data[index] = 1e8
            data = data*weight
            B = dPdd.T.dot(dPdd) + csc_matrix((data,(row,col)),shape=(num_v*3,num_v*3))
            row = np.outer(np.arange(num_v),np.ones(3)).flatten()
            col = np.arange(num_v*3).flatten()
            data = n
            A = csc_matrix((data,(row,col)),shape=(num_v,num_v*3))
            M = vstack((B,A))
            At = vstack((A.T,csc_matrix((num_v,num_v))))
            M = csc_matrix(hstack((M,At)))
            # print(M.toarray())


            d = scipy.sparse.linalg.splu(M)
            d = d.solve(b)
            d = d[:num_v*3]
            if i%24==0:
                print('interation '+str(i),'normd:',np.linalg.norm(d))
            vertlist_flatten = vertlist_flatten + d
            vertlist = vertlist_flatten.reshape(num_v,3).astype(np.float32)
            tree = ot.PyOctree(self.refv,self.reft,self.refq)
            projlist = tree.findprojectionlist(vertlist.astype(np.float32)) 
            vertlist =projlist.proj.reshape(num_v,3).astype(np.float32)
            vertlist_flatten = vertlist.flatten()
            if (np.linalg.norm(d)<tol):
                print('normd:',np.linalg.norm(d))
                print('it takes '+str(i)+' iterations to finsh qp')
                break
            # print('d',d.shape)
            # print('vertlist_flatten',vertlist_flatten.shape)
            # exit()         
        return vertlist       

    # convert to fully quad mesh
    def fullyquad(self,plot):
        print('splitting quads and tris:')
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3)     
        fullyquadmesh = meshopt.fullyquadmesh(self.fixedvert)
        
        self.trilist = np.empty((0,3),dtype=np.int32)
        self.quadlist = fullyquadmesh.quadlist.astype(np.int32)
        self.vertexCoords = fullyquadmesh.vertlist.astype(np.float32)
        self.fixedvert = fullyquadmesh.fixedvert.astype(np.int32)
        vertexCoords = self.vertexCoords
        trilist = self.trilist
        quadlist = self.quadlist
        fixed = self.vertexCoords[self.fixedvert]
        print('fully quad mesh generated',len(trilist),len(quadlist))


        
        if plot:
            # mesh_tri = vedo.Mesh([vertexCoords, trilist], alpha=0.5)
            # mesh_tri.backColor().lineColor('blue').lineWidth(6) 
            mesh_quad = vedo.Mesh([vertexCoords, quadlist], alpha=0.9)
            mesh_quad.backColor().lineColor('red').lineWidth(6) 
            # mesh_ini = vedo.Mesh([self.vertexCoords, self.trilist], alpha=0.2)
            # mesh_ini.backColor().lineColor('green').lineWidth(3) 
            points = vedo.Points(vertexCoords[self.fixedvert,:],r = 20,c ='green')
            # vd_ini = vedo.Plotter(axes=1)
            # vd_ini.show(mesh_ini,'Ini_mesh', viewup="z", interactive=False) 
            vd_test = vedo.Plotter(axes=1)
            print('plot finished', len(trilist),len(quadlist))
            vd_test.show(mesh_quad,points,'optimizie_mesh',interactive=True)    
            #exit()

    # save as stl (not used since stl only contains triangles)
    def saveasstl(self,stl_file_name):
        trilist = self.trilist
        quadlist = self.quadlist
        vertexCoords = self.vertexCoords
        quad2tri = np.concatenate((quadlist[:,:-1],quadlist[:,2:],quadlist[:,0].reshape(len(quadlist),1)),axis=1).reshape(len(quadlist)*2,3)
        tris = np.concatenate((trilist,quad2tri),axis=0)
        data = vertexCoords[(tris.flatten())].reshape(len(tris),3,3)
        cube = mesh.Mesh(np.zeros(data.shape[0], dtype=mesh.Mesh.dtype))
        cube.vectors[:, :, :] = data
        print('save as stl file...')
        cube.save('{}.stl'.format(stl_file_name))

    # save as vtk (mostly often used)
    def saveasvtk(self,vtk_file_name):
        '''
        tri: True => save as triangles, False => save as quads
        '''
        print('save as vtk')
        points = self.vertexCoords
        cells = []
        # print(len(points))
        for i in range(len(self.quadlist)):
            tuple = ("quad",self.quadlist[i,:].reshape(1,4))
            cells.append(tuple)
        meshio.write_points_cells(
            vtk_file_name+"_quad.vtk",
            points,
            cells
            )
        
        # triindex1 = np.array([0,1,2],dtype=np.int32)
        # triindex2 = np.array([2,3,0],dtype=np.int32)
        # cells = []
        # for i in range(len(self.quadlist)):
        #     tuple1 = ("triangle",self.quadlist[i,triindex1].reshape(1,3))
        #     tuple2 = ("triangle",self.quadlist[i,triindex2].reshape(1,3))
        #     cells.append(tuple1)
        #     cells.append(tuple2)
        # meshio.write_points_cells(
        #     vtk_file_name+"_tri.vtk",
        #     points,
        #     cells
        #     )