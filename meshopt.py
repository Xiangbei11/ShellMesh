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
        for i in self.itr:
            if i == 0:
                self.vertexCoords,self.trilist,self.quadlist = self.splitting()
            elif i == 1:
                self.trilist,self.quadlist = self.merging()
            elif i == 2:
                self.vertexCoords = self.qpmanual()
            elif i == 3:
                self.fullyquad()

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
        print(x[0])
        opt.setObjective(Cs @ x, GRB.MINIMIZE)
        opt.addMConstrs(A,x,'=',B)
        print('num_polygons',num_polygons, 'self.vertexCoords', self.vertexCoords.shape[0])
        print('Cs',np.shape(Cs))
        print('x', np.shape(x))
        print('A', np.shape(A))
        print('B', np.shape(B))
        print('A1', np.shape(A1))
        print('B1', np.shape(B1))
        print('A2', np.shape(A2))
        opt.optimize()
        print(x.X)
        #
        a=(np.argwhere(x.X > 0.5).flatten()-np.linspace(0,11*(num_polygons-1),num_polygons)).astype(int)
        print('a', a.shape, a)
        print('x.X', np.shape(x.X), x.X)
        print('np.argwhere(x.X > 0.5)', np.shape(np.argwhere(x.X > 0.5)), np.argwhere(x.X > 0.5))
        print('np.argwhere(x.X > 0.5).flatten()', np.shape(np.argwhere(x.X > 0.5).flatten()),np.linspace(0,11*(num_polygons-1),num_polygons))
        print('np.linspace(0,11*(num_polygons-1),num_polygons)', np.shape(np.linspace(0,11*(num_polygons-1),num_polygons)))
        exit()
        output = meshopt.splitupdatemesh(a.astype(np.int32),self.fixedvert)
        vertexCoords = output.vertlist
        trilist = output.trilist
        quadlist = output.quadlist
        self.fixedvert = output.fixedvert.astype(np.int32)
        print('splitting computenormal...')
        vertexCoords=vertexCoords.astype(np.float32)
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)
        if True:
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
    def merging(self):
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3)
        C = meshopt.create_merging_vector()
        print('C finished', C.shape)
        A_ub = meshopt.ubcon()
        print('A_ub finished')
        b_ub = np.ones(len(self.trilist))
        print('creatint bound...')
        # print()
        # print('self.fixedvert',self.fixedvert)
        # print()
        bound = meshopt.createboundsmerge(self.fixedvert)
        print('end bound...')
        opt = gp.Model("merge")
        x = opt.addMVar(shape=len(C), lb=bound[0,:], ub=bound[1,:], vtype=GRB.BINARY, name="x")
        opt.setObjective(C @ x, GRB.MINIMIZE)
        opt.addMConstrs(A_ub,x,'<',b_ub)
        a = np.eye(len(C))
        # opt.addMConstrs(a,x,'<',bound[1,:])
        # opt.addMConstrs(a,x,'>',bound[0,:])
        print('merging optimizing...')
        opt.optimize()
        res = x.X
        a = res+0.5
        output = meshopt.mergeupdatemesh(a.astype(np.int32))
        trilist = output.trilist
        quadlist = output.quadlist
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)

        if 0:
            mesh_tri = vedo.Mesh([self.vertexCoords, trilist], alpha=0.5)
            mesh_tri.backColor().lineColor('blue').lineWidth(6) 
            mesh_quad = vedo.Mesh([self.vertexCoords, quadlist], alpha=0.9)
            mesh_quad.backColor().lineColor('red').lineWidth(6) 
            mesh_ini = vedo.Mesh([self.vertexCoords, self.trilist], alpha=0.2)
            mesh_ini.backColor().lineColor('green').lineWidth(3) 
            vd_ini = vedo.Plotter(axes=1)
            vd_ini.show(mesh_ini,'Ini_mesh', viewup="z", interactive=False) 
            vd_test = vedo.Plotter(axes=1)
            print('plot finished', len(trilist),len(quadlist))
            vd_test.show(mesh_tri, mesh_quad,'optimizie_mesh', viewup="z",interactive=True)    
            exit()  
        return trilist, quadlist

    # apply the KKT conditions and get the following work (faster and is used)
    def qpmanual(self,num_itr=50, tol=1e-4, w=20, limit=0.7):
        '''
        num_itr : number of times of doing the smoothing and projection methods
        tol : if the norm of the displacement vector is smaller than tol, for loop stops
        w : weight coefficient for the regularization term
        limit : limit for the measure of the local curvature, if weight is smaller than limit, the corresponding vertex is fixed
        '''
        dupvert = self.dupvert
        vertlist1 = self.vertexCoords
        vertlist = self.vertexCoords.flatten()

        meshopt = pymeshopt.Pymeshopt(vertlist1,self.trilist,self.quadlist,self.w1,self.w2,self.w3)
        optprep = meshopt.computenormal()
        n = optprep.normal.flatten()
        weight = optprep.weight
        weight[weight>=limit] = 1/weight[weight>=limit]
        print((np.where(weight<limit)[0]).shape,'fsadfasdfasfasd')
        extraweight = np.zeros(weight.shape[0])
        extraweight[dupvert] = 1e8
        weight = weight+extraweight
        weight[weight<limit] = 1e20
        weight = np.outer(weight,np.ones(3)).flatten()
        for i in range(num_itr):
            
            num_v = len(self.vertexCoords)
            
            edges = meshopt.uniedges().astype(np.int32)
            num_e = len(edges)
            row = np.outer(np.arange(num_e*3),np.ones(2)).flatten()
            col = np.linspace(edges*3,edges*3+2,num=3,axis=1).flatten()
            data = np.repeat(np.array([[1,-1]]),num_e*3,axis=0).flatten()
            dPdd = csc_matrix((data,(row,col)),shape=(num_e*3,num_v*3))
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
            M = hstack((M,At))
            # print(M.toarray())
            v0 = vertlist1[edges[:,0]].flatten()
            v1 = vertlist1[edges[:,1]].flatten()
            v = csc_matrix(v0-v1)
            # print(v.toarray())
            b = dPdd.T.dot(v.T)
            b = vstack((b,csc_matrix((num_v,1))))
            b = b.toarray().reshape(num_v*4)*(-0.5)
            d = scipy.sparse.linalg.splu(M)
            d = d.solve(b)[:num_v*3]
            print('interation '+str(i),np.linalg.norm(d))
            vertlist = vertlist + d
            vertlist1 = vertlist.reshape(num_v,3).astype(np.float32)
            tree = ot.PyOctree(self.refv,self.reft,self.refq)
            projlist = tree.findprojectionlist(vertlist1.astype(np.float32)) 
            vertlist1 =projlist.proj.reshape(num_v,3).astype(np.float32)
            vertlist = vertlist1.flatten()
            if (np.linalg.norm(d)<tol):
                print('normd:',np.linalg.norm(d))
                print('it takes '+str(i)+' iterations to finsh qp')
                break

        fixed = vertlist1[self.fixedvert]            
        return vertlist1       

    # convert to fully quad mesh
    def fullyquad(self):
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


        
        if 0:
            # mesh_tri = vedo.Mesh([vertexCoords, trilist], alpha=0.5)
            # mesh_tri.backColor().lineColor('blue').lineWidth(6) 
            mesh_quad = vedo.Mesh([vertexCoords, quadlist], alpha=0.9)
            mesh_quad.backColor().lineColor('red').lineWidth(6) 
            # mesh_ini = vedo.Mesh([self.vertexCoords, self.trilist], alpha=0.2)
            # mesh_ini.backColor().lineColor('green').lineWidth(3) 
            # vd_ini = vedo.Plotter(axes=1)
            # vd_ini.show(mesh_ini,'Ini_mesh', viewup="z", interactive=False) 
            vd_test = vedo.Plotter(axes=1)
            print('plot finished', len(trilist),len(quadlist))
            vd_test.show(mesh_quad,'optimizie_mesh', viewup="z",interactive=True)    
            exit()

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
        
        triindex1 = np.array([0,1,2],dtype=np.int32)
        triindex2 = np.array([2,3,0],dtype=np.int32)
        cells = []
        for i in range(len(self.quadlist)):
            tuple1 = ("triangle",self.quadlist[i,triindex1].reshape(1,3))
            tuple2 = ("triangle",self.quadlist[i,triindex2].reshape(1,3))
            cells.append(tuple1)
            cells.append(tuple2)
        meshio.write_points_cells(
            vtk_file_name+"_tri.vtk",
            points,
            cells
            )