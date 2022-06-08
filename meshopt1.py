import numpy as np
import scipy
from scipy.sparse import csc_matrix, coo_matrix, vstack, hstack, csr_matrix
import vtk
import meshio
import vedo

import gurobipy as gp
from gurobipy import GRB

import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/pymeshopt')
import pymeshopt 
from pyoctree import pyoctree as ot

class meshopt(object):
    def __init__(self,vertexCoords,trilist,quadlist,fixedvert,w1,w2,w3,w4,itr,plot = False):
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
        self.w4 = w4
        self.itr = itr
        self.plot = plot
        self.fixedvert = fixedvert
        self.dupvert = np.empty([0],dtype=np.int32)
        self.refq = quadlist.astype(np.int32)
        self.refv = vertexCoords.astype(np.float64)
        self.reft = trilist.astype(np.int32)
    
    # arrange the optimzation 
    def optimization(self):
        for i in self.itr:
            if i == 0:
                print('START splitting')
                self.vertexCoords,self.trilist,self.quadlist = self.splitting()   
            elif i == 1:
                print('START merging')
                self.trilist, self.quadlist = self.merging()
            elif i == 2:
                print('START smoothing')
                self.vertexCoords = self.qpmanual()
            elif i == 3:
                print('START quad transformation')
                self.fullyquad()
        if self.plot:
            pass

    # merging optimization
    def merging(self):
        vertexCoords = self.vertexCoords
        trilist = self.trilist
        quadlist = self.quadlist
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist, self.w1,self.w2,self.w3)
        C = meshopt.creatematrixmerge()
        print('C finished')
        A_ub = meshopt.ubcon()
        print('A_ub finished')
        b_ub = np.ones(len(trilist))
        print('creatint bound...')
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
        edges = meshopt.uniedges()
        edges = edges.astype(int)
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist, self.w1,self.w2,self.w3)
        return trilist, quadlist

    # splitting optimization
    def splitting(self):
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3) 
        optprep = meshopt.computenormal()
        normal = optprep.normal# normal vector on each vertex
        cs = meshopt.computecs()# local coordinate system on each vertex
        # print out normals to check if they are valid
        print('normal',np.argwhere(np.isnan(normal)))
        print('normal',np.isnan(np.min(normal)))
        print('cs',np.argwhere(np.isnan(cs)))
        print('splitting setup...')
        vertexCoords = self.vertexCoords
        fixed = vertexCoords[self.fixedvert]
        trilist = self.trilist
        quadlist = self.quadlist
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist,self.w1,self.w2,self.w3)
        matrix = meshopt.creatematrix()## coefficient matrix
        nel = matrix.shape[0]
        print('nel',nel)
        C = matrix.flatten()
        num_poly = len(trilist)+len(quadlist)
        num_tri = len(trilist)
        num_quad = len(quadlist)
        # create linear constraint
        data = np.ones((11*nel,1)).flatten()
        row = np.einsum('i,j->ij',np.arange(nel),np.ones(11)).flatten()
        col = np.arange(len(data))
        A1 = coo_matrix((data,(row,col)), shape=[nel,nel*11])
        B1 = np.ones((num_poly,1)).flatten()
        equcon=meshopt.equcon()
        data=equcon.data
        row=equcon.row
        col=equcon.col
        A2 = coo_matrix((data,(row,col)), shape=[equcon.sharededges*3,nel*11])
        B2=equcon.b
        A = scipy.sparse.vstack([A1,A2])
        B = np.concatenate((B1,B2)).flatten()
        bound = meshopt.createbounds(self.fixedvert)
        print('bound',bound)
        print('splitting optimizing...')
        opt = gp.Model("split")
        print('shape(C)',np.shape(C))
        print('len(C)',len(C))
        #exit()
        x = opt.addMVar(shape=len(C),lb=bound[0,:], ub=bound[1,:], vtype=GRB.BINARY, name="x")
        opt.setObjective(C @ x, GRB.MINIMIZE)
        opt.addMConstrs(A,x,'=',B)
        a = np.eye(len(C))
        print('splitting optimizing start...')
        opt.optimize()
        print('splitting optimize finish')
        a=(np.argwhere(x.X > 0.5).flatten()-np.linspace(0,11*(nel-1),nel)).astype(int)
        output = meshopt.splitupdatemesh(a.astype(np.int32),self.fixedvert)
        vertexCoords = output.vertlist
        splitcoord = output.splitcoord
        trilist = output.trilist
        quadlist = output.quadlist
        self.fixedvert = output.fixedvert.astype(np.int32)
        print('splitting computenormal...')
        vertexCoords=vertexCoords.astype(np.float32)
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)
        return vertexCoords, trilist, quadlist



    # apply the KKT conditions and get the following work (faster and is used)
    def qpmanual(self,num_itr=8,tol=1e-6, w=0.01,limit=0.7):#  limit=0.7,tol=1e-4, w=20, 
        '''
        num_itr : number of times of doing the smoothing and projection methods
        tol : if the norm of the displacement vector is smaller than tol, for loop stops
        w : weight coefficient for the regularization term
        limit : limit for the measure of the local curvature, if weight is smaller than limit, the corresponding vertex is fixed
        '''
        dupvert = self.dupvert
        vertlist1 = self.vertexCoords
        print('self.vertexCoords',np.shape(vertlist1))
        vertlist = self.vertexCoords.flatten()

        print('self.trilist',np.shape(self.trilist))
        print('self.quadlist',np.shape(self.quadlist),self.quadlist)
        for z in range(len(self.trilist)):
            for zz in range(3):
                if self.trilist[z,zz] == 4932:
                    print(self.trilist[z,:])

        meshopt = pymeshopt.Pymeshopt(vertlist1,self.trilist,self.quadlist,self.w1,self.w2,self.w3)
        optprep = meshopt.computenormal()
        n = optprep.normal.flatten()
        
        # for z in range(len(self.vertexCoords)):
        #     if np.isnan(np.sum(optprep.normal[z,:])):
        #         print('Nan element in',z,optprep.normal[z,:],vertlist1[z,:])
        #         plot = z
        # print('n',np.shape(n),np.sum(n))
        # print('optprep.normal',np.shape(optprep.normal))
        # point11 = vedo.Points([vertlist1[plot,:]], r= 8, c='red')
        # mesh = vedo.Mesh([self.vertexCoords, self.trilist])#,alpha = 0.1
        # mesh.backColor().lineColor('green').lineWidth(2)
        # tri_mesh = vedo.Plotter()
        # tri_mesh.show('Structual triangular mesh',mesh, point11, viewup='z', axes=1, interactive = True)
        # exit()
        
        weight = optprep.weight
        weight[weight>=limit] = 1/weight[weight>=limit]
        #print((np.where(weight<limit)[0]).shape,'fsadfasdfasfasd')
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
            #print(np.shape(A.toarray))    
            #print('A',np.shape(A),np.sum(A))
            #print('B',np.shape(B),np.sum(B))
            #print(np.shape(B.todense))
            M = vstack((B,A))
            At = vstack((A.T,csc_matrix((num_v,num_v))))
            print(np.shape(At))
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
            vertlist = vertlist + d
            vertlist1 = vertlist.reshape(num_v,3).astype(np.float32)     
            tree = ot.PyOctree(self.refv,self.reft,self.refq)
            projlist = tree.findprojectionlist(vertlist1.astype(np.float32)) 
            vertlist1 =projlist.proj.reshape(num_v,3).astype(np.float32)
            vertlist = vertlist1.flatten()
            norm_d = np.linalg.norm(d)
            print('interation '+str(i), 'norm_d',norm_d)
            if (norm_d<tol):
                print('normd:',np.linalg.norm(d))
                print('it takes '+str(i)+' iterations to finsh qp')
                break
        # ###plot
        # mesh = vedo.Mesh([vertlist1, self.trilist], alpha=0.1)
        # mesh.backColor().lineColor('blue').lineWidth(3)
        # tri_mesh2 = vedo.Plotter()
        # tri_mesh2.show(mesh, 'After smoothing',viewup='z', axes=1, interactive = True)
    
        return vertlist1      

    # convert to fully quad mesh
    def fullyquad(self):
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3)     
        fullyquadmesh = meshopt.fullyquadmesh(self.fixedvert)
        self.trilist = np.empty((0,3),dtype=np.int32)
        self.quadlist = fullyquadmesh.quadlist.astype(np.int32)
        self.vertexCoords = fullyquadmesh.vertlist.astype(np.float32)
        self.fixedvert = fullyquadmesh.fixedvert.astype(np.int32)

    def saveasvtk(self, vtk_file_name):
        points = self.vertexCoords
        quadlist = self.quadlist
        cells = []
        triindex1 = np.array([0,1,2],dtype=np.int32)
        triindex2 = np.array([2,3,0],dtype=np.int32)
        for i in range(len(self.quadlist)):
            tuple1 = ('triangle',quadlist[i,triindex1].reshape(1,3))
            tuple2 = ('triangle',quadlist[i,triindex2].reshape(1,3))
            cells.append(tuple1)
            cells.append(tuple2)
        meshio.write_points_cells(
            vtk_file_name+'_tri.vtk',
            points,
            cells
            ) 
        cells = []
        for i in range(len(quadlist)):
            tuple = ('quad',quadlist[i,:].reshape(1,4))
            cells.append(tuple)
        meshio.write_points_cells(
            vtk_file_name+'_quad.vtk',
            points,
            cells
            ) 
        print('Save as vtk finished:',vtk_file_name)  

