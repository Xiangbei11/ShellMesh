import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/pymeshopt')
import pymeshopt 
import numpy as np
import scipy
from scipy.sparse import csc_matrix
from scipy.optimize import linprog
import matplotlib.pyplot as plt
import sys, vtk
from mpl_toolkits import mplot3d
import time
import gurobipy as gp
from gurobipy import GRB
from pyoctree import pyoctree as ot
from openmdao.api import Problem, ScipyOptimizeDriver,ExecComp, IndepVarComp, ExplicitComponent, OptionsDictionary, Group, pyOptSparseDriver
from scipy.sparse import coo_matrix
# from verticescomp import VerticesComp
# from projectioncomp import ProjectionComp
# from energycomp import EnergyComp
# from regularizationcomp import RegularizationComp
# from costcomp import CostComp
# from dvertcomp import DvertComp
from stl import mesh
import os.path
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix
import meshio

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
        if self.plot:
            plt.show()

    # splitting optimization
    def splitting(self):
        '''
        filename : (string) file name to save the plotted figure
        '''
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
        w1 = self.w1
        w2 = self.w2
        w3 = self.w3
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist,w1,w2,w3)
        matrix = meshopt.creatematrix()# coefficient matrix
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
        if self.plot:
            # plt.figure()

            # for i in range(len(trilist)):
            #     coord = np.concatenate((vertexCoords[trilist[i,:]],np.array([vertexCoords[trilist[i,0]]])))
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'b')

            # for i in range(len(quadlist)):
            #     coord = np.concatenate((vertexCoords[quadlist[i,:]],np.array([vertexCoords[quadlist[i,0]]])))
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'b')
            
            # for i in range(int(len(splitcoord)/2)):
            #     coord = splitcoord[i*2:i*2+2,:]
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,color='r')
            # plt.plot(fixed[:,0],fixed[:,1],'k-o')
            # plt.axis('equal')

            # plt.axis('off')
            # plt.savefig("splitting.pdf", bbox_inches='tight')
            print('splitting plotting...')
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            num_tri = len(trilist)
            num_quad = len(quadlist)
            for i in range(num_tri):
                coord = np.concatenate((vertexCoords[trilist[i]],np.array([vertexCoords[trilist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')
            for i in range(num_quad):
                coord = np.concatenate((vertexCoords[quadlist[i,:]],np.array([vertexCoords[quadlist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')
            for i in range(int(len(splitcoord)/2)):
                coord = splitcoord[i*2:i*2+2,:]
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='r')
            for i in range(len(vertexCoords)):
                coord = vertexCoords[i,:]
                x1 = coord[0]
                y1 = coord[1]
                z1 = coord[2]
                coord = normal[i]
                x2 = coord[0]
                y2 = coord[1]
                z2 = coord[2]
                # ax.quiver(x1,y1,z1,x2,y2,z2)
            print('---------------------------------',fixed.shape[0])
            ax.plot(fixed[:,0],fixed[:,1],fixed[:,2],'k-o')
            X = vertexCoords[:,0]
            Y = vertexCoords[:,1]
            Z = vertexCoords[:,2]
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
            # plt.savefig(filename)
            print('splitting plot finish')
        vertexCoords=vertexCoords.astype(np.float32)
        trilist=trilist.astype(np.int32)
        quadlist=quadlist.astype(np.int32)
        # print('----------------------------',len(self.fixedvert))

        # meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist,self.w1,self.w2,self.w3)    
        
        # optprep = meshopt.computenormal()
        # normal = optprep.normal
        # cs = meshopt.computecs()   
        # print('normal',np.argwhere(np.isnan(normal)))
        # print('normal',np.isnan(np.min(normal)))
        # print('cs',np.argwhere(np.isnan(cs)))
        return vertexCoords, trilist, quadlist

    # merging optimization
    def merging(self):
        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,self.trilist,self.quadlist,self.w1,self.w2,self.w3) 
        optprep = meshopt.computenormal()
        normal = optprep.normal
        cs = meshopt.computecs() 
        print('normal',np.argwhere(np.isnan(normal)))
        print('normal',np.isnan(np.min(normal)))
        print('cs',np.argwhere(np.isnan(cs)))
        print('MMMMMmerging setup...')
        vertexCoords = self.vertexCoords
        trilist = self.trilist
        quadlist = self.quadlist
        w1 = self.w1
        w2 = self.w2
        w3 = self.w3
        if len(trilist)==0:
            return trilist,quadlist
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist,w1,w2,w3)
        C = meshopt.creatematrixmerge()
        print('C finished', C.shape)
        #exit()
        A_ub = meshopt.ubcon()
        print('A_ub finished')
        b_ub = np.ones(len(trilist))
        print('creatint bound...')
        print()
        print('self.fixedvert',self.fixedvert)
        print()
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
        meshopt = pymeshopt.Pymeshopt(vertexCoords,trilist,quadlist,w1,w2,w3)
        # optprep = meshopt.computenormal()
        # normal = optprep.normal
        fixed = vertexCoords[self.fixedvert]

        meshopt = pymeshopt.Pymeshopt(self.vertexCoords,trilist,quadlist,self.w1,self.w2,self.w3)    
        optprep = meshopt.computenormal()
        normal = optprep.normal
        cs = meshopt.computecs()   
        print('normal',np.argwhere(np.isnan(normal)))
        print('normal',np.isnan(np.min(normal)))
        print('cs',np.argwhere(np.isnan(cs)))
        return trilist, quadlist

    # apply the KKT conditions and get the following work (faster and is used)
    def qpmanual(self,num_itr=15, tol=1e-4, w=20, limit=0.7):
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
            print('interation '+str(i))
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
        if self.plot:
            # plt.figure()
            # for i in range(len(edges)):
            #     coord = self.vertexCoords[edges[i]]
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'b')

            # for i in range(len(self.trilist)):
            #     coord = np.concatenate((vertlist1[self.trilist[i,:]],np.array([vertlist1[self.trilist[i,0]]])))
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'r')

            # for i in range(len(self.quadlist)):
            #     coord = np.concatenate((vertlist1[self.quadlist[i,:]],np.array([vertlist1[self.quadlist[i,0]]])))
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'r')
            # plt.plot(fixed[:,0],fixed[:,1],'k-o')
            # plt.axis('equal')

            # plt.axis('off')
            # plt.savefig("smoothing.pdf", bbox_inches='tight')
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            for i in range(len(edges)):
                coord = self.vertexCoords[edges[i]]
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='k')
            for i in range(len(self.trilist)):
                coord = np.concatenate((vertlist1[self.trilist[i,:]],np.array([vertlist1[self.trilist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')

            for i in range(len(self.quadlist)):
                coord = np.concatenate((vertlist1[self.quadlist[i,:]],np.array([vertlist1[self.quadlist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')
            ax.plot(fixed[:,0],fixed[:,1],fixed[:,2],'k-o')

            # plt.savefig('preliminary_results/optimization.png')
            X = self.vertexCoords[:,0]
            Y = self.vertexCoords[:,1]
            Z = self.vertexCoords[:,2]
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
        if self.plot:
            # plt.figure()
            # # for i in range(len(edges)):
            # #     coord = self.vertexCoords[edges[i]]
            # #     x, y, z = zip(*coord)
            # #     plt.plot(x,y,'b')

            # # for i in range(len(self.trilist)):
            # #     coord = np.concatenate((vertlist1[self.trilist[i,:]],np.array([vertlist1[self.trilist[i,0]]])))
            # #     x, y, z = zip(*coord)
            # #     plt.plot(x,y,'r')

            # for i in range(len(quadlist)):
            #     coord = np.concatenate((vertexCoords[quadlist[i,:]],np.array([vertexCoords[quadlist[i,0]]])))
            #     x, y, z = zip(*coord)
            #     plt.plot(x,y,'b')
            # plt.plot(fixed[:,0],fixed[:,1],'k-o')
            # plt.axis('equal')

            # plt.axis('off')
            # plt.savefig("fully quad.pdf", bbox_inches='tight')
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            for i in range(len(trilist)):
                coord = np.concatenate((vertexCoords[trilist[i]],np.array([vertexCoords[trilist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')

            for i in range(len(quadlist)):
                coord = np.concatenate((vertexCoords[quadlist[i,:]],np.array([vertexCoords[quadlist[i,0]]])))
                x, y, z = zip(*coord)
                ax.plot(x,y,z,color='b')
            ax.plot(fixed[:,0],fixed[:,1],fixed[:,2],'k-o')
            # plt.savefig('preliminary_results/optimization.png')
            X = vertexCoords[:,0]
            Y = vertexCoords[:,1]
            Z = vertexCoords[:,2]
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
        print('fully quad mesh generated')

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
    def saveasvtk(self,vtk_file_name,tri=True):
        '''
        tri: True => save as triangles, False => save as quads
        '''
        print('save as vtk')
        points = self.vertexCoords
        cells = []
        # print(len(points))
        triindex1 = np.array([0,1,2],dtype=np.int32)
        triindex2 = np.array([2,3,0],dtype=np.int32)
        if (tri):
            for i in range(len(self.quadlist)):
                tuple1 = ("triangle",self.quadlist[i,triindex1].reshape(1,3))
                tuple2 = ("triangle",self.quadlist[i,triindex2].reshape(1,3))
                cells.append(tuple1)
                cells.append(tuple2)
        else:
            for i in range(len(self.quadlist)):
                tuple = ("quad",self.quadlist[i,:].reshape(1,4))
                cells.append(tuple)
        meshio.write_points_cells(
            vtk_file_name+".vtk",
            points,
            cells
            )
