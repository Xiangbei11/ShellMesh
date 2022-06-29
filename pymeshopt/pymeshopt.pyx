
from __future__ import print_function
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.string cimport string
from cython.operator cimport dereference as deref, preincrement as inc, predecrement as dec
ctypedef np.float64_t float64
ctypedef np.int32_t int32

cdef extern from "cmeshopt.h":
    cdef cppclass EqualityConstraint:
        vector[int] data
        vector[int] row
        vector[int] col
        vector[int] b
        int sharededges

    cdef cppclass Indices:
        vector[vector[int]] triindex
        vector[int] edgeindex

    cdef cppclass FullyQuad:
        vector[vector[int]] quadlist
        vector[vector[double]] vertlist
        vector[int] fixedvert

    cdef cppclass SplitOutput:
        vector[vector[double]] vertlist
        vector[vector[double]] splitcoord
        vector[vector[int]] trilist
        vector[vector[int]] quadlist
        vector[int] fixedvert
    
    cdef cppclass MergeOutput:
        vector[vector[int]] trilist
        vector[vector[int]] quadlist

    cdef cppclass OptPrep:
        vector[vector[double]] normal 
        vector[double] weight

    cdef cppclass TriInt:
        vector[vector[int]] itris
        vector[vector[int]] jtris
        vector[vector[double]] icoords
        vector[vector[double]] jcoords
        vector[int] ifixed
        vector[int] jfixed

    cdef cppclass Uni:
        vector[vector[double]] vertlist
        vector[vector[int]] quadlist

    cdef cppclass opt:
        opt(vector[vector[double]] _vertexCoords, vector[vector[int]] _triConnectivity, vector[vector[int]] _quadConnectivity, double _w1, double _w2, double _w3)
        vector[double] calculate_aspect_ratio()
        vector[vector[double]] calculate_internal_angles()
        vector[vector[double]] buildupmatrix()
        vector[vector[double]] buildbounds(vector[int] fixedvert)
        EqualityConstraint equcon()
        vector[double] mergetriangles()
        vector[vector[int]] inequncon()
        vector[vector[double]] mergingbounds(vector[int] fixedvert)
        vector[vector[int]] countedges()
        Indices counttrisharededges()
        FullyQuad fullyquadmesh(vector[int] fixed)
        OptPrep computenormal()
        vector[vector[vector[double]]] computecs()
        SplitOutput splitupdatemesh(vector[int] optres, vector[int] fixed)
        MergeOutput mergeupdatemesh(vector[int] optres)
        vector[vector[int]] connectTFI(vector[vector[double]] boundaries, vector[vector[double]] vertices, vector[int] num_b, vector[int] num_v, vector[int] index_b, vector[int] index_v)
        double averagesize()
        vector[double] computeradius()
        TriInt member_intersection(vector[vector[int]] itris, vector[vector[int]] jtris, vector[vector[double]] icoords, vector[vector[double]] jcoords, vector[int] ifixed, vector[int] jfixed)
        vector[int] removeduplicatetris()
        Uni mergeduppts(double tol)

# generally speaking, the only thing that is done in pyx file is transformation of data types bewtween python and cpp variables, in our code, no calculation is done in here
cdef class Pymeshopt:
    cdef opt *thisopt

    def __cinit__(self, np.ndarray[float,ndim=2] _vertexCoords, np.ndarray[int,ndim=2] _triConnectivity, np.ndarray[int,ndim=2] _quadConnectivity, float _w1, float _w2, float _w3):

        cdef int i, j
        cdef vector[double] coords
        cdef vector[vector[double]] vertexCoords
        vertexCoords.reserve(_vertexCoords.shape[0])
        coords.resize(3)
        for i in range(_vertexCoords.shape[0]):
            for j in range(3):
                coords[j] = _vertexCoords[i,j]
            vertexCoords.push_back(coords)

        cdef vector[int] triconnect
        cdef vector[vector[int]] triConnectivity
        # tri connectivity
        triConnectivity.reserve(_triConnectivity.shape[0])
        triconnect.resize(3)
        for i in range(_triConnectivity.shape[0]):
            for j in range(3):
                triconnect[j] = _triConnectivity[i,j]
            triConnectivity.push_back(triconnect)
        
        cdef vector[int] quadconnect
        cdef vector[vector[int]] quadConnectivity
        quadConnectivity.reserve(_quadConnectivity.shape[0])
        quadconnect.resize(4)
        for i in range(_quadConnectivity.shape[0]):
            for j in range(4):
                quadconnect[j] = _quadConnectivity[i,j]
            quadConnectivity.push_back(quadconnect)
        cdef double w1 = _w1
        cdef double w2 = _w2
        cdef double w3 = _w3
        self.thisopt = new opt(vertexCoords, triConnectivity, quadConnectivity, w1, w2, w3)
    
    def __dealloc__(self):
        del self.thisopt

    def calculate_aspect_ratio(self):
        cdef int i
        cdef vector[double] aspect_ratio = self.thisopt.calculate_aspect_ratio()
        py_aspect_ratio = np.zeros((aspect_ratio.size(),1))
        for i in range(aspect_ratio.size()):
            py_aspect_ratio[i] = aspect_ratio[i]
        return py_aspect_ratio

    def calculate_internal_angles(self):
        cdef int i, j
        cdef vector[vector[double]] internal_angles = self.thisopt.calculate_internal_angles()
        py_internal_angles = np.zeros((internal_angles.size(),4))
        for i in range(internal_angles.size()):
            for j in range(4):
                py_internal_angles[i,j]= internal_angles[i][j]
        return py_internal_angles

    def create_splitting_matrix(self):
        cdef int i, j
        cdef vector[vector[double]] matrix = self.thisopt.buildupmatrix()
        pymatrix = np.zeros((matrix.size(),11))
        for i in range(matrix.size()):
            for j in range(11):
                pymatrix[i,j]=matrix[i][j]
        return pymatrix

    def createbounds(self, np.ndarray[int,ndim=1] _fixedvert):
        cdef int i,j
        cdef vector[int] fixedvert
        fixedvert.resize(_fixedvert.shape[0])
        for i in range(_fixedvert.shape[0]):
            fixedvert[i] = _fixedvert[i]
        cdef vector[vector[double]] bounds = self.thisopt.buildbounds(fixedvert)
        pybounds = np.zeros((bounds.size(),bounds[0].size()))
        for i in range(bounds.size()):
            for j in range(bounds[0].size()):
                pybounds[i,j] = bounds[i][j]
        
        return pybounds

    def equcon(self):
        cdef EqualityConstraint equalityconstraint = self.thisopt.equcon()
        cdef equalitycon = EqualityCon(equalityconstraint.data.size(),equalityconstraint.b.size())
        equalitycon.SetValues(equalityconstraint.data,equalityconstraint.row,equalityconstraint.col,equalityconstraint.b,equalityconstraint.sharededges)
        return equalitycon

    def create_merging_vector(self):
        cdef int i
        cdef vector[double] cvec = self.thisopt.mergetriangles()
        pyvec = np.zeros(cvec.size())
        for i in range(cvec.size()):
            pyvec[i] = cvec[i]
        return pyvec
    
    def createboundsmerge(self, np.ndarray[int,ndim=1] _fixedvert):
        cdef int i,j
        cdef vector[int] fixedvert
        fixedvert.resize(_fixedvert.shape[0])
        for i in range(_fixedvert.shape[0]):
            fixedvert[i] = _fixedvert[i]
        cdef vector[vector[double]] bounds = self.thisopt.mergingbounds(fixedvert)
        pybounds = np.zeros((bounds.size(),bounds[0].size()))
        for i in range(bounds.size()):
            for j in range(bounds[0].size()):
                pybounds[i,j] = bounds[i][j]
        return pybounds
    
    def ubcon(self):
        cdef int i,j
        cdef vector[vector[int]] ubcon = self.thisopt.inequncon()
        pyubcon = np.zeros((ubcon.size(),ubcon[0].size()))
        for i in range(ubcon.size()):
            for j in range(ubcon[0].size()):
                pyubcon[i,j] = ubcon[i][j]
        return pyubcon

    def uniedges(self):
        cdef int i, j
        cdef vector[vector[int]] edges = self.thisopt.countedges()
        pyedges = np.zeros((edges.size(),edges[0].size()))
        for i in range(edges.size()):
            for j in range(edges[0].size()):
                pyedges[i,j] = edges[i][j]
        return pyedges

    def counttrisharededges(self):
        cdef int inequncon
        cdef Indices indices = self.thisopt.counttrisharededges()
        cdef pyindices = PyIndices(indices.triindex.size())
        pyindices.SetValues(indices.triindex, indices.edgeindex)
        return pyindices
        
    def fullyquadmesh(self, np.ndarray[int,ndim=1] _fixed):
        cdef int i
        cdef vector[int] fixed
        fixed.resize(_fixed.shape[0])
        for i in range(_fixed.shape[0]):
            fixed[i] = _fixed[i]
        cdef FullyQuad fullyquad = self.thisopt.fullyquadmesh(fixed)
        cdef pyfullyquad = PyFullyQuad(fullyquad.quadlist.size(),fullyquad.vertlist.size(),fullyquad.fixedvert.size())
        pyfullyquad.SetValues(fullyquad.quadlist, fullyquad.vertlist, fullyquad.fixedvert)
        return pyfullyquad

    def computenormal(self):
        cdef int i, j
        cdef OptPrep optprep = self.thisopt.computenormal()
        cdef pyoptprep = PyOptPrep(optprep.normal.size())
        pyoptprep.SetValues(optprep.normal,optprep.weight)
        return pyoptprep

    def computecs(self):
        cdef int i, j, k
        cdef vector[vector[vector[double]]] cs = self.thisopt.computecs()
        cdef pycs = np.zeros((2,cs[0].size(),3),dtype=float)
        for i in range(2):
            for j in range(cs[0].size()):
                for k in range(3):
                    pycs[i][j][k] = cs[i][j][k]
        return pycs

    def splitupdatemesh(self, np.ndarray[int,ndim=1] _optres, np.ndarray[int,ndim=1] _fixed):
        cdef vector[int] optres
        cdef int i
        optres.resize(_optres.shape[0])
        for i in range(_optres.shape[0]):
            optres[i] = _optres[i]
        cdef vector[int] fixed
        fixed.resize(_fixed.shape[0])
        for i in range(_fixed.shape[0]):
            fixed[i] = _fixed[i]
        cdef SplitOutput splitoutput = self.thisopt.splitupdatemesh(optres,fixed)
        cdef pysplitoutput = PySplitOutput(splitoutput.vertlist.size(),splitoutput.splitcoord.size(),splitoutput.trilist.size(),splitoutput.quadlist.size(),splitoutput.fixedvert.size())
        pysplitoutput.SetValues(splitoutput.vertlist,splitoutput.splitcoord,splitoutput.trilist,splitoutput.quadlist,splitoutput.fixedvert)
        return pysplitoutput
    
    def mergeupdatemesh(self, np.ndarray[int,ndim=1] _optres):
        cdef vector[int] optres
        cdef int i
        optres.resize(_optres.shape[0])
        for i in range(_optres.shape[0]):
            optres[i] = _optres[i]
        cdef MergeOutput mergeoutput = self.thisopt.mergeupdatemesh(optres)
        cdef pymergeoutput = PyMergeOutput(mergeoutput.trilist.size(),mergeoutput.quadlist.size())
        pymergeoutput.SetValues(mergeoutput.trilist,mergeoutput.quadlist)
        return pymergeoutput

    def connectTFI(self, np.ndarray[float,ndim=2] _boundaries, np.ndarray[float,ndim=2] _vertices, np.ndarray[int,ndim=1] _num_b, np.ndarray[int,ndim=1] _num_v, np.ndarray[int,ndim=1] _index_b, np.ndarray[int,ndim=1] _index_v):
        cdef int i, j
        cdef vector[double] boundary
        cdef vector[vector[double]] boundaries
        boundaries.reserve(_boundaries.shape[0])
        boundary.resize(3)
        for i in range(_boundaries.shape[0]):
            for j in range(3):
                boundary[j] = _boundaries[i,j]
            boundaries.push_back(boundary)

        cdef vector[double] vertex
        cdef vector[vector[double]] vertices
        vertices.reserve(_vertices.shape[0])
        vertex.resize(3)
        for i in range(_vertices.shape[0]):
            for j in range(3):
                vertex[j] = _vertices[i,j]
            vertices.push_back(vertex)

        cdef vector[int] num_b
        num_b.resize(_num_b.shape[0])
        for i in range(_num_b.shape[0]):
            num_b[i] = _num_b[i]
        
        cdef vector[int] num_v
        num_v.resize(_num_v.shape[0])
        for i in range(_num_v.shape[0]):
            num_v[i] = _num_v[i]
        
        cdef vector[int] index_b
        index_b.resize(_index_b.shape[0])
        for i in range(_index_b.shape[0]):
            index_b[i] = _index_b[i]

        cdef vector[int] index_v
        index_v.resize(_index_v.shape[0])
        for i in range(_index_v.shape[0]):
            index_v[i] = _index_v[i]
        cdef vector[vector[int]] tris = self.thisopt.connectTFI(boundaries,vertices,num_b,num_v,index_b,index_v)
        pytris = np.zeros((tris.size(),3))
        for i in range(tris.size()):
            for j in range(3):
                pytris[i][j] = tris[i][j]
        return pytris
    
    def averagesize(self):
        cdef double distance = self.thisopt.averagesize()
        return distance

    def computeradius(self):
        cdef vector[double] radius = self.thisopt.computeradius()
        pyradius = np.zeros(radius.size())
        for i in range(radius.size()):
            pyradius[i] = radius[i]
        return pyradius
    
    def memberintersection(self, np.ndarray[int,ndim=2] _itrilist, np.ndarray[int,ndim=2] _jtrilist, np.ndarray[float,ndim=2] _icoords, np.ndarray[float,ndim=2] _jcoords, np.ndarray[int,ndim=1] _ifixed, np.ndarray[int,ndim=1] _jfixed):
        cdef int i, j
        # itrilist
        cdef vector[int] itri
        cdef vector[vector[int]] itrilist
        itrilist.reserve(_itrilist.shape[0])
        itri.resize(3)
        for i in range(_itrilist.shape[0]):
            for j in range(3):
                itri[j] = _itrilist[i,j]
            itrilist.push_back(itri)
        # jtrilist
        cdef vector[int] jtri
        cdef vector[vector[int]] jtrilist
        jtrilist.reserve(_jtrilist.shape[0])
        jtri.resize(3)
        for i in range(_jtrilist.shape[0]):
            for j in range(3):
                jtri[j] = _jtrilist[i,j]
            jtrilist.push_back(jtri)
        # icoords
        cdef vector[double] icoord
        cdef vector[vector[double]] icoords
        icoords.reserve(_icoords.shape[0])
        icoord.resize(3)
        for i in range(_icoords.shape[0]):
            for j in range(3):
                icoord[j] = _icoords[i,j]
            icoords.push_back(icoord)
        # jcoords
        cdef vector[double] jcoord
        cdef vector[vector[double]] jcoords
        jcoords.reserve(_jcoords.shape[0])
        jcoord.resize(3)
        for i in range(_jcoords.shape[0]):
            for j in range(3):
                jcoord[j] = _jcoords[i,j]
            jcoords.push_back(jcoord)
        # ifixed
        cdef vector[int] ifixed
        ifixed.resize(_ifixed.shape[0])
        for i in range(_ifixed.shape[0]):
            ifixed[i] = _ifixed[i]
        # jfixed
        cdef vector[int] jfixed
        jfixed.resize(_jfixed.shape[0])
        for i in range(_jfixed.shape[0]):
            jfixed[i] = _jfixed[i]
        
        cdef TriInt triint = self.thisopt.member_intersection(itrilist,jtrilist,icoords,jcoords,ifixed,jfixed)
        pytriint = PyTriInt(triint.itris.size(),triint.jtris.size(),triint.icoords.size(),triint.jcoords.size(),triint.ifixed.size(),triint.jfixed.size())
        pytriint.SetValues(triint.itris,triint.jtris,triint.icoords,triint.jcoords,triint.ifixed,triint.jfixed)
        return pytriint
            
    def pyremoveduplicatetris(self):
        cdef vector[int] output = self.thisopt.removeduplicatetris()
        pyoutput = np.zeros(output.size(),dtype=np.int32)
        for i in range(output.size()):
            pyoutput[i] = output[i]
        return pyoutput
        
    def pymergeduppts(self):
        cdef double tol = 1e-5
        cdef Uni uni = self.thisopt.mergeduppts(tol)
        pyuni = PyUni(uni.vertlist.size(),uni.quadlist.size())
        pyuni.SetValues(uni.vertlist,uni.quadlist)
        return pyuni



# ------------------------------------------------------------------------------------------------------------------------------------------------
cdef class EqualityCon:
    cdef public np.ndarray data
    cdef public np.ndarray row
    cdef public np.ndarray col
    cdef public np.ndarray b
    cdef public int sharededges
    cdef public int len1
    cdef public int len2
    def __init__(self, len1, len2):
        self.data = np.zeros(len1,dtype=int)
        self.row = np.zeros(len1,dtype=int)
        self.col = np.zeros(len1,dtype=int)
        self.b = np.zeros(len2,dtype=int)
        self.sharededges = 0
        self.len1 = len1
        self.len2 = len2
    def SetValues(self, vector[int] data, vector[int] row, vector[int] col, vector[double] b, int sharededges):
        cdef int i, j
        for i in range(self.len1):
            self.data[i] = data[i]
            self.row[i] = row[i]
            self.col[i] = col[i]
        for j in range(self.len2):
            self.b[j] = b[j]
        self.sharededges = sharededges

cdef class PyIndices:
    cdef public np.ndarray triindex
    cdef public np.ndarray edgeindex
    cdef public int len
    def __init__(self, len):
        self.triindex = np.zeros((len,2),dtype=int)
        self.edgeindex = np.zeros(len,dtype=int)
        self.len = len
    def SetValues(self, vector[vector[int]] triindex, vector[int] edgeindex):
        cdef int i, j
        for i in range(self.len):
            for j in range(2):
                self.triindex[i][j] = triindex[i][j]
            self.edgeindex[i] = edgeindex[i]

cdef class PyFullyQuad:
    cdef public np.ndarray quadlist
    cdef public np.ndarray vertlist
    cdef public np.ndarray fixedvert
    cdef public int len1
    cdef public int len2
    cdef public int len3
    def __init__(self, len1, len2, len3):
        self.quadlist = np.zeros((len1,4),dtype=int)
        self.vertlist = np.zeros((len2,3),dtype=float)
        self.fixedvert = np.zeros(len3,dtype=int)
        self.len1 = len1
        self.len2 = len2
        self.len3 = len3
    def SetValues(self, vector[vector[int]] quadlist, vector[vector[double]] vertlist, vector[int] fixedvert):
        cdef int i, j
        for i in range(self.len1):
            for j in range(4):
                self.quadlist[i][j] = quadlist[i][j]
        for i in range(self.len2):
            for j in range(3):
                self.vertlist[i][j] = vertlist[i][j]
        for i in range(self.len3):
            self.fixedvert[i] = fixedvert[i]

cdef class PySplitOutput:
    cdef public np.ndarray vertlist
    cdef public np.ndarray splitcoord
    cdef public np.ndarray trilist
    cdef public np.ndarray quadlist
    cdef public np.ndarray fixedvert
    cdef public int len1
    cdef public int len2
    cdef public int len3
    cdef public int len4
    cdef public int len5
    def __init__(self,len1,len2,len3,len4,len5):
        self.vertlist = np.zeros((len1,3),dtype=float)
        self.splitcoord = np.zeros((len2,3),dtype=float)
        self.trilist = np.zeros((len3,3),dtype=int)
        self.quadlist = np.zeros((len4,4),dtype=int)
        self.fixedvert = np.zeros(len5,dtype=int)
        self.len1 = len1    
        self.len2 = len2
        self.len3 = len3
        self.len4 = len4
        self.len5 = len5

    def SetValues(self, vector[vector[double]] vertlist, vector[vector[double]] splitcoord, vector[vector[int]] trilist, vector[vector[int]] quadlist, vector[int] fixedvert):
        cdef int i, j
        for i in range(self.len1):
            for j in range(3):
                self.vertlist[i][j] = vertlist[i][j]
        for i in range(self.len2):
            for j in range(3):
                self.splitcoord[i][j] = splitcoord[i][j]
        for i in range(self.len3):
            for j in range(3):
                self.trilist[i][j] = trilist[i][j]
        for i in range(self.len4):
            for j in range(4):
                self.quadlist[i][j] = quadlist[i][j]
        for i in range(self.len5):
            self.fixedvert[i] = fixedvert[i]


cdef class PyMergeOutput:
    cdef public np.ndarray trilist
    cdef public np.ndarray quadlist
    cdef public int len1
    cdef public int len2
    def __init__(self,len1,len2):
        self.trilist = np.zeros((len1,3),dtype=int)
        self.quadlist = np.zeros((len2,4),dtype=int)
        self.len1 = len1
        self.len2 = len2
    def SetValues(self, vector[vector[int]] trilist, vector[vector[int]] quadlist):
        cdef int i, j
        for i in range(self.len1):
            for j in range(3):
                self.trilist[i][j] = trilist[i][j]
        for i in range(self.len2):
            for j in range(4):
                self.quadlist[i][j] = quadlist[i][j]

cdef class PyOptPrep:
    cdef public np.ndarray normal
    cdef public np.ndarray weight
    def __init__(self, len):
        # print(len)
        self.normal = np.zeros((len,3),dtype=float)
        self.weight = np.zeros(len,dtype=float)
    def SetValues(self, vector[vector[double]] normal, vector[double] weight):
        for i in range(normal.size()):
            for j in range(3):
                self.normal[i][j] = normal[i][j]
        # print('weight.size:',weight.size())
        # print('self.weight.size:',self.weight.shape[0])
        # print('self.weight',self.weight)
        # print('self.normal.size:',self.normal.shape[0])
        # print('self.normal',self.weight)
        for i in range(weight.size()):
            self.weight[i] = weight[i]     

cdef class PyTriInt:
    cdef public np.ndarray itris
    cdef public np.ndarray jtris
    cdef public np.ndarray icoords
    cdef public np.ndarray jcoords
    cdef public np.ndarray ifixed
    cdef public np.ndarray jfixed
    cdef public int len1
    cdef public int len2
    cdef public int len3
    cdef public int len4
    cdef public int len5
    cdef public int len6
    def __init__(self,len1,len2,len3,len4,len5,len6):
        self.itris = np.zeros((len1,3),dtype=int)
        self.jtris = np.zeros((len2,3),dtype=int)
        self.icoords = np.zeros((len3,3),dtype=float)
        self.jcoords = np.zeros((len4,3),dtype=float)
        self.ifixed = np.zeros(len5,dtype=int)
        self.jfixed = np.zeros(len6,dtype=int)
        self.len1 = len1
        self.len2 = len2
        self.len3 = len3
        self.len4 = len4
        self.len5 = len5
        self.len6 = len6

    def SetValues(self, vector[vector[int]] itris, vector[vector[int]] jtris, vector[vector[double]] icoords, vector[vector[double]] jcoords, vector[int] ifixed, vector[int] jfixed):
        for i in range(self.len1):
            for j in range(3):
                self.itris[i][j] = itris[i][j]
        for i in range(self.len2):
            for j in range(3):
                self.jtris[i][j] = jtris[i][j]
        for i in range(self.len3):
            for j in range(3):
                self.icoords[i][j] = icoords[i][j]
        for i in range(self.len4):
            for j in range(3):
                self.jcoords[i][j] = jcoords[i][j]
        for i in range(self.len5):
            self.ifixed[i] = ifixed[i]
        for i in range(self.len6):
            self.jfixed[i] = jfixed[i]
        

cdef class PyUni:
    cdef public np.ndarray vertlist
    cdef public np.ndarray quadlist
    cdef public int len1
    cdef public int len2
    def __init__(self,len1,len2):
        self.vertlist = np.zeros((len1,3),dtype=float)
        self.quadlist = np.zeros((len2,4),dtype=int)
        self.len1 = len1
        self.len2 = len2
    def SetValues(self, vector[vector[double]] vertlist, vector[vector[int]] quadlist):
        cdef int i, j
        for i in range(self.len1):
            for j in range(3):
                self.vertlist[i][j] = vertlist[i][j]
        for i in range(self.len2):
            for j in range(4):
                self.quadlist[i][j] = quadlist[i][j]
