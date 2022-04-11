
// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm> // find, sort, max
#include <utility>   // pair
#include <functional>  // 
#include <numeric>  

// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

typedef struct intersection
{
    int triLabel;
    double s;
    vector<double> p;
    intersection() { triLabel=0; p.resize(3,0.0); s=0.0; }
    intersection(int _triLabel, vector<double> _p, double _s) { triLabel=_triLabel; p=_p; s=_s;}
    // Overload operator < so we can use sort function of a std::vector
    bool operator < (const intersection& intersect) const { 
        return (s < intersect.s); }
} Intersection;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct indexanddist
{
    unsigned int branchnum;
    double distance;
    indexanddist() { branchnum=0; distance= 0.0;}
    indexanddist(unsigned int _branchnum, double _distance) { branchnum=_branchnum; distance=_distance;}
}Indexanddist;


typedef struct distandpos
{
    double distance;
    vector<double> position;
    distandpos() { distance=0.0; position.resize(3,0.0);}
    distandpos(double _distance, vector<double> _position) { distance=_distance; position=_position;}
}Distandpos;

typedef struct coordandconnect
{
    vector<vector<double> > coordinate;
    vector<vector<int> > connectivity;
    vector<int> constrained;
    coordandconnect() { coordinate.resize(3, vector<double>(3)); connectivity.resize(3, vector<int>(3));constrained.resize(3,1);}
    coordandconnect(vector<vector<double> > _coordinate, vector<vector<int> > _connectivity, vector<int> _constrained) { coordinate=_coordinate; connectivity=_connectivity; constrained=_constrained;}
}CoordandConnect;

typedef struct projection
{   
    int list;
    int label;
    double distance;
    vector<double> phycoord;
    vector<double> paracoord;
    vector<vector<double> > gradient;
    projection() { list = 0; label=0; distance=0.0; phycoord.resize(3,0.0); paracoord.resize(2,0.0); gradient.resize(3,vector<double>(3));}
    projection(int _list, int _label, double _distance, vector<double> _phycoord, vector<double> _paracoord, vector<vector<double> > _gradient) { list=_list; label=_label; distance=_distance; phycoord=_phycoord; paracoord=_paracoord; gradient=_gradient;}
}Projection;

typedef struct projflatten
{
    vector<double> proj;
    vector<double> gradient;
    projflatten() { proj.resize(3,0.0); gradient.resize(3,0.0);}
    projflatten(vector<double> _proj, vector<double> _gradient) {proj=_proj; gradient=_gradient;}
}ProjFlatten;

typedef struct reconnection
{
    vector<vector<double> > vertlist;
    vector<vector<int> > trilist;
    vector<int> projlist;
    reconnection() { vertlist.resize(3,vector<double>(3,0.0)); trilist.resize(3,vector<int>(3,0)); projlist.resize(3,0);}
    reconnection(vector<vector<double> > _vertlist, vector<vector<int> > _trilist, vector<int> _projlist){ vertlist = _vertlist; trilist = _trilist; projlist = _projlist;}
}Reconnection;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class cLine {
public:

    vector<double> p0,p1,dir;
    cLine();
    cLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~cLine();
    void getDir();
    void getP1();
};

class cOctNode {
public:

    static const int MAX_OCTNODE_OBJECTS  = 200;
    static const int NUM_BRANCHES_OCTNODE = 8;  
    double size;
    int level;
    string nid;
    vector<double> position;
    vector<cOctNode> branches;
    vector<int> data;
    vector<double> low, upp;
    cOctNode();
    cOctNode(int _level, string _nid, vector<double> _position, double _size);
    ~cOctNode();
    bool isLeafNode();
    int numPolys();
    void addPoly(int _indx);
    void addNode(int _level, string _nid, vector<double> _position, double _size);
    void getLowUppVerts();
    bool boxRayIntersect(cLine &ray);
    bool sphereRayIntersect(cLine &ray);
};

class cTri {
public:

    double D;
    int label;
    vector<vector<double> > vertices;
    vector<double> N;
    vector<double> lowVert, uppVert;
    cTri();
    cTri(int _label, vector<vector<double> > _vertices);
    ~cTri();
    bool isInNode(cOctNode &node);
    bool isPointInTri(vector<double> &p);
    bool rayPlaneIntersectPoint(cLine &ray, bool entryOnly);
    bool rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s);
    void getN();
    void getD();
    void getLowerVert();
    void getUpperVert();
};

class cQuad {
public:

    vector<vector<double> > vertices;
    int label;
    double lower, upper;
    cQuad();
    cQuad(int _label, vector<vector<double> > _vertices, double _lower, double _upper);
    ~cQuad();
    vector<double> evaluate(double u, double v);
    vector<double> gradient(double u, double v, vector<double> p0);
    vector<vector<double> > hessian(double u, double v, vector<double> p0);
    double distance(double u, double v, vector<double> p0);
    vector<vector<double> > dfdx(double u, double v, vector<double> p0);
};

class cOctree {
public:

    static const int MAX_OCTREE_LEVELS = 10;
    int branchOffsets[8][3];
    cOctNode root;
    vector<vector<double> > vertexCoords3D;
    vector<vector<int> > polyConnectivity;
    vector<vector<int> > quadlist;
    vector<cTri> polyList;
    vector<cQuad> quadgroup;
    cOctree(vector<vector<double> > _vertexCoords3D, vector<vector<int> > _polyConnectivity, vector<vector<int> > quadlist);
    ~cOctree();    
    double getSizeRoot();
    // ---------------------------------------------------------
    double getavgpolysize(cOctNode &node);
    double getmaxpolysize();
    double getdistance(cOctNode &node, vector<double> projpoint);
    Intersection findmindist(cOctNode &node, vector<double> projpoint);
    vector<Intersection> findmindists(cOctNode &node, vector<vector<double> > projpoints);
    Distandpos findtridist(vector<double> projpoint, cTri &poly);
    CoordandConnect retriangulation(vector<vector<double> > projectedpoint, vector<int> trilist, vector<vector<double> > vertexCoords3D, vector<vector<int> > polyConnectivity);
    vector<vector<int> > countedges();
    vector<vector<int> > countedges(vector<vector<int> > trilist);
    vector<vector<int> > countedges(vector<vector<int> > trilist, vector<vector<double> > vertlist,vector<double> revprojdir);
    vector<vector<double> > computevertexnorm(vector<vector<int> > polyConnectivity, vector<vector<double> > vertexCoords3D);
    // ---------------------------------------------------------
    int numPolys();
    cOctNode* getNodeFromId(string nodeId);
    cOctNode* findBranchById(string nodeId, cOctNode &node);
    set<int> getListPolysToCheck(cLine &ray);    
    vector<double> getPositionRoot();	
    vector<Intersection> findRayIntersect(cLine &ray);    
    // vector<int> findRayIntersects(vector<cLine> &rayList);
    // ---------------------------------------------------------
    // Intersection findRayIntersect(cLine &ray);    
    vector<Intersection> findRayIntersects(vector<cLine> &rayList);
    Intersection findRaymindist(cLine &ray);
    vector<Projection> project(cOctNode &node, vector<vector<double> > projpoints);
    Projection quadproject(vector<double> projpoint, int i);
    ProjFlatten projlist(cOctNode &node, vector<vector<double> > projpoints);
    Reconnection reconnect(cOctNode &node, vector<vector<double> > projpoints, vector<vector<double> > projdirts);
    vector<Intersection> ProjectionOfPoints(cOctNode &node, vector<vector<double> > projpoints, vector<vector<double> > projdirts);
    // ---------------------------------------------------------

    vector<int> findRayIntersectsSorted(vector<cLine> &rayList);		
    vector<cOctNode*> getNodesFromLabel(int polyLabel);	
    vector<cOctNode*> getSortedNodesToCheck(cLine &ray);
    void insertPoly(cOctNode &node, cTri &poly);
    void insertPolys();
    void setupPolyList();
    void setupQuadList();
    void splitNodeAndReallocate(cOctNode &node);
    void findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList);
    void getPolysToCheck(cOctNode &node, cLine &ray, set<int> &intTestPolys);
    void getNodesToCheck(cOctNode &node, cLine &ray, vector<pair<cOctNode*,double> > &nodeList);
};

// Function prototypes
double vectNorm(vector<double> &a);
vector<double> normalize(vector<double> &v);
bool sortNodes(const pair<cOctNode*,double>&i, const pair<cOctNode*,double>&j);
double dotProduct( vector<double> &v1, vector<double> &v2 );
double distBetweenPoints(vector<double> &p1, vector<double> &p2);
string NumberToString( int Number );
vector<double> crossProduct( vector<double> &v1, vector<double> &v2 );
vector<double> vectAdd( vector<double> &a, vector<double> &b);
vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf);
vector<double> vectSubtract( vector<double> &a, vector<double> &b );
vector<double> vectMuldou( vector<double> &a, double sf);
vector<vector<double> > vectMultiply(vector<vector<double> > m1, vector<vector<double> > m2);
vector<vector<double> > vectInverse(vector<vector<double> > &m);
vector<vector<double> > vectneg(vector<vector<double> > m);
vector<double> vectneg(vector<double> &m);
vector<double> midPoint(vector<double> &a, vector<double> &b);
double angleOfVecs(vector<double> &a, vector<double> &b);