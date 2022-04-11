
// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include "cOctree.h"

// ------------------------------------------------------

cLine::cLine() 
{
    // Default constructor for cLine 
    // Default line is unit vector along x-axis
    p0.resize(3,0.0); p1.resize(3,0.0); dir.resize(3,0.0);
    p1[0]=1.0; dir[0]=1.0;
}

cLine::cLine(vector<double> &_p0, vector<double> &p1_dir, int isP1orDir)
{
    // cLine constructor with p0 and p1 or dir
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = _p0;
    if (isP1orDir==0) {
        p1 = p1_dir;
        getDir(); }
    else if (isP1orDir==1) {
        dir = p1_dir;
        getP1(); }
}

// cLine destructor
cLine::~cLine() {}

void cLine::getDir()
{
    // Get unit vector defining direction of cLine
    vector<double> p0p1(3); double dmag=0.0;
    for (unsigned int i=0; i<3; i++) {
        p0p1[i] = p1[i]-p0[i];
        dmag   += pow(p0p1[i],2.0); }
    dmag = sqrt(dmag);
    dir  = p0p1;
    for (vector<double>::iterator it=dir.begin(); it!=dir.end(); ++it)
        *it /= dmag; 
}

void cLine::getP1()
{
    // Get a point on the cLine, p1, located 1.0 units away from the origin, p0
    vector<double> p1(3);
    for (unsigned int i=0; i<3; i++)
        p1[i] = p0[i]+dir[i];
}

// ------------------------------------------------------

cTri::cTri()
{
    // Default cTri constructor
    label = 0;
    vertices.resize(3);
    for (vector<vector<double> >::iterator it=vertices.begin(); it!=vertices.end(); ++it)
        (*it).resize(3,0.0);
    vertices[1][0]=1.0; vertices[2][1]=1.0;
    getN();
    getD();
    getLowerVert();
    getUpperVert();    
}

cTri::cTri(int _label, vector<vector<double> > _vertices)
{
    // cTri constructor with label and vertices
    label    = _label;
    vertices = _vertices;
    getN();
    getD();
    getLowerVert();
    getUpperVert();
}

// cTri destructor
cTri::~cTri() {
    //cout << "Destroying cTri" << endl;
}

void cTri::getN()
{
    // Get cTri face normal
    vector<vector<double> > v = vertices;
    vector<double> v1(3),v2(3),v3;
    for (unsigned int i=0; i<3; i++) {
        v1[i] = v[0][i] - v[1][i];
        v2[i] = v[0][i] - v[2][i]; }
    v3 = crossProduct(v1,v2);
    double v3mag = sqrt(dotProduct(v3,v3));
    N = v3;
    for (vector<double>::iterator it=N.begin(); it!=N.end(); ++it)
        *it /= v3mag;
}

void cTri::getD()
{
    // Perp distance of cTri face from origin which, along with the face normal,
    // defines the plane of the cTri face
    D = dotProduct(vertices[0],N);
}

void cTri::getLowerVert()
{
    // Lower vertices of cTri bounding box
    lowVert.resize(3,1.0e+30);
    for (unsigned int j=0; j<3; j++) {
        for (unsigned int i=0; i<3; i++) {
            if (vertices[i][j] < lowVert[j]) 
            {
                lowVert[j] = vertices[i][j]; 
            }
        }
    } 
}

void cTri::getUpperVert()
{
    // Upper vertices of cTri bounding box
    uppVert.resize(3,-1.0e+30);
    for (unsigned int j=0; j<3; j++) {
        for (unsigned int i=0; i<3; i++) {
            if (vertices[i][j] > uppVert[j]) 
            {
                uppVert[j] = vertices[i][j]; 
            }
        }
    } 
} 

bool cTri::isInNode(cOctNode &node)
{
    // Tests if bounding box of cTri is inside of or overlapping the given cOctNode
    // This is a simple test and even if bounding box is found to be inside the
    // cOctNode, the cTri itself may not be
    if (lowVert[0] > node.upp[0]) return false;
    if (lowVert[1] > node.upp[1]) return false;
    if (lowVert[2] > node.upp[2]) return false;
    if (uppVert[0] < node.low[0]) return false;
    if (uppVert[1] < node.low[1]) return false;
    if (uppVert[2] < node.low[2]) return false;
    return true;
}

bool cTri::isPointInTri(vector<double> &p)
{
    // Determines if point p is within the cTri by computing and
    // testing the barycentric coordinates (u, v, w) of p

    // Find Barycentric coordinates of point (u,v,w)
    vector<double> v0 = vectSubtract(vertices[1],vertices[0]);
    vector<double> v1 = vectSubtract(vertices[2],vertices[0]);
    vector<double> v2 = vectSubtract(p,vertices[0]);
    double d00 = dotProduct(v0, v0);
    double d01 = dotProduct(v0, v1);
    double d11 = dotProduct(v1, v1);
    double d20 = dotProduct(v2, v0);
    double d21 = dotProduct(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    // Use Barycentric coordinates to work out if point lies within the cTri element
    double tol = 0.0;	
    return ((v>=tol) && (w>=tol) && (u>=tol));
}

bool cTri::rayPlaneIntersectPoint(cLine &ray, bool entryOnly=false)
{
    // Tests if ray intersects with the cTri face
    // NOTE: Provide option to just check for entry intersections, not both 
    //       entries/exits. This should cut down checking somewhat.
    double tol  = 1.0e-06;
    double sDen = dotProduct(ray.dir,N);
    if ((entryOnly && sDen>tol) || (!entryOnly && fabs(sDen)>tol))
    {
        double sNum = D - dotProduct(ray.p0,N);
        double s = sNum / sDen;
        vector<double> p = vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    } 
    return false;
}

bool cTri::rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s)
{
    // Tests if ray intersects with the cTri face 
    // Returns the coordinates of the intersection point and the distance, s, from
    // the origin of the ray   
    double tol  = 1.0e-06;
    double sDen = dotProduct(ray.dir,N);
    if (fabs(sDen)> tol) // Normals cannot be perpendicular such that dot product equals 0
    {
        double sNum = D - dotProduct(ray.p0,N);
        s = sNum / sDen;
        p = vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    } 
    return false;
}
// ------------------------------------------------------
// Ning's code start here: this is the definition of cquad
cQuad::cQuad()
{
    label = 0;
    vertices.resize(4);
    for (vector<vector<double> >:: iterator it=vertices.begin(); it != vertices.end(); it++){
        (*it).resize(3,0.0);
    }
    lower = 1.0;
    upper = 1.0;
}

cQuad::cQuad(int _label, vector<vector<double> > _vertices, double _lower, double _upper)
{
    label = _label;
    vertices = _vertices;
    lower = _lower;
    upper = _upper;
}


cQuad::~cQuad(){

}

// calculate the physical coordinates p given the parametric coordinates u and v
vector<double> cQuad::evaluate(double u, double v)
{   
    unsigned int i;
    vector<double> p(3,0.0);
    for (i=0;i<3;i++){
        p[i] = (1-u)*(1-v)*vertices[0][i] + u*(1-v)*vertices[1][i] + u*v*vertices[2][i] + (1-u)*v*vertices[3][i];
    }
    return p;
}

// calculate the gradient pu, pv, puv, and projection direction between p and p0 given the parametric coordinates u and v
vector<double> cQuad::gradient(double u, double v, vector<double> p0)
{
    unsigned int i;
    vector<double> gradient(2,0.0);
    vector<double> pu(3,0.0), pv(3,0.0), puv(3,0.0), f(3,0.0);
    vector<double> p = evaluate(u,v);
    for (i=0;i<3;i++){
        pu[i] = (-1)*(1-v)*vertices[0][i] + 1*(1-v)*vertices[1][i] + 1*v*vertices[2][i] + (-1)*v*vertices[3][i];
        pv[i] = (1-u)*(-1)*vertices[0][i] + u*(-1)*vertices[1][i] + u*1*vertices[2][i] + (1-u)*1*vertices[3][i];
        puv[i] = 1*vertices[0][i] - 1*vertices[1][i] + 1*vertices[2][i] - 1*vertices[3][i];
    }
    f = vectSubtract(p,p0);
    gradient[0] = 2*dotProduct(f,pu);
    gradient[1] = 2*dotProduct(f,pv);
    return gradient;
}
// calculate the hessian matrix 
vector<vector<double> > cQuad::hessian(double u, double v, vector<double> p0)
{
    unsigned int i;
    vector<vector<double> > hessian(2,vector<double>(2,0.0));
    vector<double> pu(3,0.0), pv(3,0.0), puv(3,0.0), f(3,0.0);
    vector<double> p = evaluate(u,v);
    for (i=0;i<3;i++){
        pu[i] = (-1)*(1-v)*vertices[0][i] + 1*(1-v)*vertices[1][i] + 1*v*vertices[2][i] + (-1)*v*vertices[3][i];
        pv[i] = (1-u)*(-1)*vertices[0][i] + u*(-1)*vertices[1][i] + u*1*vertices[2][i] + (1-u)*1*vertices[3][i];
        puv[i] = 1*vertices[0][i] - 1*vertices[1][i] + 1*vertices[2][i] - 1*vertices[3][i];
    }
    f = vectSubtract(p,p0);
    hessian[0][0] = 2*dotProduct(pu,pu);
    hessian[0][1] = 2*dotProduct(pu,pv) + 2*dotProduct(f,puv);
    hessian[1][0] = 2*dotProduct(pu,pv) + 2*dotProduct(f,puv);
    hessian[1][1] = 2*dotProduct(pv,pv);
    return hessian;
}

// distance between p0 and p
double cQuad::distance(double u, double v, vector<double> p0)
{
    vector<double> p = evaluate(u,v);
    vector<double> f = vectSubtract(p,p0);
    double distance = sqrt(dotProduct(f,f));
    return distance;
}

// adjoint method to compute the devrivative of f wrt x
vector<vector<double> > cQuad::dfdx(double u, double v, vector<double> p0)
{
    unsigned int i;
    vector<double> pu(3,0.0);
    vector<double> pv(3,0.0);
    vector<double> puv(3,0.0);
    vector<double> p = evaluate(u,v);
    vector<double> f = vectSubtract(p,p0);
    for (i=0;i<3;i++){
        pu[i] = (-1)*(1-v)*vertices[0][i] + 1*(1-v)*vertices[1][i] + 1*v*vertices[2][i] + (-1)*v*vertices[3][i];
        pv[i] = (1-u)*(-1)*vertices[0][i] + u*(-1)*vertices[1][i] + u*1*vertices[2][i] + (1-u)*1*vertices[3][i];
        puv[i] = 1*vertices[0][i] - 1*vertices[1][i] + 1*vertices[2][i] - 1*vertices[3][i];
    }
    vector<vector<double> > dRdx(2,vector<double>(3,0.0));
    vector<vector<double> > dRdy(2,vector<double>(2,0.0));
    vector<vector<double> > dFdy(3,vector<double>(2,0.0));
    vector<vector<double> > dydx(2,vector<double>(3,0.0));
    vector<vector<double> > dfdx(3,vector<double>(3,0.0));
    vector<vector<double> > invdRdy(2,vector<double>(2,0.0));
    if (((u<1e-8)||(u>1-1e-8))&&((v<1e-8)||(v>1-1e-8))){// projection locates on vertex
        
    }else if ((u<1e-8)||(u>1-1e-8)){// projection locates on u edge
        dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
        dRdy[0][0] = 1; dRdy[1][0] = dotProduct(puv,f) + dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
        dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
        invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
        dfdx = vectMultiply(dFdy,dydx);
    }else if ((v<1e-8)||(v>1-1e-8)){
        dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2];
        dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(puv,f) + dotProduct(pu,pv); dRdy[1][1] = 1;
        dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
        invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
        dfdx = vectMultiply(dFdy,dydx);
    }else{
        dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2]; dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
        dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(puv,f) + dotProduct(pu,pv);
        dRdy[1][0] = dotProduct(puv,f) + dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
        dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
        invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
        dfdx = vectMultiply(dFdy,dydx);
    }
    return dfdx;
}

// Ning's code ends here
// ------------------------------------------------------

cOctNode::cOctNode()
{  
    // Default octNode constructor
    level = 0;
    nid   = "";
    size  = 1.0;
    position.resize(3,0.0);
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);
}

cOctNode::cOctNode(int _level, string _nid, vector<double> _position, double _size)
{
    // octNode constructor with level, node id (nid), position and size
    level    = _level;
    nid      = _nid;
    position = _position;
    size     = _size;
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);    
}

// octNode destructor
cOctNode::~cOctNode() {
    //cout << "Calling destructor for cOctnode " << nid << endl;
}

bool cOctNode::isLeafNode() 
{ 
    // Checks if cOctNode is a leaf node by counting the number of branches. A 
    // leaf node has no branches
    return branches.size()==0; 
}

void cOctNode::getLowUppVerts() 
{
    // Get coordinates of the lower and upper vertices of the cOctNode
    low.resize(3);
    upp.resize(3);
    double halfSize = size/2.0;
    for (unsigned int i=0; i<3; i++) {
        low[i] = position[i] - halfSize;
        upp[i] = position[i] + halfSize;
    }
}

void cOctNode::addPoly(int _indx) { data.push_back(_indx); }

int cOctNode::numPolys() { return (int)(data.size()); }

void cOctNode::addNode(int _level, string _nid, vector<double> _position, double _size)
{
    branches.push_back(cOctNode(_level,_nid,_position,_size));
}

bool cOctNode::sphereRayIntersect(cLine &ray)
{
    // Quick test for determining if a ray is *likely* to intersect a given node

    // Radius of sphere that contains node
    double radius = distBetweenPoints(low,position);
    
    // Project centre of sphere (node.position) onto ray
    vector<double> oc = vectSubtract(position, ray.p0);
    double s = dotProduct(oc,ray.dir);
    vector<double> projpnt = vectAdd(ray.p0, ray.dir, s);
    double dist = distBetweenPoints(projpnt,position);
    
    // If distance between spherical centre and projected point is 
    // less than the radius of the sphere, then an intersection is
    // *possible*
    return (dist<=radius);
}

bool cOctNode::boxRayIntersect(cLine &ray)
{
    // An accurate test for determining if a ray will intersect a given node.
    // Tests for intersections between the ray and all 6 faces of the node.
    
    vector<double> p; double D, sDen, sNum, s, tol = 1.0e-06; int i, j;
    
    for (unsigned int faceId=1; faceId<=6; faceId++)
    {
        // Get D (distance of plane to origin) and N (face normal) of node face
        vector<double> N(3,0.0);
        switch(faceId) {
            case 1: {
                D = -low[0]; N[0] = -1.0; break; } // -x face
            case 2: {
                D = -low[1]; N[1] = -1.0; break; } // -y face
            case 3: {
                D = -low[2]; N[2] = -1.0; break; } // -z face
            case 4: {
                D =  upp[0]; N[0] =  1.0; break; } // +x face
            case 5: {
                D =  upp[1]; N[1] =  1.0; break; } // +y face
            case 6: {
                D =  upp[2]; N[2] =  1.0; }        // +z face    
        }
        
        // Get intersection point between face plane and ray. If no intersection is 
        // possible (i.e. the normal of the face is perp. to the line) then skip face
        sDen = dotProduct(ray.dir,N);
        if (fabs(sDen)>tol) {
        
            // Find intersection point p
            sNum = D - dotProduct(ray.p0,N);
            s    = sNum / sDen;
            p    = vectAdd(ray.p0,ray.dir,s);
            
            // Check if intersection point is within bounds of face. If so, then 
            // return true. If not, then skip face
            if      (faceId==1 || faceId==4) { i=1; j=2; } // -x,+x
            else if (faceId==2 || faceId==5) { i=0; j=2; } // -y,+y
            else if (faceId==3 || faceId==6) { i=0; j=1; } // -z,+z
            if ((p[i]>=low[i] && p[i]<=upp[i]) && (p[j]>=low[j] && p[j]<=upp[j])) {
                return true; }
        }
    }
    return false;
}

// ------------------------------------------------------


cOctree::cOctree(vector<vector<double> > _vertexCoords3D, vector<vector<int> > _polyConnectivity, vector<vector<int> > _quadlist)
{
    vertexCoords3D    = _vertexCoords3D;
    polyConnectivity  = _polyConnectivity;
    quadlist = _quadlist;
    int _offsets[][3] = {{-1,-1,-1},{+1,-1,-1},{-1,+1,-1},{+1,+1,-1},
                         {-1,-1,+1},{+1,-1,+1},{-1,+1,+1},{+1,+1,+1}};
    
    for (unsigned int i=0; i<8; i++) {
        for (unsigned int j=0; j<3; j++) {
            branchOffsets[i][j] = _offsets[i][j];
        }
    }

    setupPolyList(); 
    if (quadlist.size()>0){   
        setupQuadList();
    }
    vector<double> position = getPositionRoot();
    double size = getSizeRoot();
    root = cOctNode(1,"0", position, size);
    insertPolys();
}

void cOctree::setupPolyList()
{
    int indx;
    vector<vector<double> > vertices(3,vector<double>(3,0.0));
    
    polyList.reserve(polyConnectivity.size());
    for (unsigned int i=0; i<polyConnectivity.size(); i++) {
        for (unsigned int j=0; j<3; j++) {
            indx = polyConnectivity[i][j];
            vertices[j] = vertexCoords3D[indx]; }
        polyList.push_back(cTri(i,vertices)); 
    }
}

void cOctree::setupQuadList()
{
    unsigned int i, j;
    vector<vector<double> > vertices(4,vector<double>(3,0.0));
    quadgroup.reserve(quadlist.size());
    for (i=0;i<quadlist.size();i++){
        for(j=0;j<4;j++){
            vertices[j] = vertexCoords3D[quadlist[i][j]];
        }
        quadgroup.push_back(cQuad(i,vertices,0.0,1.0));
    }
}

int cOctree::numPolys() { return (int)(polyList.size()); }

void cOctree::insertPoly(cOctNode &node, cTri &poly)
{
    if (node.isLeafNode()) {
    
        if (poly.isInNode(node)) {
        
            if (node.numPolys() < node.MAX_OCTNODE_OBJECTS) {
                node.addPoly(poly.label);
            } else {
                node.addPoly(poly.label);
                if (node.level < MAX_OCTREE_LEVELS) {
                    splitNodeAndReallocate(node);
                }
            }
        }
        
    } else {
      
        for (unsigned int i=0; i<node.branches.size(); i++) {
            insertPoly(node.branches[i],poly);
        }
        
    }
}

void cOctree::insertPolys()
{
    for (int i=0; i<numPolys(); i++) {
        insertPoly(root,polyList[i]);
    }
}

vector<double> cOctree::getPositionRoot() {

    // Get low and upp
    vector<double> low, upp, position(3);
    low = vertexCoords3D[0];
    upp = vertexCoords3D[0];
    for (unsigned int i=1; i<vertexCoords3D.size(); i++) {
        for (unsigned int j=0; j<3; j++) {
            if (vertexCoords3D[i][j] < low[j]) { low[j] = vertexCoords3D[i][j]; }
            if (vertexCoords3D[i][j] > upp[j]) { upp[j] = vertexCoords3D[i][j]; }
        }
    }
    // Center of node is average of low and upp
    for (unsigned int i=0; i<3; i++) {
        position[i] = 0.5 * (low[i]+upp[i]);
    }
    return position;
}

double cOctree::getSizeRoot() {

    // Get low and upp
    vector<double> low, upp, range;
    low = vertexCoords3D[0];
    upp = vertexCoords3D[0];
    for (unsigned int i=1; i<vertexCoords3D.size(); i++) {
        for (unsigned int j=0; j<3; j++) {
            if (vertexCoords3D[i][j] < low[j]) { low[j] = vertexCoords3D[i][j]; }
            if (vertexCoords3D[i][j] > upp[j]) { upp[j] = vertexCoords3D[i][j]; }
        }
    }
    // Range is the size of the node in each coord direction
    range = vectSubtract(upp,low);
    double size = range[0];
    for (unsigned int i=1; i<3; i++) {
        if (range[i] > size) { size = range[i]; }
    }
    // Scale up size of node by 5%
    size *= 1.05;
    return size;
}

void cOctree::splitNodeAndReallocate(cOctNode &node)
{
    // Split node into 8 branches
    vector<double> position(3);
    for (unsigned int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
        for (unsigned int j=0; j<3; j++) {
            position[j] = node.position[j] + 0.25*node.size*branchOffsets[i][j]; }
        string nid = node.nid + "-" + NumberToString(i);
        node.addNode(node.level+1,nid,position,0.5*node.size); 
    }
    
    // Reallocate date from node to branches
    for (unsigned int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
        for (int j=0; j<node.numPolys(); j++) {
            int indx = node.data[j];
            if (polyList[indx].isInNode(node.branches[i])) {
                if (node.branches[i].numPolys() < node.MAX_OCTNODE_OBJECTS) {
                    node.branches[i].addPoly(indx);
                } else {
                    splitNodeAndReallocate(node.branches[i]);
                }
            }
        }
    }
    node.data.resize(0);
}

vector<cOctNode*> cOctree::getNodesFromLabel(int polyLabel)
{
    // Function for finding all the nodes that contains tri with given label 
    vector<cOctNode*> nodeList;
    findBranchesByLabel(polyLabel,root,nodeList);
    return nodeList;
}

void cOctree::findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList)
{
    // Recursive function used by getNodesFromLabel
    if (node.isLeafNode()) {
        vector<int>::iterator it;
        it = find(node.data.begin(),node.data.end(),polyLabel);
        if (it != node.data.end()) { nodeList.push_back(&node); }
    } else {
        for (unsigned int i=0; i<node.branches.size(); i++) {
            findBranchesByLabel(polyLabel, node.branches[i], nodeList);
        }
    }
}

cOctNode* cOctree::getNodeFromId(string nodeId)
{
    // cout << nodeId;
    return findBranchById(nodeId,root);
}

cOctNode* cOctree::findBranchById(string nodeId, cOctNode &node)
{
    if (nodeId.compare(node.nid)==0) {
        return &node;
    } else {
        for (unsigned int i=0; i<node.branches.size(); i++) {
            cOctNode *branch = findBranchById(nodeId, node.branches[i]);
            if (branch != NULL) { return branch; }
        }
    }
    return NULL;
}

cOctree::~cOctree() 
{
    //cout << "Destroying the cOctree" << endl;
}
//----------------------------------------------------------------------------------------------------
// Ning's code starts here
// get average polygon size
double cOctree::getavgpolysize(cOctNode &node)
{
    double sum = 0.0;
    for (unsigned int i=0; i<polyConnectivity.size(); i++){
        vector<vector<double> > v = polyList[i].vertices;
        double size1 = sqrt((v[0][0] - v[1][0])*(v[0][0] - v[1][0])+(v[0][1] - v[1][1])*(v[0][1] - v[1][1])+(v[0][2] - v[1][2])*(v[0][2] - v[1][2]));
        double size2 = sqrt((v[2][0] - v[1][0])*(v[2][0] - v[1][0])+(v[2][1] - v[1][1])*(v[2][1] - v[1][1])+(v[2][2] - v[1][2])*(v[2][2] - v[1][2]));
        double size3 = sqrt((v[0][0] - v[2][0])*(v[0][0] - v[2][0])+(v[0][1] - v[2][1])*(v[0][1] - v[2][1])+(v[0][2] - v[2][2])*(v[0][2] - v[2][2]));
        double size = (size1+size2+size3)/3;
        sum += size;
    }
    double avg = sum/polyConnectivity.size();
    // double avg=1;
    return avg;
}
// get the maximum ploygon size
double cOctree::getmaxpolysize()
{
    // compute the maximum element size of all the triangles
    double maxsize = 1e-6;
    for (unsigned int i=0; i<root.data.size(); i++){
        vector<vector<double> > v = polyList[*(root.data.begin()+i)].vertices;
        // vector<double> v1 = vectSubtract(v[0], v[1]);
        // vector<double> v2 = vertSubtract(v[1], v[2]);
        // vector<double> v3 = vertSubtract(v[0], v[2]);
        for (unsigned int i=0; i<3; i++){
            double size1 = sqrt((v[0][0] - v[1][0])*(v[0][0] - v[1][0])+(v[0][1] - v[1][1])*(v[0][1] - v[1][1])+(v[0][2] - v[1][2])*(v[0][2] - v[1][2]));
            double size2 = sqrt((v[2][0] - v[1][0])*(v[2][0] - v[1][0])+(v[2][1] - v[1][1])*(v[2][1] - v[1][1])+(v[2][2] - v[1][2])*(v[2][2] - v[1][2]));
            double size3 = sqrt((v[0][0] - v[2][0])*(v[0][0] - v[2][0])+(v[0][1] - v[2][1])*(v[0][1] - v[2][1])+(v[0][2] - v[2][2])*(v[0][2] - v[2][2]));
            double size = max(max(size1, size2), size3);
            if (maxsize < size){
                maxsize = size;
            }   
        }
    }
    return maxsize;
}
// calculate distance between a node and a projection point
double cOctree::getdistance(cOctNode &node, vector<double> projpoint)
{
    node.getLowUppVerts();//get low and upp vector of the node
    double distance;
    double dx = max(max(node.low[0]-projpoint[0], 0.0), projpoint[0]-node.upp[0]);
    double dy = max(max(node.low[1]-projpoint[1], 0.0), projpoint[1]-node.upp[1]);
    double dz = max(max(node.low[2]-projpoint[2], 0.0), projpoint[2]-node.upp[2]);   
    distance = sqrt(dx*dx+dy*dy+dz*dz) - getmaxpolysize();
    if (distance < 0){
        distance = 0;
    }
    return distance;
}

// find the closest point on the mesh to the projection point (recursive function, if unable to understand, ask Ning)
Intersection cOctree::findmindist(cOctNode &node, vector<double> projpoint)
{   
    if (node.isLeafNode()){// if current node is a leafnode great we just need to loop through all the triangles in the node
        if (node.numPolys()>0){// some leafnode has zero trianlges which are the nodes around the corner of the octree
            Distandpos mindistandpos = findtridist(projpoint, polyList[*(node.data.begin())]);
            double mindistri = mindistandpos.distance;
            vector<double> PP0 = mindistandpos.position;
            int k=0;
            for (unsigned int i=0; i<node.data.size(); i++){// find the trianlge that gives the smallest distance and corresponding location
                Distandpos distandpos = findtridist(projpoint, polyList[*(node.data.begin()+i)]);
                if (distandpos.distance < mindistri){
                    mindistri = distandpos.distance;
                    PP0 = distandpos.position;
                    k=*(node.data.begin()+i);
                }
            }
            Intersection mindistandposall(k,PP0,mindistri);
            return mindistandposall;// for recursive function, the result may return to the upper hierachy
        }else{// if the node has zero triangles, just return a very large value of distance
            double mindistri = 1e10;
            vector<double> PP0(3,0.);
            int k=0;
            Intersection mindistandposall(k,PP0,mindistri);
            return mindistandposall;
        }
    }else{// if not, we need to go deeper until the leaf node
        double mindistbox = getdistance(node.branches[0], projpoint);
        vector<Indexanddist> distances;
        unsigned int minbranchnum=0;

        // find the cube that gives the smallest distance to the projpoint
        for (unsigned int i=1; i<node.branches.size(); i++) {
            double distboxi = getdistance(node.branches[i], projpoint);
            distances.push_back(Indexanddist(i,distboxi));
            if (distboxi < mindistbox){
                mindistbox = distboxi; 
                minbranchnum = i;                               
            }
        }
        Intersection mindistandpos = findmindist(node.branches[minbranchnum], projpoint);// recursively use the function to find the smallest distance in the node
        double mindistri = mindistandpos.s;
        vector<double> PP0 = mindistandpos.p;
        int k = mindistandpos.triLabel;
        // check if the found distance is smaller than the distance to the other nodes, if not, we need to check the other nodes for result too
        for (unsigned int i=0; i<node.branches.size(); i++) {
            if (distances[i].distance < mindistri){
                Intersection mindistandposi = findmindist(node.branches[i], projpoint);
                double disttrii = mindistandposi.s;
                if (disttrii < mindistri){
                    mindistri = disttrii;
                    PP0 = mindistandposi.p;
                    k = mindistandposi.triLabel;
                }
            }
        
        }
        Intersection mindistandposall(k,PP0,mindistri);
        return mindistandposall;
    }
}

// for a list of points find their corresponding cloestest points
vector<Intersection> cOctree::findmindists(cOctNode &node, vector<vector<double> > projpoints)
{
    vector<Intersection> distandpositions;
    for (unsigned int i=0; i<projpoints.size(); i++){
        vector<double> projpoint = projpoints[i];
        Intersection distandpos = findmindist(node, projpoint);
        distandpositions.push_back(distandpos);
    }
    return distandpositions;
}

// find the smallest distance between a triangle and a point
Distandpos cOctree::findtridist(vector<double> projpoint, cTri &poly)
{
    // The algorithm is based on
    // "David Eberly, 'Distance Between Point and Triangle in 3D',
    // Geometric Tools, LLC, (1999)"
    // http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    vector<vector<double> > vertices = poly.vertices;
    vector<double> P = projpoint;
    vector<double> B = vertices[0];
    vector<double> E0 = vectSubtract(vertices[1], B);
    vector<double> E1 = vectSubtract(vertices[2], B);
    vector<double> D = vectSubtract(B, P);
    double a = dotProduct(E0, E0);
    double b = dotProduct(E0, E1);
    double c = dotProduct(E1, E1);
    double d = dotProduct(E0, D);
    double e = dotProduct(E1, D);
    double f = dotProduct(D, D);
    double det = abs(a*c -b*b);
    double s = b*e - c*d;
    double t = b*d - a*e;
    double sqrdistance;
    // tree of conditions to determine in which region of the diagram the projection of the point into th triangle-plane lies.
    if ((s+t) <= det){
        if (s < 0.0){
            if (t < 0.0){
                // region 4
                if (d<0){
                    t = 0.0;
                    if ((-d) >= a){
                        s = 1.0;
                        sqrdistance = a + 2.0*d +f;
                    }else{
                        s = -d / a;
                        sqrdistance = d*s + f;
                    }
                }else{
                    s = 0.0;
                    if (e >= 0.0){
                        t = 0.0;
                        sqrdistance = f;
                    }else if ((-e) >= c){
                        t=1.0;
                        sqrdistance = c + 2.0*e +f;
                    }else{
                        t = -e/c;
                        sqrdistance = e*t + f;

                        // region 4
                    }
                }
            }else{
                // region 3
                s = 0.0;
                if (e >= 0.0){
                    t = 0.0;
                    sqrdistance = f;
                }else if((-e) >= c){
                    t = 1.0;
                    sqrdistance = c + 2.0*e + f;
                }else{
                    t = -e/c;
                    sqrdistance = e*t + f;
                    // region 3
                }
            }
        }else if (t < 0.0){
            // region 5
            t = 0.0;
            if (d >= 0.0){
                s = 0.0;
                sqrdistance = f;
            }else if ((-d) >= a){
                s = 1.0;
                sqrdistance = a + 2.0*d + f;
            }else{
                s = -d/a;
                sqrdistance = d*s +f;
            }
        }else{
            // region 0
            double invDet = 1.0/det;
            s = s*invDet;
            t = t*invDet;
            sqrdistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f;
        }
    }else{
        if (s < 0.0){
            // region 2
            double tmp0 = b + d;
            double tmp1 = c + e;
            if (tmp1 > tmp0){
                double numer = tmp1 - tmp0;
                double denom = a - 2.0*b +c;
                if (numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqrdistance = a + 2.0*d + f;
                }else{
                    s =numer/denom;
                    t = 1 -s;
                    sqrdistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2*e) + f;
                }
            }else{
                s = 0.0;
                if (tmp1 <= 0.0){
                    t =1.0;
                    sqrdistance = c + 2.0*e +f;
                }else if (e >= 0.0){
                    t = 0.0;
                    sqrdistance = f;
                }else{
                    t = -e/c;
                    sqrdistance = e*t + f;
                    // region 2
                }
            }
        }else if (t < 0.0){
            // region 6
            double tmp0 = b + e;
            double tmp1 = a + d;
            if (tmp1 > tmp0){
                double numer = tmp1 - tmp0;
                double denom = a - 2.0*b + c;
                if (numer >= denom){
                    t = 1.0;
                    s = 0.0;
                    sqrdistance = c + 2.0*e +f;
                }else{
                    t = numer/denom;
                    s = 1.0 - t;
                    sqrdistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f;
                }
            }else{
                t = 0.0;
                if (tmp1 <= 0.0){
                    s = 1.0;
                    sqrdistance = a+ 2.0*d + f;
                }else if (d >= 0.0){
                    s = 0.0;
                    sqrdistance =f;
                }else{
                    s  = -d/a;
                    sqrdistance = d*s + f;
                }
            }
        }else{
            // region 1
            double numer = c + e - b - d;
            if (numer <= 0.0){
                s = 0.0;
                t = 1.0;
                sqrdistance = c + 2.0*e + f;
            }else{
                double denom = a - 2.0*b +c;
                if (numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqrdistance = a + 2.0*d + f; 
                }else{
                    s = numer/denom;
                    t = 1-s;
                    sqrdistance = s*(a*s + b*t + 2.0*d) + t *(b*s + c*t +2.0*e) + f;
                }
            }
        }
        
    }
    sqrdistance = s*(a*s + b*t + 2.0*d) + t *(b*s + c*t +2.0*e) + f;
    if (sqrdistance < 0.0){
        sqrdistance = 0.0;
    }
    double dist = sqrt(sqrdistance);
    // transform(E0.begin(), E0.end(), E0.begin(), bind(multiplies<1>(s), placeholders::_1, 3));
    // transform(E1.begin(), E1.end(), E1.begin(), bind(multiplies<1>(t), placeholders::_1, 3));
    vector<double> E00 = vectMuldou(E0, s);
    vector<double> E11 = vectMuldou(E1, t);

    vector<double> PP = vectAdd(B, E00);
    vector<double> PP0 = vectAdd(PP, E11);
    Distandpos distanceandposition(dist, PP0);

    return distanceandposition;
} 

// retriangulation given the projected points
CoordandConnect cOctree::retriangulation(vector<vector<double> > projectedpoint, vector<int> trilist, vector<vector<double> > vertexCoords3D, vector<vector<int> > polyConnectivity)
{
    vector<int> ilist, a1, a2, a3, constrained;
    vector<double> Pi, P1, P2, P3, PiP1, PiP2, PiP3, P1P2, P1P3, cross1, cross2, cross3, cross;
    double area1, area2, area3, area, dist;
    unsigned int i, j, k, l, m, n, o, p;
    vector<vector<double> > projectedlist;
    vector<vector<int> > newtrilist;
    for (i=0; i<projectedpoint.size(); i++){
        if (i==0){
            projectedlist.push_back(projectedpoint[i]);
            ilist.push_back(trilist[i]);
            for (j=i+1; j<projectedpoint.size();j++){
                if (trilist[j]==trilist[i]){
                    projectedlist.push_back(projectedpoint[j]);
                    // ilist.push_back(j);
                }
            }
            newtrilist.push_back(polyConnectivity[trilist[i]]);
            // cout<< k << "   " << l;
            // o=newtrilist.size();
            for (k=0; k<projectedlist.size();k++){
                for (l=0; l<newtrilist.size();l++){
                    Pi=projectedlist[k];
                    P1=vertexCoords3D[newtrilist[l][0]];
                    P2=vertexCoords3D[newtrilist[l][1]];
                    P3=vertexCoords3D[newtrilist[l][2]];    
                    PiP1=vectSubtract(Pi,P1);
                    PiP2=vectSubtract(Pi,P2);
                    PiP3=vectSubtract(Pi,P3);
                    P1P2=vectSubtract(P1,P2);
                    P1P3=vectSubtract(P1,P3);
                    cross1=crossProduct(PiP1,PiP2);
                    cross2=crossProduct(PiP2,PiP3);
                    cross3=crossProduct(PiP3,PiP1);
                    cross=crossProduct(P1P2,P1P3);
                    area1=0.5*sqrt(dotProduct(cross1,cross1));
                    area2=0.5*sqrt(dotProduct(cross2,cross2));
                    area3=0.5*sqrt(dotProduct(cross3,cross3));
                    area=0.5*sqrt(dotProduct(cross,cross));
                    if (abs(area-area1-area2-area3)<1e-10){
                        for (p=0;p<3;p++){
                            // P1, P2, P3
                            dist=distBetweenPoints(vertexCoords3D[newtrilist[l][p]],Pi);
                            if (dist<1e-10){
                                constrained.push_back(newtrilist[l][p]);
                                goto here1;
                            }
                        }
                        constrained.push_back(vertexCoords3D.size()-1);
                        vertexCoords3D.push_back(Pi);
                        a1.clear();
                        a1.push_back(newtrilist[l][0]);
                        a1.push_back(newtrilist[l][1]);
                        a1.push_back((int) vertexCoords3D.size()-1);
                        a2.clear();
                        a2.push_back(newtrilist[l][1]);
                        a2.push_back(newtrilist[l][2]);
                        a2.push_back((int) vertexCoords3D.size()-1);
                        a3.clear();
                        a3.push_back(newtrilist[l][2]);
                        a3.push_back(newtrilist[l][0]);
                        a3.push_back((int) vertexCoords3D.size()-1);
                        // cout<<a1[0]<<a1[1]<<a1[2];
                        // cout<<newtrilist[0][0];
                        newtrilist.push_back(a1);
                        newtrilist.push_back(a2);
                        newtrilist.push_back(a3);
                        newtrilist.erase(newtrilist.begin()+l);
                        here1:;
                        break;
                    }
                }
            }
            // polyConnectivity.erase(polyConnectivity.begin()+trilist[i]);
            for (m=0; m<newtrilist.size();m++){
                polyConnectivity.push_back(newtrilist[m]);
            }
        }else{
            for (n=0; n<ilist.size();n++){
                if (trilist[i]==ilist[n]){
                    goto endloop;
                }
            }
            projectedlist.clear();
            projectedlist.push_back(projectedpoint[i]);
            ilist.push_back(trilist[i]);
            for (j=i+1; j<projectedpoint.size();j++){
                if (trilist[j]==trilist[i]){
                    projectedlist.push_back(projectedpoint[j]);
                }
            }
            newtrilist.clear();
            newtrilist.push_back(polyConnectivity[trilist[i]]);
            // o=newtrilist.size();
            for (k=0; k<projectedlist.size();k++){
                for (l=0; l<newtrilist.size();l++){
                    // cout<<i;
                    Pi=projectedlist[k];
                    P1=vertexCoords3D[newtrilist[l][0]];
                    P2=vertexCoords3D[newtrilist[l][1]];
                    P3=vertexCoords3D[newtrilist[l][2]];    
                    PiP1=vectSubtract(Pi,P1);
                    PiP2=vectSubtract(Pi,P2);
                    PiP3=vectSubtract(Pi,P3);
                    P1P2=vectSubtract(P1,P2);
                    P1P3=vectSubtract(P1,P3);
                    cross1=crossProduct(PiP1,PiP2);
                    cross2=crossProduct(PiP2,PiP3);
                    cross3=crossProduct(PiP3,PiP1);
                    cross=crossProduct(P1P2,P1P3);
                    area1=0.5*sqrt(dotProduct(cross1,cross1));
                    area2=0.5*sqrt(dotProduct(cross2,cross2));
                    area3=0.5*sqrt(dotProduct(cross3,cross3));
                    area=0.5*sqrt(dotProduct(cross,cross));
                    if (abs(area-area1-area2-area3)<1e-9){
                        for (p=0;p<3;p++){
                            // P1, P2, P3
                            dist=distBetweenPoints(vertexCoords3D[newtrilist[l][p]],Pi);
                            if (dist<1e-6){
                                constrained.push_back(newtrilist[l][p]);
                                goto here2;
                            }
                        }
                        constrained.push_back(vertexCoords3D.size()-1);
                        vertexCoords3D.push_back(Pi);
                        a1.clear();
                        a1.push_back(newtrilist[l][0]);
                        a1.push_back(newtrilist[l][1]);
                        a1.push_back((int) vertexCoords3D.size()-1);
                        a2.clear();
                        a2.push_back(newtrilist[l][1]);
                        a2.push_back(newtrilist[l][2]);
                        a2.push_back((int) vertexCoords3D.size()-1);
                        a3.clear();
                        a3.push_back(newtrilist[l][2]);
                        a3.push_back(newtrilist[l][0]);
                        a3.push_back((int) vertexCoords3D.size()-1);
                        newtrilist.push_back(a1);
                        newtrilist.push_back(a2);
                        newtrilist.push_back(a3);
                        newtrilist.erase(newtrilist.begin()+l);
                        here2:;
                        break;   
                    }
                }
            }
            // polyConnectivity.erase(polyConnectivity.begin()+trilist[i]);
            for (m=0; m<newtrilist.size();m++){
                polyConnectivity.push_back(newtrilist[m]);
            }
        endloop:;
        }
    }
    sort(ilist.begin(), ilist.end(), greater<int>());
    for (o=0;o<ilist.size();o++){
        polyConnectivity.erase(polyConnectivity.begin()+ilist[o]);
    }
    CoordandConnect coordandconnect(vertexCoords3D,polyConnectivity,constrained);
    return coordandconnect;
}

// count all the edges in the triangular mesh
vector<vector<int> > cOctree::countedges()
{   unsigned int i;
    vector<int> v1, v2, v3;
    vector<vector<int> > edges;
    for (i=0;i<polyConnectivity.size();i++){
        // cout<<polyConnectivity[i][0]<<'\n';
        // cout<<polyConnectivity[i][1]<<'\n';
        // cout<<polyConnectivity[i][2]<<'\n';
        v1.clear();
        v2.clear();
        v3.clear();
        v1.push_back(polyConnectivity[i][0]);
        v1.push_back(polyConnectivity[i][1]);
        v2.push_back(polyConnectivity[i][1]);
        v2.push_back(polyConnectivity[i][2]);
        v3.push_back(polyConnectivity[i][2]);
        v3.push_back(polyConnectivity[i][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        edges.push_back(v1);
        edges.push_back(v2);
        edges.push_back(v3);
    }
    sort(edges.begin(),edges.end());
    edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return edges;
}

// count all edges given the trilist
vector<vector<int> > cOctree::countedges(vector<vector<int> > trilist)
{   
    unsigned int i;
    vector<int> v1, v2, v3;
    vector<vector<int> > edges;
    for (i=0;i<trilist.size();i++){
        // cout<<trilist[i][0]<<'\n';
        // cout<<trilist[i][1]<<'\n';
        // cout<<trilist[i][2]<<'\n';
        v1.clear();
        v2.clear();
        v3.clear();
        v1.push_back(trilist[i][0]);
        v1.push_back(trilist[i][1]);
        v2.push_back(trilist[i][1]);
        v2.push_back(trilist[i][2]);
        v3.push_back(trilist[i][2]);
        v3.push_back(trilist[i][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        edges.push_back(v1);
        edges.push_back(v2);
        edges.push_back(v3);
    }
    sort(edges.begin(),edges.end());
    edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return edges;
}

// count all edges considering considering the normal vectors of the triangles that they are located in
vector<vector<int> > cOctree::countedges(vector<vector<int> > trilist, vector<vector<double> > vertlist,vector<double> revprojdir)
{
    unsigned int i, j;
    vector<int> v1, v2, v3;
    vector<vector<int> > edges;
    for (i=0;i<trilist.size();i++){
        vector<vector<double> > v;
        for (j=0;j<3;j++){ 
            v.push_back(vertlist[trilist[i][j]]);
        }
        vector<double> vv1(3),vv2(3),vv3;
        for (unsigned int k=0; k<3; k++) {
            vv1[k] = v[0][k] - v[1][k];
            vv2[k] = v[0][k] - v[2][k]; }
        vv3 = crossProduct(vv1,vv2);
        double v3mag = sqrt(dotProduct(vv3,vv3));
        vector<double> N = vv3;
        for (vector<double>::iterator it=N.begin(); it!=N.end(); ++it)
            *it /= v3mag;
        if (dotProduct(N, revprojdir)>0){
            v1.clear();
            v2.clear();
            v3.clear();
            v1.push_back(trilist[i][0]);
            v1.push_back(trilist[i][1]);
            v2.push_back(trilist[i][1]);
            v2.push_back(trilist[i][2]);
            v3.push_back(trilist[i][2]);
            v3.push_back(trilist[i][0]);
            sort(v1.begin(), v1.end());
            sort(v2.begin(), v2.end());
            sort(v3.begin(), v3.end());
            edges.push_back(v1);
            edges.push_back(v2);
            edges.push_back(v3);
        }
    }
    sort(edges.begin(),edges.end());
    edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return edges;
}

// compute the normal vector of each vertex
vector<vector<double> > cOctree::computevertexnorm(vector<vector<int> > polyConnectivity, vector<vector<double> > vertexCoords3D)
{
    unsigned int i;
    vector<vector<double> > norms(vertexCoords3D.size()); 
    for (i=0;i<vertexCoords3D.size();i++){
        norms[i].resize(3,0.0);
    }
    for (i=0;i<polyConnectivity.size();i++){
        vector<double> v1 = vertexCoords3D[polyConnectivity[i][0]];
        vector<double> v2 = vertexCoords3D[polyConnectivity[i][1]];
        vector<double> v3 = vertexCoords3D[polyConnectivity[i][2]];
        vector<double> edge1 = vectSubtract(v2,v1);
        vector<double> edge2 = vectSubtract(v3,v2);
        vector<double> norm = crossProduct(edge1,edge2);
        norm = normalize(norm);
        norms[polyConnectivity[i][0]]= vectAdd(norms[polyConnectivity[i][0]],norm);
        norms[polyConnectivity[i][1]]= vectAdd(norms[polyConnectivity[i][1]],norm);
        norms[polyConnectivity[i][2]]= vectAdd(norms[polyConnectivity[i][2]],norm);
    }
    for (i=0;i<vertexCoords3D.size();i++){
        norms[i]=normalize(norms[i]);
    }
    return norms;
}

// used for openmdao nlp, no longer useful, but may be applied in the future
ProjFlatten cOctree::projlist(cOctNode &node, vector<vector<double> > projpoints)
{
    vector<double> proj;
    vector<double> gradient;
    vector<Projection> projections = project(node, projpoints);
    for (unsigned int i=0;i<projections.size();i++){
        for (unsigned int j=0;j<3;j++){
            proj.push_back(projections[i].phycoord[j]);
            for (unsigned int k=0;k<3;k++){
                gradient.push_back(projections[i].gradient[j][k]);
            }
        }
    }
    ProjFlatten projflatten(proj,gradient);
    return projflatten;
}

// used for openmdao nlp, no longer useful, but may be applied in the future
vector<Projection> cOctree::project(cOctNode &node, vector<vector<double> > projpoints)
{
    unsigned int i;
    vector<Projection> projections;
    for (i=0;i<projpoints.size();i++){
        // cout<<i<<'\n';
        vector<double> projpoint = projpoints[i];
        Intersection intersect = findmindist(node,projpoint);
        if (quadlist.size()>0){
            Projection projection = quadproject(projpoint,i);
            if (intersect.s<projection.distance){
                // cout<<'1';
                int list = 0;
                int label = intersect.triLabel;
                double distance = intersect.s;
                vector<double> phycoord = intersect.p;
                vector<double> P0=vertexCoords3D[polyConnectivity[label][0]];
                vector<double> P1=vertexCoords3D[polyConnectivity[label][1]];
                vector<double> P2=vertexCoords3D[polyConnectivity[label][2]];  
                vector<double> PiP0=vectSubtract(phycoord,P0);  
                vector<double> PiP1=vectSubtract(phycoord,P1);
                vector<double> PiP2=vectSubtract(phycoord,P2);
                vector<double> P0P1=vectSubtract(P0,P1);
                vector<double> P0P2=vectSubtract(P0,P2);
                vector<double> cross0=crossProduct(PiP1,PiP2);
                vector<double> cross1=crossProduct(PiP0,PiP2);
                vector<double> cross=crossProduct(P0P1,P0P2);
                double area0=0.5*sqrt(dotProduct(cross0,cross0));
                double area1=0.5*sqrt(dotProduct(cross1,cross1));
                double area=0.5*sqrt(dotProduct(cross,cross));
                vector<double> paracoord;
                double u = area0/area; double v = area1/area; double w = 1-u-v;
                paracoord.push_back(u);paracoord.push_back(v);
                vector<double> f = vectSubtract(phycoord,projpoint);
                vector<vector<double> > dRdx(2,vector<double>(3,0.0));
                vector<vector<double> > dRdy(2,vector<double>(2,0.0));
                vector<vector<double> > dFdy(3,vector<double>(2,0.0));
                vector<vector<double> > dydx(2,vector<double>(3,0.0));
                vector<vector<double> > dfdx(3,vector<double>(3,0.0));
                vector<vector<double> > invdRdy(2,vector<double>(2,0.0));
                vector<double> pv = vectSubtract(P1,P2);
                vector<double> pu = vectSubtract(P0,P2);
                vector<double> pw = vectSubtract(P1,P0);
                if ((u+v<1e-8)||(v+w<1e-8)||(w+u<1e-8)){// vertice

                }else if (u<1e-8){// edge 0
                    dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
                    dRdy[0][0] = 1; dRdy[1][0] = dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
                    dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                    invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                    dfdx = vectMultiply(dFdy,dydx);
                }else if (v<1e-8){
                    dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2];
                    dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(pu,pv); dRdy[1][1] = 1;
                    dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                    invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                    dfdx = vectMultiply(dFdy,dydx);
                }else if (w<1e-8){
                    dRdx[0][0] = -pw[0]; dRdx[0][1] = -pw[1]; dRdx[0][2] = -pw[2]; vector<double> negpu = vectneg(pu);
                    dRdy[0][0] = dotProduct(pw,pw); dRdy[0][1] = dotProduct(pw,negpu); dRdy[1][1] = 1;
                    dFdy[0][0] = pw[0]; dFdy[1][0] = pw[1]; dFdy[2][0] = pw[2]; dFdy[0][1] = negpu[0]; dFdy[1][1] = negpu[1]; dFdy[2][1] = negpu[2];
                    // dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                    invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                    dfdx = vectMultiply(dFdy,dydx);
                }else{
                    dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2]; dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
                    dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(pu,pv); dRdy[1][0] = dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
                    dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                    invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                    dfdx = vectMultiply(dFdy,dydx);
                }
                Projection projection(list,label,distance,phycoord,paracoord,dfdx);
                projections.push_back(projection);
            }else{
                // cout<<'2';
                projections.push_back(projection);
            }
        }else{
            int list = 0;
            int label = intersect.triLabel;
            double distance = intersect.s;
            vector<double> phycoord = intersect.p;
            vector<double> P0=vertexCoords3D[polyConnectivity[label][0]];
            vector<double> P1=vertexCoords3D[polyConnectivity[label][1]];
            vector<double> P2=vertexCoords3D[polyConnectivity[label][2]];  
            vector<double> PiP0=vectSubtract(phycoord,P0);  
            vector<double> PiP1=vectSubtract(phycoord,P1);
            vector<double> PiP2=vectSubtract(phycoord,P2);
            vector<double> P0P1=vectSubtract(P0,P1);
            vector<double> P0P2=vectSubtract(P0,P2);
            vector<double> cross0=crossProduct(PiP1,PiP2);
            vector<double> cross1=crossProduct(PiP0,PiP2);
            vector<double> cross=crossProduct(P0P1,P0P2);
            double area0=0.5*sqrt(dotProduct(cross0,cross0));
            double area1=0.5*sqrt(dotProduct(cross1,cross1));
            double area=0.5*sqrt(dotProduct(cross,cross));
            vector<double> paracoord;
            double u = area0/area; double v = area1/area; double w = 1-u-v;
            paracoord.push_back(u);paracoord.push_back(v);
            vector<double> f = vectSubtract(phycoord,projpoint);
            vector<vector<double> > dRdx(2,vector<double>(3,0.0));
            vector<vector<double> > dRdy(2,vector<double>(2,0.0));
            vector<vector<double> > dFdy(3,vector<double>(2,0.0));
            vector<vector<double> > dydx(2,vector<double>(3,0.0));
            vector<vector<double> > dfdx(3,vector<double>(3,0.0));
            vector<vector<double> > invdRdy(2,vector<double>(2,0.0));
            vector<double> pv = vectSubtract(P1,P2);
            vector<double> pu = vectSubtract(P0,P2);
            vector<double> pw = vectSubtract(P1,P0);
            if ((u+v<1e-8)||(v+w<1e-8)||(w+u<1e-8)){// vertice

            }else if (u<1e-8){// edge 0
                dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
                dRdy[0][0] = 1; dRdy[1][0] = dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
                dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                dfdx = vectMultiply(dFdy,dydx);
            }else if (v<1e-8){
                dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2];
                dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(pu,pv); dRdy[1][1] = 1;
                dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                dfdx = vectMultiply(dFdy,dydx);
            }else if (w<1e-8){
                dRdx[0][0] = -pw[0]; dRdx[0][1] = -pw[1]; dRdx[0][2] = -pw[2]; vector<double> negpu = vectneg(pu);
                dRdy[0][0] = dotProduct(pw,pw); dRdy[0][1] = dotProduct(pw,negpu); dRdy[1][1] = 1;
                dFdy[0][0] = pw[0]; dFdy[1][0] = pw[1]; dFdy[2][0] = pw[2]; dFdy[0][1] = negpu[0]; dFdy[1][1] = negpu[1]; dFdy[2][1] = negpu[2];
                // dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                dfdx = vectMultiply(dFdy,dydx);
            }else{
                dRdx[0][0] = -pu[0]; dRdx[0][1] = -pu[1]; dRdx[0][2] = -pu[2]; dRdx[1][0] = -pv[0]; dRdx[1][1] = -pv[1]; dRdx[1][2] = -pv[2];
                dRdy[0][0] = dotProduct(pu,pu); dRdy[0][1] = dotProduct(pu,pv); dRdy[1][0] = dotProduct(pu,pv); dRdy[1][1] = dotProduct(pv,pv);
                dFdy[0][0] = pu[0]; dFdy[1][0] = pu[1]; dFdy[2][0] = pu[2]; dFdy[0][1] = pv[0]; dFdy[1][1] = pv[1]; dFdy[2][1] = pv[2]; 
                invdRdy = vectInverse(dRdy); dydx = vectMultiply(invdRdy,dRdx); dydx = vectneg(dydx);
                dfdx = vectMultiply(dFdy,dydx);
            }
            Projection projection(list,label,distance,phycoord,paracoord,dfdx);
            projections.push_back(projection);
        }
    }
    // cout<<projections.size();
    return projections;
}

// used for openmdao nlp, no longer useful, but may be applied in the future
Projection cOctree::quadproject(vector<double> p0, int index)
{
    // cout<<index<<'\n';
    unsigned int i, j, k, l, m, n, o, count;
    int list = 1;
    double min_d_local = 1e8;
    int min_i_local = 0;
    double min_u_local = 0;
    double min_v_local = 0;
    vector<double> min_coord_local(3,0.0);
    vector<double> x(2,0.0);
    for (i=0;i<quadgroup.size();i++){
        double min_d = 1e8;
        double min_u = 0, min_v = 0, u = 0, v = 0;
        vector<double> min_coord(3,0.0);
        double nrefine = 51;
        for (j=0;j<nrefine;j++){
            for (k=0;k<nrefine;k++){
                u = j/(nrefine-1);
                v = k/(nrefine-1);
                vector<double> p = quadgroup[i].evaluate(u,v);
                vector<double> diff = vectSubtract(p,p0);
                double d = sqrt(dotProduct(diff,diff));
                // cout<<"\nk: "<<k<<'\n';
                // cout<<"u: "<<u<<'\n';
                // cout<<"v: "<<v<<'\n';
                // cout<<"min_d: "<<min_d<<'\n';
                // cout<<"d: "<<d<<'\n';
                if (d<min_d){
                    min_d = d;
                    min_u = u;
                    min_v = v;
                    min_coord = p;
                }
            }
        }
        if (min_d<min_d_local){
            min_d_local = min_d;
            min_i_local = i;
            min_u_local = min_u;
            min_v_local = min_v;
            min_coord_local = min_coord;
        }
        x[0] = min_u;
        int fail = 1;
        double det = 0;
        for (count=0;count<200;count++){
            vector<double> g = quadgroup[i].gradient(x[0],x[1],p0);
            vector<vector<double> > h = quadgroup[i].hessian(x[0],x[1],p0);
            vector<vector<double> > w(2,vector<double>(2,0.0));
            vector<double> dx(2,0.0);
            for (m=0;m<2;m++){
                if (((x[m]==quadgroup[i].lower)&&(g[m]>0))||((x[m]==quadgroup[i].upper)&&(g[m]<0))){
                    g[m] = 0;
                    h[0][1] = 0;
                    h[1][0] = 0;
                    h[m][m] = 1;
                }
            }
            double norm_g = vectNorm(g);
            det = h[0][0]*h[1][1] - h[0][1]*h[1][0];
            if ((det<1e-10)||(isnan(det))||(det==0)){
                dx = g;
            }else{
                w[0][0] = h[1][1]/det;
                w[0][1] = -h[0][1]/det;
                w[1][0] = -h[1][0]/det;
                w[1][1] = h[0][0]/det;
                dx[0] = dotProduct(w[0],g); dx[0] = -dx[0];
                dx[1] = dotProduct(w[1],g); dx[1] = -dx[1];
            }
        
            for (n=0;n<2;n++){
                if (x[n]+dx[n]<quadgroup[i].lower){
                    dx[n] = quadgroup[i].lower-x[n];
                }else if(x[n]+dx[n]>quadgroup[i].upper){
                    dx[n] = quadgroup[i].upper-x[n];
                }
            }
            double norm_dx = vectNorm(dx);
            
            x = vectAdd(x,dx);
            if ((norm_g<1e-10)||(norm_dx<1e-10)){
                fail = 0;
                break;
            }
        }
        if (fail==0){
            if (quadgroup[i].distance(x[0],x[1],p0)<min_d_local){
                min_d_local = quadgroup[i].distance(x[0],x[1],p0);
                min_i_local = i;
                min_u_local = x[0];
                min_v_local = x[1];
                min_coord_local = quadgroup[i].evaluate(x[0],x[1]);
            }
        }else{
            cout<<x[0]<<' '<<x[1]<<'\n';
            cout<<det<<'\n';
            cout<<index<<'\n';
        }
    }
    vector<vector<double> > dfdx = quadgroup[min_i_local].dfdx(min_u_local,min_v_local,p0);
    vector<double> min_para_local;
    min_para_local.push_back(min_u_local);min_para_local.push_back(min_v_local);
    Projection projection(list,min_i_local,min_d_local,min_coord_local,min_para_local,dfdx);
    return projection;
}





// retriangulation
vector<Intersection>  cOctree::ProjectionOfPoints(cOctNode &node, vector<vector<double> > projpoints, vector<vector<double> > projdirts)
{
    unsigned int i, j, k;
    vector<vector<double> > vertlist = vertexCoords3D;
    vector<int> projlist;
    vector<vector<int> > trilist = polyConnectivity;
    vector<Intersection> proj;
    if (projdirts.size()==0){
        proj = findmindists(node,projpoints);
        cout<<"findmindists(node,projpoints);"<<'\n';
    }else{
        vector<cLine> raylist;
        for (i=0;i<projdirts.size();i++){
            cLine ray(projpoints[i],projdirts[i],0);
            raylist.push_back(ray);
        }
        proj = findRayIntersects(raylist);
        //cout<<"findRayIntersects(raylist)"<<'\n';
    }
    // cout<<"update trilist"<<'\n';
    // cout<<"proj"<<'\n';
    // cout<<proj[0].triLabel<<' ';
    // cout<<proj[1].triLabel<<' ';
    // cout<<proj[2].triLabel<<'\n';
    //cout<<"projpoints.size()"<<'\n';
    //cout<<projpoints.size()<<'\n';
    
    return proj;
}


// retriangulation
Reconnection cOctree::reconnect(cOctNode &node, vector<vector<double> > projpoints, vector<vector<double> > projdirts)
{
    unsigned int i, j, k;
    vector<vector<double> > vertlist = vertexCoords3D;
    vector<int> projlist;
    vector<vector<int> > trilist = polyConnectivity;
    vector<Intersection> proj;
    if (projdirts.size()==0){
        proj = findmindists(node,projpoints);
        cout<<"findmindists(node,projpoints);"<<'\n';
    }else{
        vector<cLine> raylist;
        for (i=0;i<projdirts.size();i++){
            cLine ray(projpoints[i],projdirts[i],0);
            raylist.push_back(ray);
        }
        proj = findRayIntersects(raylist);
        //cout<<"findRayIntersects(raylist)"<<'\n';
    }
    cout<<"update trilist"<<'\n';
    cout<<"proj"<<'\n';
    cout<<proj[0].triLabel<<'\n';
    cout<<proj[1].triLabel<<'\n';
    //cout<<"projpoints.size()"<<'\n';
    //cout<<projpoints.size()<<'\n';
    for (i=0;i<projpoints.size()-1;i++){// loop through all the projected points
        vector<int> projlisti;
        vector<vector<double> > inti;
        vector<int> V;
        vector<double> mid = midPoint(projpoints[i],projpoints[i+1]);// find the midpoint of two projection points
        vector<double> midp = midPoint(proj[i].p,proj[i+1].p);// find the midpoint of two projected points
        vector<double> projvector = vectSubtract(proj[i+1].p,proj[i].p);
        vector<double> revprojdir = vectSubtract(mid,midp);
        vector<vector<int> > edges = countedges(trilist, vertlist, revprojdir);
        vector<double> e0, e1;
        vector<double> p0 = proj[i].p;
        vector<double> p1 = proj[i+1].p;
        int count = 0;
        //"<<'\n';
        //cout<<vertlist.size()<<'\n';
        for (j=0;j<vertlist.size();j++){
            if (distBetweenPoints(p0,vertlist[j])<1e-5){
                projlisti.push_back(j);
                break;
            }
            if (j==vertlist.size()-1){
                count = 1;
            }
        }
        if (count==1){
            vertlist.push_back(p0);
            projlisti.push_back(vertlist.size()-1);
        }
        vector<double> normal = crossProduct(revprojdir,projvector); normal = normalize(normal);
        for (j=0;j<edges.size();j++){// loop through all the edges to find the intersection
            e0 = vertlist[edges[j][0]];// the starting vertex of the edge
            e1 = vertlist[edges[j][1]];// the ending vertex of the edge
            vector<double> line = vectSubtract(e0,e1);// vector from e1 to e0
            vector<double> temp1 = vectSubtract(p0,e0);
            double temp2 = dotProduct(temp1,normal);
            double temp3 = dotProduct(line,normal);
            double d = temp2/temp3;
            vector<double> intersect = vectAdd(e0,line,d);// find the intersection point between the plane and the line
            vector<double> ie0 = vectSubtract(e0,intersect), ie1 = vectSubtract(e1,intersect);
            if (dotProduct(ie0,ie1)<0){// see if edge intersects with the plane 
                vector<double> p01 = vectSubtract(p1,p0);
                vector<double> p0int = vectSubtract(intersect,p0);
                double pos = dotProduct(p01,p0int)/dotProduct(p01,p01);
                if (pos<=1 && pos>=0){
                    inti.push_back(intersect);// save intersection point
                    int count = 0;
                    for (k=0;k<vertlist.size();k++){
                        if (distBetweenPoints(intersect,vertlist[k])<1e-5){
                            V.push_back(k);
                            break;
                        }
                        if (k==vertlist.size()-1){
                            count = 1;
                        }
                    }
                    if (count==1){
                        vertlist.push_back(intersect);
                        V.push_back(vertlist.size()-1);
                    }
                }
            }
        }
        vector<double> dists(inti.size(),0.0);
        for (j=0;j<inti.size();j++){
            dists[j] = distBetweenPoints(inti[j],p0);// calculate the distances between the intersection points and projected point
        }
        vector<double> reference; reference.assign(dists.begin(), dists.end());// save the distances as reference
        
        sort(dists.begin(),dists.end());// sort distance in the increasing order
        
        for (j=0;j<inti.size();j++){// sort index by distances 
            for (k=0;k<inti.size();k++){
                if (reference[k]==dists[j]){
                    projlisti.push_back(V[k]);
                }
            }
        }
        if (i==projpoints.size()-2){// at the second last point, add the last point
            int count = 0;
            for (k=0;k<vertlist.size();k++){
                if (distBetweenPoints(proj[i+1].p,vertlist[k])<1e-5){
                    V.push_back(k);
                    break;
                }
                if (k==vertlist.size()-1){
                    count = 1;
                }
            }
            if (count==1){
                vertlist.push_back(proj[i+1].p);
                projlisti.push_back(vertlist.size()-1);
            }
        }
        projlist.insert(projlist.end(),projlisti.begin(),projlisti.end());
        // update trilist
        // cout<<"projlist"<<'\n';
        //cout<<projlist<<'\n';
        // cout<<projlist.size()<<'\n';
        // cout<<projlisti.size()<<'\n';
        // cout<<projlist[0]<<'\n';
        for (i = 0;i< projlisti.size(); i++){
            cout<<projlisti[i];
            cout << " " ;
        }
        // cout<<"projlisti"<<'\n';
        cout<<"update trilist test2"<<'\n';
        for (j=0;j<projlisti.size();j++){
            vector<vector<int> > tripush;
            vector<double> Pi = vertlist[projlisti[j]];
            for (k=0;k<trilist.size();k++){// check the projected point has intersection with which tri
                vector<double> P1=vertlist[trilist[k][0]];
                vector<double> P2=vertlist[trilist[k][1]];
                vector<double> P3=vertlist[trilist[k][2]];    
                vector<double> PiP1=vectSubtract(Pi,P1);
                vector<double> PiP2=vectSubtract(Pi,P2);
                vector<double> PiP3=vectSubtract(Pi,P3);
                vector<double> P1P2=vectSubtract(P1,P2);
                vector<double> P1P3=vectSubtract(P1,P3);
                vector<double> cross1=crossProduct(PiP1,PiP2);
                vector<double> cross2=crossProduct(PiP2,PiP3);
                vector<double> cross3=crossProduct(PiP3,PiP1);
                vector<double> cross=crossProduct(P1P2,P1P3);
                double area1=0.5*sqrt(dotProduct(cross1,cross1));
                double area2=0.5*sqrt(dotProduct(cross2,cross2));
                double area3=0.5*sqrt(dotProduct(cross3,cross3));
                double area=0.5*sqrt(dotProduct(cross,cross));
                // cout<<area1+area2+area3-area<<'\n';
                if ((area1+area2+area3-area)<1e-5){// check area
                    // cout<<area1<<'\n';
                    // cout<<area2<<'\n';
                    // cout<<area3<<'\n';
                    // cout<<area<<'\n';
                    // projection locates in the tri
                    if((abs(area1-area)<1e-5) || (abs(area2-area)<1e-5) || (abs(area3-area)<1e-5)){
                        // cout<<"on vert"<<'\n';
                    }else if ((area1<1e-5)||(area2<1e-5)||(area3<1e-5)){
                        // cout<<"on edge"<<'\n';
                        // projection locates on the edges
                        if (area1<1e-5){
                            int tri1list[] = {trilist[k][0],projlisti[j],trilist[k][2]}; vector<int> tri1(tri1list, tri1list + sizeof(tri1list) / sizeof(int) );
                            int tri2list[] = {projlisti[j],trilist[k][1],trilist[k][2]}; vector<int> tri2(tri2list, tri2list + sizeof(tri2list) / sizeof(int) );
                            tripush.push_back(tri1);tripush.push_back(tri2);
                        }else if (area2<1e-5){
                            int tri1list[] = {trilist[k][1],projlisti[j],trilist[k][0]}; vector<int> tri1(tri1list, tri1list + sizeof(tri1list) / sizeof(int) );
                            int tri2list[] = {projlisti[j],trilist[k][2],trilist[k][0]}; vector<int> tri2(tri2list, tri2list + sizeof(tri2list) / sizeof(int) );
                            tripush.push_back(tri1);tripush.push_back(tri2);
                        }else{
                            int tri1list[] = {trilist[k][2],projlisti[j],trilist[k][1]}; vector<int> tri1(tri1list, tri1list + sizeof(tri1list) / sizeof(int) );
                            int tri2list[] = {projlisti[j],trilist[k][0],trilist[k][1]}; vector<int> tri2(tri2list, tri2list + sizeof(tri2list) / sizeof(int) );
                            tripush.push_back(tri1);tripush.push_back(tri2);
                        }
                        trilist.erase(trilist.begin()+k);
                        k -= 1;
                    }else{
                        // cout<<"on interior"<<'\n';
                        // projection locates on the interior
                        int tri1list[] = {trilist[k][0],projlisti[j],trilist[k][2]}; vector<int> tri1(tri1list, tri1list + sizeof(tri1list) / sizeof(int) );
                        int tri2list[] = {trilist[k][1],projlisti[j],trilist[k][0]}; vector<int> tri2(tri2list, tri2list + sizeof(tri2list) / sizeof(int) );
                        int tri3list[] = {trilist[k][2],projlisti[j],trilist[k][1]}; vector<int> tri3(tri3list, tri3list + sizeof(tri3list) / sizeof(int) );
                        trilist.push_back(tri1);trilist.push_back(tri2);trilist.push_back(tri3);
                        trilist.erase(trilist.begin()+k);
                        k -= 1;
                        break;
                    }
                }
            }
            
            for (k=0;k<tripush.size();k++){
                trilist.push_back(tripush[k]);
            }
        }
    }
    Reconnection reconnection(vertlist,trilist,projlist);
    return reconnection;
}
// Ning's code ends here
//----------------------------------------------------------------------------------------------------

set<int> cOctree::getListPolysToCheck(cLine &ray)
{
    // Returns a list of all polygons that are within OctNodes hit by a given ray
    set<int> intTestPolys;
    getPolysToCheck(root,ray,intTestPolys);
    return intTestPolys;
}

void cOctree::getPolysToCheck(cOctNode &node, cLine &ray, set<int> &intTestPolys)
{
    // Utility function for getListPolysToCheck. Finds all OctNodes hit by a given ray
    // and returns a list of the objects contained within
    if (node.sphereRayIntersect(ray)) {
        if (node.boxRayIntersect(ray)) {
            if (node.isLeafNode()) {
                for (int i=0; i<node.numPolys(); i++) {
                    intTestPolys.insert(node.data[i]); }
            } else {
                for (unsigned int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
                    getPolysToCheck(node.branches[i],ray,intTestPolys);
                } 
            }
        }
    }
}

vector<cOctNode*> cOctree::getSortedNodesToCheck(cLine &ray)
{
    // Finds all the nodes that intersect with given ray. Uses the nodes "position"
    // to sort the nodes by distance from the ray origin (in ascending order). 
    // Nodes that are closest to the ray origin will be checked first for poly 
    // intersections 
    vector<pair<cOctNode*,double> > nodeList;
    getNodesToCheck(root,ray,nodeList);
    sort(nodeList.begin(),nodeList.end(),sortNodes);
    vector<cOctNode*> nodeListSorted;
    vector<pair<cOctNode*,double> >::iterator it;
    for (it=nodeList.begin(); it!=nodeList.end(); it++) {
        nodeListSorted.push_back((*it).first); }
    return nodeListSorted;
}

void cOctree::getNodesToCheck(cOctNode &node, cLine &ray, vector<pair<cOctNode*,double> > &nodeList)
{
    // Utility function for getSortedNodesToCheck
    // Finds all the nodes that intersect with given ray. Projects the node "position" (node
    // centre) onto the ray to facilitate sorting of the nodes by distance from ray origin
    if (node.sphereRayIntersect(ray)) {
        if (node.boxRayIntersect(ray)) {
            if (node.isLeafNode()) {
                // Project node "position" on to ray and find distance from ray origin
                vector<double> oc = vectSubtract(node.position, ray.p0);
                double s = dotProduct(oc,ray.dir);
                // Add node and distance to list
                nodeList.push_back(std::make_pair(&node,s));
            } else {
                for (unsigned int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
                    getNodesToCheck(node.branches[i],ray,nodeList);
                } 
            }
        }
    }
}

vector<Intersection> cOctree::findRayIntersect(cLine &ray)
{   
    // Get polys to check
    set<int> polyListCheck = getListPolysToCheck(ray);
    
    // Loop through all polys in check list to find a possible intersection
    vector<Intersection> intersectList;
    set<int>::iterator it;
    vector<double> ip;
    double s;
    for (it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
        int polyLabel = *it;
        
        if (polyList[polyLabel].rayPlaneIntersectPoint(ray,ip,s)) {
            intersectList.push_back(Intersection(polyLabel,ip,s)); } 
    }
    //cout<<intersectList[1].triLabel<<'\n';
    if (intersectList.size()>0){
        // Sort list in terms of distance of the intersection from the ray origin
        sort(intersectList.begin(),intersectList.end());
        vector<Intersection> minintersect;
        // minintersect[0]=intersectList[0];
        minintersect.assign(intersectList.begin(),intersectList.begin()+1);
        // return intersectList;
        return minintersect;
    }else{
        Intersection mini;
        mini=findRaymindist(ray);
        //cout<<"mini.triLabel test2"<<'\n';
        //cout<<mini.triLabel<<'\n';
        vector<Intersection> minintersect;
        minintersect.push_back(mini);
        return minintersect;
    }
}

Intersection cOctree::findRaymindist(cLine &ray)
{
    // double distance = 1e99;
    double distanceij = 1e99;
    double area2;
    double h;
    vector<double> x1, x2, x0, x0x1, x0x2, x2x1, cross, PP;
    int k;
    for (unsigned int i=0;i!=polyList.size();i++){
        for (unsigned int j=0;j<3;j++){
            x1.assign(ray.p0.begin(),ray.p0.end());
            x2.assign(ray.p1.begin(),ray.p1.end());
            x0.assign(polyList[i].vertices[j].begin(),polyList[i].vertices[j].end());
            x0x1=vectSubtract(x0,x1);
            x0x2=vectSubtract(x0,x2);
            cross=crossProduct(x0x1,x0x2);
            area2=sqrt(dotProduct(cross,cross));
            x2x1=vectSubtract(x2,x1);
            h=area2/sqrt(dotProduct(x2x1,x2x1));
            if (distanceij>h){
                PP.assign(x0.begin(),x0.end());
                distanceij=h;    
                k=i;
            }
        }
    }
    Intersection intersect(k,PP,distanceij);
    return intersect;
}
// ------------------------------------------------------
// Intersection cOctree::findRayIntersect(cLine &ray)
// {   
//     // Get polys to check
//     set<int> polyListCheck = getListPolysToCheck(ray);
    
//     // Loop through all polys in check list to find a possible intersection
//     vector<Intersection> intersectList;
//     set<int>::iterator it;
//     vector<double> ip;
//     double s;
//     for (it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
//         int polyLabel = *it;
//         if (polyList[polyLabel].rayPlaneIntersectPoint(ray,ip,s)) {
//             intersectList.push_back(Intersection(polyLabel,ip,s)); } 
//     }
    
//     // Sort list in terms of distance of the intersection from the ray origin
//     sort(intersectList.begin(),intersectList.end());
//     Intersection intersection = intersectList[0];
//     return intersection;
// }


vector<Intersection> cOctree::findRayIntersects(vector<cLine> &rayList)
{
    // For each ray provided, determines if ray hits a poly in the tree and 
    // returns a boolean integer. Uses openmp to speed up the calculation
    // Function findRayIntersectsSorted is a similar function that sorts the
    // triangles in order of closest octNodes. For ray casting, this alternative
    // function should be faster in most cases
    
    int numRays = (int)(rayList.size());
    vector<Intersection> intersections;
    #pragma omp parallel for
    // cout<<"numRays ";
    // cout<<numRays<<'\n';
    for (int i=0; i<numRays; i++) 
    {
        cLine *ray = &rayList[i]; 
        set<int> polyListCheck = getListPolysToCheck(*ray);
        vector<Intersection> intersection = findRayIntersect(*ray);
        intersections.push_back(intersection[0]);
        //cout<<"intersection[0].triLabel"<<'\n';
        //cout<<intersection[0].triLabel<<'\n';
    }
    return intersections;
}



// ------------------------------------------------------


// vector<int> cOctree::findRayIntersects(vector<cLine> &rayList)
// {
//     // For each ray provided, determines if ray hits a poly in the tree and 
//     // returns a boolean integer. Uses openmp to speed up the calculation
//     // Function findRayIntersectsSorted is a similar function that sorts the
//     // triangles in order of closest octNodes. For ray casting, this alternative
//     // function should be faster in most cases
    
//     int numRays = (int)(rayList.size());
//     vector<int> foundIntsects(numRays,0);
//     #pragma omp parallel for
//     for (int i=0; i<numRays; i++) 
//     {
//         cLine *ray = &rayList[i]; 
//         set<int> polyListCheck = getListPolysToCheck(*ray);
//         for (set<int>::iterator it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
//             int polyLabel = *it;
//             if (polyList[polyLabel].rayPlaneIntersectPoint(*ray)) {
//                 foundIntsects[i] = 1; break; } 
//         }
//     }
//     return foundIntsects;
// }

vector<int> cOctree::findRayIntersectsSorted(vector<cLine> &rayList)
{
    // For each ray provided, determines if ray hits a poly in the tree and 
    // returns a boolean integer. 
    // Uses "getSortedNodesToCheck", which returns a list of nodes that intersect
    // with the given ray sorted in ascending order of distance from the ray origin
    // Uses openmp to speed up the calculation
    
    int numRays = (int)(rayList.size());
    vector<int> foundIntsects(numRays,0);
    
    #pragma omp parallel for
    for (int i=0; i<numRays; i++) 
    {
        // Get branches to check. Branches are sorted in ascending distance 
        // from ray origin
        cLine *ray = &rayList[i]; 
        vector<cOctNode*> nodeList = getSortedNodesToCheck(*ray);
        
        // Loop through sorted branches, checking the polys contained within each
        vector<cOctNode*>::iterator it;
        for (it=nodeList.begin(); it!=nodeList.end(); ++it) {
            cOctNode *node = *it;
            for (unsigned int j=0; j<node->data.size(); j++) {
                int polyLabel = node->data[j];
                if (polyList[polyLabel].rayPlaneIntersectPoint(*ray)) {
                    foundIntsects[i]=1; break; } 
            }
            // If any poly from current node is hit, proceed on to the next node 
            if (foundIntsects[i]==1) break;
        } 
    }
    return foundIntsects;
}

// ------------------------------------------------------
double vectNorm(vector<double> &a)
{
    double norm = sqrt(dotProduct(a,a));
    return norm;
}

vector<double> normalize(vector<double> &v)
{   double sum=0.0;
    for (unsigned int i=0;i<3;i++){
        sum+=v[i]*v[i];
    }
    sum=sqrt(sum);
    vector<double> norm(3);
    for (unsigned int i=0;i<3;i++){
        norm[i]=v[i]/sum;
    }
    return norm;
}


double dotProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates dot product v1.v2
    double dp=0.0;
    for (unsigned int i=0; i<v1.size(); i++)
        dp += v1[i]*v2[i]; 
    return dp;
}

double distBetweenPoints(vector<double> &p1, vector<double> &p2)
{
    // Calculate the distance between points p1 and p2, |p1-p2|
    double sum=0.0;
    for (unsigned int i=0; i<3; i++)
        sum += pow((p1[i]-p2[i]),2.0);
    sum = sqrt(sum);
    return sum;
}

vector<double> crossProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates cross product v1xv2
    vector<double> cp(3);
    cp[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cp[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cp[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cp;
}

vector<double> vectAdd( vector<double> &a, vector<double> &b )
{
    // Vector addition, c=a+b
    return vectAdd(a, b, 1.0);
}

vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf )
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i] + sf*b[i];
    return c;
}
vector<double> vectSubtract( vector<double> &a, vector<double> &b )
{
    // Vector subtraction, c=a-b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i]-b[i];
    return c;
}
// ---------------------------------------------
vector<double> vectMuldou( vector<double> &a, double sf)
{
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i]*sf;
    return c;
}

// ---------------------------------------------
string NumberToString( int Number )
{
    // Converts integer to string
    ostringstream ss;
    ss << Number;
    return ss.str();
}

bool sortNodes(const pair<cOctNode*,double>&i, const pair<cOctNode*,double>&j)
{
    // Function used to sort a vector of cOctnode,double pairs by the value of
    // the double. The double will typically represent distance from the ray
    // origin in a ray-node intersection test
    return i.second < j.second;
}

vector<vector<double> > vectMultiply(vector<vector<double> > m1, vector<vector<double> > m2)
{
    vector<vector<double> > res(m1.size(),vector<double>(m2[0].size(),0.0));
    for (unsigned int i=0;i<m1.size();i++){
        for (unsigned int j=0;j<m2[0].size();j++){
            for (unsigned int k=0;k<m2.size();k++){
                res[i][j] += m1[i][k]*m2[k][j];
            }
        }
    }
    return res;
}

vector<vector<double> > vectInverse(vector<vector<double> > &m)
{
    double det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
    vector<vector<double> > w(2,vector<double>(2,0.0));
    w[0][0] = m[1][1]/det;
    w[0][1] = -m[0][1]/det;
    w[1][0] = -m[1][0]/det;
    w[1][1] = m[0][0]/det;
    return w;
}

vector<vector<double> > vectneg(vector<vector<double> > m)
{
    vector<vector<double> > res(m.size(),vector<double>(m[0].size(),0.0));
    for (unsigned int i=0;i<m.size();i++){
        for (unsigned int j=0;j<m[0].size();j++){
            res[i][j] = -m[i][j];
        }
    }
    return res;
}

vector<double> vectneg(vector<double> &m)
{
    vector<double> res(m.size(),0.0);
    for (unsigned int i=0;i<m.size();i++){
            res[i] = -m[i];
    }
    return res;
}

vector<double> midPoint(vector<double> &a, vector<double> &b)
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> m(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        m[i] = 0.5*a[i] + 0.5*b[i];
    return m;
}

double angleOfVecs(vector<double> &a, vector<double> &b)
{
    double PI = atan(1)*4;
    double angle = acos(dotProduct(a,b)/(vectNorm(a)*vectNorm(b)))/PI*180;
    return angle;
} 
// ------------------------------------------------------
