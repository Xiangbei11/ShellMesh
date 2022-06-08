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
using namespace std;

typedef struct equalityconstraint
{
    vector<int> data;
    vector<int> row;
    vector<int> col;
    vector<int> b;
    int sharededges;
    equalityconstraint() { data.resize(3,0); row.resize(3,0.0); col.resize(3,0.0); b.resize(3,0.0); sharededges=3;}
    equalityconstraint(vector<int> _data, vector<int> _row, vector<int> _col, vector<int> _b, int _sharededges) { data=_data; row=_row; col=_col; b=_b; sharededges=_sharededges;}
} EqualityConstraint;

typedef struct indices
{
    vector<vector<int> > triindex;
    vector<int> edgeindex;
    vector<vector<int> > verindex;
    indices() {triindex.resize(3,vector<int>(2)); edgeindex.resize(3,0.0); verindex.resize(3,vector<int>(2));}
    indices(vector<vector<int> > _triindex, vector<int> _edgeindex, vector<vector<int> > _verindex) {triindex=_triindex; edgeindex=_edgeindex; verindex=_verindex;}
}Indices;

typedef struct updatemesh
{
    vector<vector<double> > splitcoord;
    vector<vector<int> > triadd;
    vector<vector<int> > quadadd;
    vector<int> tridelete;
    vector<int> quaddelete;
    vector<vector<double> > vertlist;
    updatemesh() {splitcoord.resize(3,vector<double>(3));triadd.resize(3,vector<int>(3));quadadd.resize(4,vector<int>(3));tridelete.resize(3,1);quaddelete.resize(3,1);vertlist.resize(3,vector<double>(3));}
    updatemesh(vector<vector<double> > _splitcoord, vector<vector<int> > _triadd, vector<vector<int> > _quadadd, vector<int> _tridelete, vector<int> _quaddelete, vector<vector<double> > _vertlist)
    {splitcoord=_splitcoord; triadd=_triadd; quadadd=_quadadd; tridelete=_tridelete; quaddelete=_quaddelete; vertlist=_vertlist;} 
}UpdateMesh;

typedef struct fullyquad
{
    vector<vector<int> > quadlist;
    vector<vector<double> > vertlist;
    vector<int> fixedvert;
    fullyquad() {quadlist.resize(3,vector<int>(2)); vertlist.resize(3,vector<double>(2)); fixedvert.resize(3,0);}
    fullyquad(vector<vector<int> > _quadlist, vector<vector<double> > _vertlist, vector<int> _fixedvert) {quadlist = _quadlist; vertlist = _vertlist; fixedvert = _fixedvert;}
}FullyQuad;

typedef struct splitoutput
{
    vector<vector<double> > vertlist;
    vector<vector<double> > splitcoord;
    vector<vector<int> > trilist;
    vector<vector<int> > quadlist;
    vector<int> fixedvert;
    splitoutput() {vertlist.resize(3,vector<double>(3));splitcoord.resize(3,vector<double>(3));trilist.resize(3,vector<int>(3));quadlist.resize(3,vector<int>(4));fixedvert.resize(3,0);}
    splitoutput(vector<vector<double> > _vertlist, vector<vector<double> > _splitcoord, vector<vector<int> > _trilist, vector<vector<int> > _quadlist, vector<int> _fixedvert)
    {vertlist=_vertlist; splitcoord=_splitcoord; trilist=_trilist; quadlist=_quadlist; fixedvert = _fixedvert;}
}SplitOutput;

typedef struct mergeoutput
{
    vector<vector<int> > trilist;
    vector<vector<int> > quadlist;
    mergeoutput() {trilist.resize(3,vector<int>(3));quadlist.resize(3,vector<int>(4));}
    mergeoutput(vector<vector<int> > _trilist, vector<vector<int> > _quadlist){trilist=_trilist; quadlist=_quadlist;}
}MergeOutput;

typedef struct optprep
{
    vector<vector<double> > normal;
    vector<double> weight;
    optprep() {normal.resize(3,vector<double>(3,0.0));weight.resize(3,0.0);}
    optprep(vector<vector<double> > _normal, vector<double> _weight){normal=_normal; weight=_weight;}
}OptPrep;

typedef struct triint
{
    vector<vector<int> > itris;
    vector<vector<int> > jtris;
    vector<vector<double> > icoords;
    vector<vector<double> > jcoords;
    vector<int> ifixed;
    vector<int> jfixed;
    triint() {itris.resize(3,vector<int>(3,0)); jtris.resize(3,vector<int>(3,0)); icoords.resize(3,vector<double>(3,0.0)); jcoords.resize(3,vector<double>(3,0.0)); ifixed.resize(3,0); jfixed.resize(3,0);}
    triint(vector<vector<int> > _itris, vector<vector<int> > _jtris, vector<vector<double> > _icoords, vector<vector<double> > _jcoords, vector<int> _ifixed, vector<int> _jfixed) {itris = _itris; jtris = _jtris; icoords = _icoords; jcoords = _jcoords; ifixed = _ifixed; jfixed = _jfixed;}
}TriInt;

typedef struct reconnection
{
    vector<vector<double> > vertlist;
    vector<vector<int> > trilist;
    vector<int> projlist;
    reconnection() { vertlist.resize(3,vector<double>(3,0.0)); trilist.resize(3,vector<int>(3,0)); projlist.resize(3,0);}
    reconnection(vector<vector<double> > _vertlist, vector<vector<int> > _trilist, vector<int> _projlist){ vertlist = _vertlist; trilist = _trilist; projlist = _projlist;}
}Reconnection;

typedef struct uni
{
    vector<vector<double> > vertlist;
    vector<vector<int> > quadlist;
    uni() {vertlist.resize(3,vector<double>(3,0.0)); quadlist.resize(3,vector<int>(4,0));}
    uni(vector<vector<double> > _vertlist, vector<vector<int> > _quadlist){ vertlist = _vertlist; quadlist = _quadlist;}
}Uni;

class opt{
public:
    
    vector<vector<double> > vertexCoords;
    vector<vector<int> > triConnectivity;
    vector<vector<int> > quadConnectivity;
    double w1;
    double w2;
    double w3;
    const double p = 2.0;
    const double theta0 = 60.0;// reference angle from equilateral triangles
    const double theta1 = 90.0;// reference angle from square
    opt(vector<vector<double> > _vertexCoords, vector<vector<int> > _triConnectivity, vector<vector<int> > _quadConnectivity, double w1, double w2, double w3);
    ~opt();
    vector<double> splittriangle(vector<vector<double> > vertices);//, double w1, double w2);
    vector<double> splitquad(vector<vector<double> > vertices);//, double w1, double w2);
    vector<vector<double> > buildupmatrix();
    vector<vector<double> > buildbounds(vector<int> fixedvert);
    vector<vector<int> > countedges();
    vector<vector<int> > countedges(vector<vector<int> > quadlist);
    EqualityConstraint equcon();
    Indices counttrisharededges();
    vector<double> mergetriangles();
    vector<vector<int> > definetriedges();
    vector<vector<int> > inequncon();
    vector<vector<double> > mergingbounds(vector<int> fixedvert);
    FullyQuad fullyquadmesh(vector<int> fixed);
    OptPrep computenormal();
    double computeweight(vector<vector<double> > directions);
    vector<vector<vector<double> > > computecs();
    SplitOutput splitupdatemesh(vector<int> optres, vector<int> fixed);
    UpdateMesh updatetrisplit(vector<int> optres, vector<vector<int> > trilist, vector<vector<double> > vertlist);
    UpdateMesh updatequadsplit(vector<int> optres, vector<vector<int> > quadlist, vector<vector<double> > vertlist);
    MergeOutput mergeupdatemesh(vector<int> optres);
    vector<vector<int> > connectTFI(vector<vector<double> > boundaries, vector<vector<double> > vertices, vector<int> num_b, vector<int> num_v, vector<int> index_b, vector<int> index_v);
    double averagesize();
    TriInt member_intersection(vector<vector<int> > itris, vector<vector<int> > jtris, vector<vector<double> > icoords, vector<vector<double> > jcoords, vector<int> ifixed, vector<int> jfixed);
    vector<double> intersectionpoint(vector<double> q1, vector<double> q2, vector<double> p1, vector<double> p2, vector<double> p3);
    vector<double> computeradius();
    vector<vector<double> > sortshared(vector<vector<double> > shared);
    vector<vector<double> > uniquevet(vector<vector<double> > sorted);
    Reconnection reconnect(vector<vector<double> > sorted, vector<vector<int> > tris, vector<vector<double> > coords, vector<int> fixed);
    vector<vector<int> > counttriedges(vector<vector<int> > trilist);
    vector<int> removeduplicatetris();
    Uni mergeduppts(double tol);
};
vector<double> vectSubtract( vector<double> &a, vector<double> &b );
vector<double> midPoint( vector<double> &a, vector<double> &b );
double angleOfVecs(vector<double> &a, vector<double> &b);
double angleRatio(vector<double> &a, vector<double> &b, double theta0);
double dotProduct(vector<double> &a, vector<double> &b);
double vectNorm(vector<double> &a);
void arrayadd(int arr[], int size, int val);
vector<double> crossProduct(vector<double> &v1, vector<double> &v2);
vector<double> normalize(vector<double> &v);
double distBetweenPoints(vector<double> &p1, vector<double> &p2);
vector<double> midPoint(vector<double> &a, vector<double> &b, vector<double> &c);
vector<double> midPoint(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d);
vector<double> vectScale(vector<double> &p, double s);
vector<double> vectMean(vector<vector<double> > p);
vector<double> vectAdd( vector<double> &a, vector<double> &b );
vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf );
double signedvolume(vector<double> a, vector<double> b, vector<double> c, vector<double> d);