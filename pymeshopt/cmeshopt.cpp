#include "cmeshopt.h"


opt::opt(vector<vector<double> > _vertexCoords, vector<vector<int> > _triConnectivity, vector<vector<int> > _quadConnectivity, double _w1, double _w2, double _w3)
{
    vertexCoords = _vertexCoords;
    triConnectivity = _triConnectivity;
    quadConnectivity = _quadConnectivity;
    w1 = _w1;
    w2 = _w2;
    w3 = _w3;
}

opt::~opt() {}

// bulid up coefficient matrix for splitting optimization
vector<double> opt::calculate_aspect_ratio()
{       
    int i;
    if (triConnectivity.size() != 0) {
        cout<<"Not fully quad"<<"\n";        
    }
    int num_poly = quadConnectivity.size(); 
    vector<double> aspect_ratio(num_poly,0.0);
    vector<double> verinit(3,0.0);
    vector<vector<double> > quadvertices(4,verinit);

    for (i=0;i<(int) quadConnectivity.size();i++){
        quadvertices[0] = vertexCoords[quadConnectivity[i][0]];
        quadvertices[1] = vertexCoords[quadConnectivity[i][1]];
        quadvertices[2] = vertexCoords[quadConnectivity[i][2]];
        quadvertices[3] = vertexCoords[quadConnectivity[i][3]];
        vector<double> e01 = vectSubtract(quadvertices[1],quadvertices[0]);
        vector<double> e12 = vectSubtract(quadvertices[2],quadvertices[1]);
        vector<double> e23 = vectSubtract(quadvertices[3],quadvertices[2]);
        vector<double> e30 = vectSubtract(quadvertices[0],quadvertices[3]);
        double opt10Lavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(e23)+vectNorm(e30))/4;
        double shortest = vectNorm(e01);
        double longest = vectNorm(e01);
        if (vectNorm(e12)>longest){longest = vectNorm(e12);}
        if (vectNorm(e23)>longest){longest = vectNorm(e23);}
        if (vectNorm(e30)>longest){longest = vectNorm(e30);}
        if (vectNorm(e12)<shortest){shortest = vectNorm(e12);}
        if (vectNorm(e23)<shortest){shortest = vectNorm(e23);}
        if (vectNorm(e30)<shortest){shortest = vectNorm(e30);}        
        //aspect_ratio[i] = (pow(vectNorm(e01)/opt10Lavg,p) + pow(vectNorm(e12)/opt10Lavg,p)
        //    + pow(vectNorm(e23)/opt10Lavg,p) + pow(vectNorm(e30)/opt10Lavg,p))/4;
        aspect_ratio[i] = longest/shortest;
    }
    return aspect_ratio;    
}

// bulid up coefficient matrix for splitting optimization
vector<vector<double> > opt::calculate_internal_angles()
{   
    int i;
    if (triConnectivity.size() != 0) {
        cout<<"Not fully quad"<<"\n";        
    }
    int num_poly = quadConnectivity.size(); 
    vector<double> init(4,0.0);
    vector<vector<double> > internal_angles(num_poly, init);
    vector<double> verinit(3,0.0);
    vector<vector<double> > quadvertices(4,verinit);
    
    int theta = 1;
    for (i=0;i<(int) quadConnectivity.size();i++){
        quadvertices[0] = vertexCoords[quadConnectivity[i][0]];
        quadvertices[1] = vertexCoords[quadConnectivity[i][1]];
        quadvertices[2] = vertexCoords[quadConnectivity[i][2]];
        quadvertices[3] = vertexCoords[quadConnectivity[i][3]];
        vector<double> angles(4,0.0);
        vector<double> e01 = vectSubtract(quadvertices[1],quadvertices[0]);
        vector<double> e03 = vectSubtract(quadvertices[3],quadvertices[0]);
        vector<double> e10 = vectSubtract(quadvertices[0],quadvertices[1]);
        vector<double> e12 = vectSubtract(quadvertices[2],quadvertices[1]);
        vector<double> e21 = vectSubtract(quadvertices[1],quadvertices[2]);
        vector<double> e23 = vectSubtract(quadvertices[3],quadvertices[2]);
        vector<double> e32 = vectSubtract(quadvertices[2],quadvertices[3]);
        vector<double> e30 = vectSubtract(quadvertices[0],quadvertices[3]);        
        angles[0] = angleRatio(e01,e03,theta);
        angles[1] = angleRatio(e10,e12,theta);
        angles[2] = angleRatio(e21,e23,theta);
        angles[3] = angleRatio(e32,e30,theta);
        internal_angles[i] = angles;
    }
    return internal_angles;
    
}


// bulid up coefficient matrix for splitting optimization
vector<vector<double> > opt::buildupmatrix()
{   
    
    int i, j;
    int num_poly = triConnectivity.size() + quadConnectivity.size(); 
    vector<double> init(11,0.0);
    vector<vector<double> > matrix(num_poly, init);
    vector<double> verinit(3,0.0);
    vector<vector<double> > trivertices(3,verinit);
    vector<vector<double> > quadvertices(4,verinit);
    
    // loop through every triangle and quad
    for (i=0;i<(int) triConnectivity.size();i++){
        trivertices[0] = vertexCoords[triConnectivity[i][0]];
        trivertices[1] = vertexCoords[triConnectivity[i][1]];
        trivertices[2] = vertexCoords[triConnectivity[i][2]];
        
        vector<double> abc = splittriangle(trivertices);//Xb: can be comment
        matrix[i] = splittriangle(trivertices);
    }
    for (j=0;j<(int) quadConnectivity.size();j++){
        quadvertices[0] = vertexCoords[quadConnectivity[j][0]];
        quadvertices[1] = vertexCoords[quadConnectivity[j][1]];
        quadvertices[2] = vertexCoords[quadConnectivity[j][2]];
        quadvertices[3] = vertexCoords[quadConnectivity[j][3]];
        matrix[triConnectivity.size()+j] = splitquad(quadvertices); 
    }
    return matrix;
    
}

// 7 splitting options for triangles and their corresponding energy
vector<double> opt::splittriangle(vector<vector<double> > vertices)
{   
    // six vectors starting from each vertice to the other two
    vector<double> e01 = vectSubtract(vertices[1],vertices[0]);
    vector<double> e02 = vectSubtract(vertices[2],vertices[0]);
    vector<double> e10 = vectSubtract(vertices[0],vertices[1]);
    vector<double> e12 = vectSubtract(vertices[2],vertices[1]);
    vector<double> e20 = vectSubtract(vertices[0],vertices[2]);
    vector<double> e21 = vectSubtract(vertices[1],vertices[2]);
    // mid point of each edge
    vector<double> m01 = midPoint(vertices[0],vertices[1]);
    vector<double> m12 = midPoint(vertices[1],vertices[2]);
    vector<double> m20 = midPoint(vertices[2],vertices[0]);
    vector<double> options(11,0.0);// initialize the options vector
    // option 0 
    vector<double> v0m12 = vectSubtract(m12,vertices[0]);
    vector<double> m12v0 = vectSubtract(vertices[0],m12);
    vector<double> v1m12 = vectSubtract(m12,vertices[1]);
    vector<double> m12v1 = vectSubtract(vertices[1],m12);
    vector<double> v2m12 = vectSubtract(m12,vertices[2]);
    vector<double> m12v2 = vectSubtract(vertices[2],m12);
    double opt0Lavg = (vectNorm(e02)+vectNorm(e01)+vectNorm(v0m12)+vectNorm(v1m12)+vectNorm(v2m12))/5;

     
    //double p = 2.0;
    options[0] = w1*(pow(angleRatio(e02,v0m12,theta0),p) + pow(angleRatio(e01,v0m12,theta0),p)
               + pow(angleRatio(e10,v1m12,theta0),p) + pow(angleRatio(e20,v2m12,theta0),p)
               + pow(angleRatio(m12v2,m12v0,theta0),p) + pow(angleRatio(m12v1,m12v0,theta0),p))/6
               + w2*(pow(vectNorm(e02)/opt0Lavg,p) + pow(vectNorm(e01)/opt0Lavg,p)
               + pow(vectNorm(v0m12)/opt0Lavg,p) + pow(vectNorm(v1m12)/opt0Lavg,p)
               + pow(vectNorm(v2m12)/opt0Lavg,p))/5;
    // option 1
    vector<double> v0m20 = vectSubtract(m20,vertices[0]);
    vector<double> m20v0 = vectSubtract(vertices[0],m20);
    vector<double> v1m20 = vectSubtract(m20,vertices[1]);
    vector<double> m20v1 = vectSubtract(vertices[1],m20);
    vector<double> v2m20 = vectSubtract(m20,vertices[2]);
    vector<double> m20v2 = vectSubtract(vertices[2],m20);
    double opt1Lavg = (vectNorm(e10)+vectNorm(e12)+vectNorm(v0m20)+vectNorm(v1m20)+vectNorm(v2m20))/5;
    options[1] = w1*(pow(angleRatio(e10,v1m20,theta0),p) + pow(angleRatio(e12,v1m20,theta0),p)
               + pow(angleRatio(e21,v2m20,theta0),p) + pow(angleRatio(e01,v0m20,theta0),p)
               + pow(angleRatio(m20v1,m20v2,theta0),p) + pow(angleRatio(m20v1,m20v0,theta0),p))/6
               + w2*(pow(vectNorm(e10)/opt1Lavg,p) + pow(vectNorm(e12)/opt1Lavg,p)
               + pow(vectNorm(v0m20)/opt1Lavg,p) + pow(vectNorm(v1m20)/opt1Lavg,p)
               + pow(vectNorm(v2m20)/opt1Lavg,p))/5;
    // option 2
    vector<double> v0m01 = vectSubtract(m01,vertices[0]);
    vector<double> m01v0 = vectSubtract(vertices[0],m01);
    vector<double> v1m01 = vectSubtract(m01,vertices[1]);
    vector<double> v2m01 = vectSubtract(m01,vertices[2]);
    vector<double> m01v1 = vectSubtract(vertices[1],m01);
    vector<double> m01v2 = vectSubtract(vertices[2],m01);
    double opt2Lavg = (vectNorm(e21)+vectNorm(e20)+vectNorm(v0m01)+vectNorm(v1m01)+vectNorm(v2m01))/5;
    options[2] = w1*(pow(angleRatio(e21,v2m01,theta0),p) + pow(angleRatio(e20,v2m01,theta0),p)
               + pow(angleRatio(e02,v0m01,theta0),p) + pow(angleRatio(e12,v1m01,theta0),p)
               + pow(angleRatio(m01v2,m01v0,theta0),p) + pow(angleRatio(m01v2,m01v1,theta0),p))/6
               + w2*(pow(vectNorm(e21)/opt2Lavg,p) + pow(vectNorm(e20)/opt2Lavg,p)
               + pow(vectNorm(v0m01)/opt2Lavg,p) + pow(vectNorm(v1m01)/opt2Lavg,p)
               + pow(vectNorm(v2m01)/opt2Lavg,p))/5;
    // option 3
    vector<double> m01m12 = vectSubtract(m12,m01);
    vector<double> m12m01 = vectSubtract(m01,m12);
    double opt3Lavg = (vectNorm(e02)+vectNorm(v0m01)+vectNorm(m01v1)+vectNorm(v1m12)+vectNorm(m12v2)+vectNorm(m01m12))/6;
    options[3] = w1*(pow(angleRatio(e02,v0m01,theta1),p) + pow(angleRatio(m01v0,m01m12,theta1),p)
               + pow(angleRatio(m12m01,m12v2,theta1),p) + pow(angleRatio(v2m12,e20,theta1),p)// quad reference angle theta1
               + pow(angleRatio(m01m12,m01v1,theta0),p) + pow(angleRatio(v1m01,v1m12,theta0),p) + pow(angleRatio(m12v1,m12m01,theta0),p))/7
               + w2*(pow(vectNorm(e02)/opt3Lavg,p) + pow(vectNorm(v0m01)/opt3Lavg,p)
               + pow(vectNorm(m01v1)/opt3Lavg,p) + pow(vectNorm(v1m12)/opt3Lavg,p)
               + pow(vectNorm(m12v2)/opt3Lavg,p) + pow(vectNorm(m01m12)/opt3Lavg,p))/6;
    // option 4
    vector<double> m12m20 = vectSubtract(m20,m12);
    vector<double> m20m12 = vectSubtract(m12,m20);
    double opt4Lavg = (vectNorm(e10)+vectNorm(v1m12)+vectNorm(m12v2)+vectNorm(v2m20)+vectNorm(m20v0)+vectNorm(m12m20))/6;
    options[4] =  w1*(pow(angleRatio(e10,v1m12,theta1),p) + pow(angleRatio(m12v1,m12m20,theta1),p)
               + pow(angleRatio(m20m12,m20v0,theta1),p) + pow(angleRatio(v0m20,e01,theta1),p)// quad reference angle theta1
               + pow(angleRatio(m12m20,m12v2,theta0),p) + pow(angleRatio(v2m12,v2m20,theta0),p) + pow(angleRatio(m20v2,m20m12,theta0),p))/7
               + w2*(pow(vectNorm(e10)/opt4Lavg,p) + pow(vectNorm(v1m12)/opt4Lavg,p)
               + pow(vectNorm(m12v2)/opt4Lavg,p) + pow(vectNorm(v2m20)/opt4Lavg,p)
               + pow(vectNorm(m20v0)/opt4Lavg,p) + pow(vectNorm(m12m20)/opt4Lavg,p))/6;
    // option 5
    vector<double> m20m01 = vectSubtract(m01,m20);
    vector<double> m01m20 = vectSubtract(m20,m01);
    double opt5Lavg = (vectNorm(e21)+vectNorm(v2m20)+vectNorm(m20v0)+vectNorm(v0m01)+vectNorm(m01v1)+vectNorm(m20m01))/6;
    options[5] =  w1*(pow(angleRatio(e21,v2m20,theta1),p) + pow(angleRatio(m20v2,m20m01,theta1),p)
               + pow(angleRatio(m01m20,m01v1,theta1),p) + pow(angleRatio(v1m01,e12,theta1),p)// quad reference angle theta1
               + pow(angleRatio(m20m01,m20v0,theta0),p) + pow(angleRatio(v0m20,v0m01,theta0),p) + pow(angleRatio(m01v0,m01m20,theta0),p))/7
               + w2*(pow(vectNorm(e21)/opt5Lavg,p) + pow(vectNorm(v2m20)/opt5Lavg,p)
               + pow(vectNorm(m20v0)/opt5Lavg,p) + pow(vectNorm(v0m01)/opt5Lavg,p)
               + pow(vectNorm(m01v1)/opt5Lavg,p) + pow(vectNorm(m20m01)/opt5Lavg,p))/6;
    // option 6
    double opt6Lavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(e20))/3;
    options[6] = w1*(pow(angleRatio(e01,e02,theta0),p) + pow(angleRatio(e10,e12,theta0),p)
               + pow(angleRatio(e20,e21,theta0),p))/3
               + w2*(pow(vectNorm(e01)/opt6Lavg,p) + pow(vectNorm(e12)/opt6Lavg,p)
               + pow(vectNorm(e20)/opt6Lavg,p))/3;
    return options;
}

vector<double> opt::splitquad(vector<vector<double> > vertices)
{
    // p-norm constant
    // double p = 2.0;
    // eight vectors starting from each vertice to the adjacent two vertices
    vector<double> e01 = vectSubtract(vertices[1],vertices[0]);
    vector<double> e03 = vectSubtract(vertices[3],vertices[0]);
    vector<double> e10 = vectSubtract(vertices[0],vertices[1]);
    vector<double> e12 = vectSubtract(vertices[2],vertices[1]);
    vector<double> e21 = vectSubtract(vertices[1],vertices[2]);
    vector<double> e23 = vectSubtract(vertices[3],vertices[2]);
    vector<double> e32 = vectSubtract(vertices[2],vertices[3]);
    vector<double> e30 = vectSubtract(vertices[0],vertices[3]);
    // mid point of each edge
    vector<double> m01 = midPoint(vertices[0],vertices[1]);
    vector<double> m12 = midPoint(vertices[1],vertices[2]);
    vector<double> m23 = midPoint(vertices[2],vertices[3]);
    vector<double> m30 = midPoint(vertices[3],vertices[0]);
    // third energy involving dihedral angle
    vector<double> edgevector = vectSubtract(vertices[2],vertices[0]);
    vector<double> diagvector = vectSubtract(vertices[3],vertices[1]);
    vector<double> cross = crossProduct(edgevector,diagvector);
    cross = normalize(cross);
    vector<double> e01n = normalize(e01), e12n = normalize(e12), e23n = normalize(e23), e30n = normalize(e30);



    double theta0 = 60.0;// reference angle from equilateral triangles
    double theta1 = 90.0;// reference angle from square
    vector<double> options(11,0.0);// initialize the options vector
    // option 0
    vector<double> v0m01 = vectSubtract(m01,vertices[0]);
    vector<double> m01v0 = vectSubtract(vertices[0],m01);
    vector<double> v1m01 = vectSubtract(m01,vertices[1]);
    vector<double> m01v1 = vectSubtract(vertices[1],m01);
    vector<double> v2m23 = vectSubtract(m23,vertices[2]);
    vector<double> m23v2 = vectSubtract(vertices[2],m23);
    vector<double> v3m23 = vectSubtract(m23,vertices[3]);
    vector<double> m23v3 = vectSubtract(vertices[3],m23);
    vector<double> m01m23 = vectSubtract(m23,m01);
    vector<double> m23m01 = vectSubtract(m01,m23);
    double opt0Lavg = (vectNorm(v0m01)+vectNorm(m01v1)+vectNorm(e12)+vectNorm(v2m23)+vectNorm(m23v3)+vectNorm(e30)+vectNorm(m01m23))/7;
    options[0] = w1*(pow(angleRatio(e03,v0m01,theta1),p) + pow(angleRatio(v1m01,e12,theta1),p)
               + pow(angleRatio(e21,v2m23,theta1),p) + pow(angleRatio(v3m23,e30,theta1),p)
               + pow(angleRatio(m01v0,m01m23,theta1),p) + pow(angleRatio(m01m23,m01v1,theta1),p)
               + pow(angleRatio(m23v2,m23m01,theta1),p) + pow(angleRatio(m23m01,m23v3,theta1),p))/8
               + w2*(pow(vectNorm(v0m01)/opt0Lavg,p) + pow(vectNorm(m01v1)/opt0Lavg,p)
               + pow(vectNorm(e12)/opt0Lavg,p) + pow(vectNorm(v2m23)/opt0Lavg,p)
               + pow(vectNorm(m23v3)/opt0Lavg,p) + pow(vectNorm(e30)/opt0Lavg,p)
               + pow(vectNorm(m01m23)/opt0Lavg,p))/7;
    // option 1
    vector<double> v0m30 = vectSubtract(m30,vertices[0]);
    vector<double> m30v0 = vectSubtract(vertices[0],m30);
    vector<double> v3m30 = vectSubtract(m30,vertices[3]);
    vector<double> m30v3 = vectSubtract(vertices[3],m30);
    vector<double> v1m12 = vectSubtract(m12,vertices[1]);
    vector<double> m12v1 = vectSubtract(vertices[1],m12);
    vector<double> v2m12 = vectSubtract(m12,vertices[2]);
    vector<double> m12v2 = vectSubtract(vertices[2],m12);
    vector<double> m12m30 = vectSubtract(m30,m12);
    vector<double> m30m12 = vectSubtract(m12,m30);   
    double opt1Lavg = (vectNorm(e01)+vectNorm(v1m12)+vectNorm(m12v2)+vectNorm(e23)+vectNorm(v3m30)+vectNorm(m30v0)+vectNorm(m12m30))/7;
    options[1] = w1*(pow(angleRatio(v0m30,e01,theta1),p) + pow(angleRatio(e10,v1m12,theta1),p)
               + pow(angleRatio(v2m12,e23,theta1),p) + pow(angleRatio(e32,v3m30,theta1),p)
               + pow(angleRatio(m12v1,m12m30,theta1),p) + pow(angleRatio(m12m30,m12v2,theta1),p)
               + pow(angleRatio(m30v3,m30m12,theta1),p) + pow(angleRatio(m30m12,m30v0,theta1),p))/8
               + w2*(pow(vectNorm(e01)/opt1Lavg,p) + pow(vectNorm(v1m12)/opt1Lavg,p)
               + pow(vectNorm(m12v2)/opt1Lavg,p) + pow(vectNorm(e23)/opt1Lavg,p)
               + pow(vectNorm(v3m30)/opt1Lavg,p) + pow(vectNorm(m30v0)/opt1Lavg,p)
               + pow(vectNorm(m12m30)/opt1Lavg,p))/7; 
    // option 2
    vector<double> m23v1 = vectSubtract(vertices[1],m23);
    vector<double> v1m23 = vectSubtract(m23,vertices[1]);
    double opt2Lavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(v2m23)+vectNorm(m23v3)+vectNorm(e30)+vectNorm(v1m23))/6;
    options[2] = w1*(pow(angleRatio(e03,e01,theta1),p) + pow(angleRatio(e10,v1m23,theta1),p)
               + pow(angleRatio(m23v1,m23v3,theta1),p) + pow(angleRatio(v3m23,e30,theta1),p)// quad
               + pow(angleRatio(v1m23,e12,theta0),p) + pow(angleRatio(e21,v2m23,theta0),p)
               + pow(angleRatio(m23v1,m23v2,theta0),p))/7
               + w2*(pow(vectNorm(e01)/opt2Lavg,p) + pow(vectNorm(e12)/opt2Lavg,p)
               + pow(vectNorm(v2m23)/opt2Lavg,p) + pow(vectNorm(m23v3)/opt2Lavg,p)
               + pow(vectNorm(e30)/opt2Lavg,p) + pow(vectNorm(v1m23)/opt2Lavg,p))/6; 
    // option 3
    vector<double> m30v2 = vectSubtract(vertices[2],m30);
    vector<double> v2m30 = vectSubtract(m30,vertices[2]);
    double opt3Lavg = (vectNorm(e12)+vectNorm(e23)+vectNorm(v3m30)+vectNorm(m30v0)+vectNorm(e01)+vectNorm(v2m30))/6;
    options[3] = w1*(pow(angleRatio(e10,e12,theta1),p) + pow(angleRatio(e21,v2m30,theta1),p)
               + pow(angleRatio(m30v2,m30v0,theta1),p) + pow(angleRatio(v0m30,e01,theta1),p)// quad
               + pow(angleRatio(v2m30,e23,theta0),p) + pow(angleRatio(e32,v3m30,theta0),p)
               + pow(angleRatio(m30v2,m30v3,theta0),p))/7
               + w2*(pow(vectNorm(e12)/opt3Lavg,p) + pow(vectNorm(e23)/opt3Lavg,p)
               + pow(vectNorm(v3m30)/opt3Lavg,p) + pow(vectNorm(m30v0)/opt3Lavg,p)
               + pow(vectNorm(e01)/opt3Lavg,p) + pow(vectNorm(v2m30)/opt3Lavg,p))/6; 
    // option 4
    vector<double> m01v3 = vectSubtract(vertices[3],m01);
    vector<double> v3m01 = vectSubtract(m01,vertices[3]);
    double opt4Lavg = (vectNorm(e23)+vectNorm(e30)+vectNorm(v0m01)+vectNorm(m01v1)+vectNorm(e12)+vectNorm(v3m01))/6;
    options[4] = w1*(pow(angleRatio(e21,e23,theta1),p) + pow(angleRatio(e32,v3m01,theta1),p)
               + pow(angleRatio(m01v3,m01v1,theta1),p) + pow(angleRatio(v1m01,e12,theta1),p)// quad
               + pow(angleRatio(v3m01,e30,theta0),p) + pow(angleRatio(e03,v0m01,theta0),p)
               + pow(angleRatio(m01v3,m01v0,theta0),p))/7
               + w2*(pow(vectNorm(e23)/opt4Lavg,p) + pow(vectNorm(e30)/opt4Lavg,p)
               + pow(vectNorm(v0m01)/opt4Lavg,p) + pow(vectNorm(m01v1)/opt4Lavg,p)
               + pow(vectNorm(e12)/opt4Lavg,p) + pow(vectNorm(v3m01)/opt4Lavg,p))/6; 
    // option 5
    vector<double> m12v0 = vectSubtract(vertices[0],m12);
    vector<double> v0m12 = vectSubtract(m12,vertices[0]);
    double opt5Lavg = (vectNorm(e30)+vectNorm(e01)+vectNorm(v1m12)+vectNorm(m12v2)+vectNorm(e23)+vectNorm(v0m12))/6;
    options[5] = w1*(pow(angleRatio(e32,e30,theta1),p) + pow(angleRatio(e03,v0m12,theta1),p)
               + pow(angleRatio(m12v0,m12v2,theta1),p) + pow(angleRatio(v2m12,e23,theta1),p)// quad
               + pow(angleRatio(v0m12,e01,theta0),p) + pow(angleRatio(e10,v1m12,theta0),p)
               + pow(angleRatio(m12v0,m12v1,theta0),p))/7
               + w2*(pow(vectNorm(e30)/opt5Lavg,p) + pow(vectNorm(e01)/opt5Lavg,p)
               + pow(vectNorm(v1m12)/opt5Lavg,p) + pow(vectNorm(m12v2)/opt5Lavg,p)
               + pow(vectNorm(e23)/opt5Lavg,p) + pow(vectNorm(v0m12)/opt5Lavg,p))/6; 
    // option 6
    vector<double> m23v0 = vectSubtract(vertices[0],m23);
    vector<double> v0m23 = vectSubtract(m23,vertices[0]);
    double opt6Lavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(v2m23)+vectNorm(m23v3)+vectNorm(e30)+vectNorm(v0m23))/6;
    options[6] = w1*(pow(angleRatio(e12,e10,theta1),p) + pow(angleRatio(e01,v0m23,theta1),p)
               + pow(angleRatio(m23v0,m23v2,theta1),p) + pow(angleRatio(v2m23,e21,theta1),p)// quad
               + pow(angleRatio(v0m23,e03,theta0),p) + pow(angleRatio(e30,v3m23,theta0),p)
               + pow(angleRatio(m23v0,m23v3,theta0),p))/7
               + w2*(pow(vectNorm(e01)/opt6Lavg,p) + pow(vectNorm(e12)/opt6Lavg,p)
               + pow(vectNorm(v2m23)/opt6Lavg,p) + pow(vectNorm(m23v3)/opt6Lavg,p)
               + pow(vectNorm(e30)/opt6Lavg,p) + pow(vectNorm(v0m23)/opt6Lavg,p))/6; 
    // option 7
    vector<double> m30v1 = vectSubtract(vertices[1],m30);
    vector<double> v1m30 = vectSubtract(m30,vertices[1]);
    double opt7Lavg = (vectNorm(e12)+vectNorm(e23)+vectNorm(v3m30)+vectNorm(m30v0)+vectNorm(e01)+vectNorm(v1m30))/6;
    options[7] = w1*(pow(angleRatio(e23,e21,theta1),p) + pow(angleRatio(e12,v1m30,theta1),p)
               + pow(angleRatio(m30v1,m30v3,theta1),p) + pow(angleRatio(v3m30,e32,theta1),p)// quad
               + pow(angleRatio(v1m30,e10,theta0),p) + pow(angleRatio(e01,v0m30,theta0),p)
               + pow(angleRatio(m30v1,m30v0,theta0),p))/7
               + w2*(pow(vectNorm(e12)/opt7Lavg,p) + pow(vectNorm(e23)/opt7Lavg,p)
               + pow(vectNorm(v3m30)/opt7Lavg,p) + pow(vectNorm(m30v0)/opt7Lavg,p)
               + pow(vectNorm(e01)/opt7Lavg,p) + pow(vectNorm(v1m30)/opt7Lavg,p))/6; 
    // option 8
    vector<double> m01v2 = vectSubtract(vertices[2],m01);
    vector<double> v2m01 = vectSubtract(m01,vertices[2]);
    double opt8Lavg = (vectNorm(e23)+vectNorm(e30)+vectNorm(v0m01)+vectNorm(m01v1)+vectNorm(e12)+vectNorm(v2m01))/6;
    options[8] = w1*(pow(angleRatio(e30,e32,theta1),p) + pow(angleRatio(e23,v2m01,theta1),p)
               + pow(angleRatio(m01v2,m01v0,theta1),p) + pow(angleRatio(v0m01,e03,theta1),p)// quad
               + pow(angleRatio(v2m01,e21,theta0),p) + pow(angleRatio(e12,v1m01,theta0),p)
               + pow(angleRatio(m01v2,m01v1,theta0),p))/7
               + w2*(pow(vectNorm(e23)/opt8Lavg,p) + pow(vectNorm(e30)/opt8Lavg,p)
               + pow(vectNorm(v0m01)/opt8Lavg,p) + pow(vectNorm(m01v1)/opt8Lavg,p)
               + pow(vectNorm(e12)/opt8Lavg,p) + pow(vectNorm(v2m01)/opt8Lavg,p))/6; 
    // option 9
    vector<double> m12v3 = vectSubtract(vertices[3],m12);
    vector<double> v3m12 = vectSubtract(m12,vertices[3]);
    double opt9Lavg = (vectNorm(e30)+vectNorm(e01)+vectNorm(v1m12)+vectNorm(m12v2)+vectNorm(e23)+vectNorm(v3m12))/6;
    options[9] = w1*(pow(angleRatio(e01,e03,theta1),p) + pow(angleRatio(e30,v3m12,theta1),p)
               + pow(angleRatio(m12v3,m12v1,theta1),p) + pow(angleRatio(v1m12,e10,theta1),p)// quad
               + pow(angleRatio(v3m12,e32,theta0),p) + pow(angleRatio(e23,v2m12,theta0),p)
               + pow(angleRatio(m12v3,m12v2,theta0),p))/7
               + w2*(pow(vectNorm(e30)/opt9Lavg,p) + pow(vectNorm(e01)/opt9Lavg,p)
               + pow(vectNorm(v1m12)/opt9Lavg,p) + pow(vectNorm(m12v2)/opt9Lavg,p)
               + pow(vectNorm(e23)/opt9Lavg,p) + pow(vectNorm(v3m12)/opt9Lavg,p))/6; 
    // option 10
    double opt10Lavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(e23)+vectNorm(e30))/4;
    options[10] = w1*(pow(angleRatio(e01,e03,theta1),p) + pow(angleRatio(e10,e12,theta1),p)
               + pow(angleRatio(e21,e23,theta1),p) + pow(angleRatio(e32,e30,theta1),p))/4
               + w2*(pow(vectNorm(e01)/opt10Lavg,p) + pow(vectNorm(e12)/opt10Lavg,p)
               + pow(vectNorm(e23)/opt10Lavg,p) + pow(vectNorm(e30)/opt10Lavg,p))/4
               + w3*(abs(dotProduct(cross,e01n)) + abs(dotProduct(cross,e12n)) + abs(dotProduct(cross,e23n)) + abs(dotProduct(cross,e30n)))/4;
    return options;
}

// create the lower and upper bounds for splitting optimization
vector<vector<double> > opt::buildbounds(vector<int> fixedvert)
{
    // fixedvert is an array of points, if an edge from the mesh has the same indices as two neighouring points from teh fixedvert, then the edge is fixed
    // for instance, edge (10,100) is fixed if the fixedvert is (0,2,3,4,5,100,10,324,4324) sice 10 and 100 are neighbouring points in fixedvert
    vector<vector<double> > bounds;
    vector<double> lowerbound;
    double arrlower[] = {0,0,0,0,0,0,0,0,0,0,0};
    vector<double> lower(arrlower,arrlower+11);
    vector<double> upperbound;
    for (unsigned int i=0;i<triConnectivity.size();i++){
        int v0 = triConnectivity[i][0], v1 = triConnectivity[i][1], v2 = triConnectivity[i][2];
        int f0 = 0, f1 = 0, f2 = 0; 
        if (fixedvert.size()>0){
            for (unsigned int j=0;j<fixedvert.size()-1;j++){
                if ((v0 == fixedvert[j] && v1 == fixedvert[j+1]) || (v1 == fixedvert[j] && v0 == fixedvert[j+1])){// v0 and v1 get fixed meaning the edge 0 is fixed
                    f0 = 1;
                }else if ((v1 == fixedvert[j] && v2 == fixedvert[j+1]) || (v2 == fixedvert[j] && v1 == fixedvert[j+1])){
                    f1 = 1;
                }else if ((v2 == fixedvert[j] && v0 == fixedvert[j+1]) || (v0 == fixedvert[j] && v2 == fixedvert[j+1])){
                    f2 = 1;
                }
            }
        }
        if (f2==1 && f0==0 && f1==0){
            double upperlist[] = {1,0,1,1,0,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// 1, 4, 5 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==0 && f2==0){
            double upperlist[] = {1,1,0,0,1,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// 2, 3, 5 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f1==1 && f2==0 && f0==0){
            double upperlist[] = {0,1,1,0,0,1,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// 0, 3, 4 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==0){
            double upperlist[] = {0,1,0,0,0,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// only 1 and 6 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f1==1 && f2==1 && f0==0){
            double upperlist[] = {0,0,1,0,0,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// only 2 and 6 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f2==1 && f1==0){
            double upperlist[] = {1,0,0,0,0,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// only 0 and 6 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==1){// should almost never happen
            double upperlist[] = {0,0,0,0,0,0,0,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// all not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else{
            double upperlist[] = {1,1,1,1,1,1,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// all not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }

        lowerbound.insert(lowerbound.end(),lower.begin(),lower.end());
    }
    for (unsigned int i=0;i<quadConnectivity.size();i++){
        int v0 = quadConnectivity[i][0], v1 = quadConnectivity[i][1], v2 = quadConnectivity[i][2], v3 = quadConnectivity[i][3];
        int f0 = 0, f1 = 0, f2 = 0, f3 = 0; 
        if (fixedvert.size()>0){
            for (unsigned int j=0;j<fixedvert.size()-1;j++){
                if ((v0 == fixedvert[j] && v1 == fixedvert[j+1]) || (v1 == fixedvert[j] && v0 == fixedvert[j+1])){
                    f0 = 1;
                }else if ((v1 == fixedvert[j] && v2 == fixedvert[j+1]) || (v2 == fixedvert[j] && v1 == fixedvert[j+1])){
                    f1 = 1;
                }else if ((v2 == fixedvert[j] && v3 == fixedvert[j+1]) || (v3 == fixedvert[j] && v2 == fixedvert[j+1])){
                    f2 = 1;
                }else if ((v3 == fixedvert[j] && v0 == fixedvert[j+1]) || (v0 == fixedvert[j] && v3 == fixedvert[j+1])){
                    f3 = 1;
                }
            }
        }
        if (f0==1 && f1==0 && f2==0 && f3==0){
            double upperlist[] = {0,1,1,1,0,1,1,1,0,1,1}; vector<double> upper(upperlist,upperlist+11);// 0, 4, 8 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f1==1 && f2==0 && f0==0 && f3==0){
            double upperlist[] = {1,0,1,1,1,0,1,1,1,0,1}; vector<double> upper(upperlist,upperlist+11);// 1, 5, 9 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f2==1 && f3==0 && f0==0 && f1==0){
            double upperlist[] = {0,1,0,1,1,1,0,1,1,1,1}; vector<double> upper(upperlist,upperlist+11);// 0, 2, 6 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f3==1 && f0==0 && f1==0 && f2==0){
            double upperlist[] = {1,0,1,0,1,1,1,0,1,1,1}; vector<double> upper(upperlist,upperlist+11);// 1, 3, 7 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==0 && f1==1 && f2==1 && f3==0){
            double upperlist[] = {0,0,0,1,1,0,0,1,1,0,1}; vector<double> upper(upperlist,upperlist+11);// 0, 1, 2, 5, 6, 9 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==0 && f1==0 && f2==1 && f3==1){
            double upperlist[] = {0,0,0,0,1,1,0,0,1,1,1}; vector<double> upper(upperlist,upperlist+11);// 0, 1, 2, 3, 6, 7 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==0 && f2==0 && f3==1){
            double upperlist[] = {0,0,1,0,0,1,1,0,0,1,1}; vector<double> upper(upperlist,upperlist+11);// 0, 1, 3, 4, 7, 8 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==0 && f3==0){
            double upperlist[] = {0,0,1,1,0,0,1,1,0,0,1}; vector<double> upper(upperlist,upperlist+11);// 0, 1, 4, 5, 8, 9 not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==1 && f3==0){
            double upperlist[] = {0,0,0,1,0,0,0,1,0,0,0}; vector<double> upper(upperlist,upperlist+11);// 3 and 7 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==0 && f1==1 && f2==1 && f3==1){
            double upperlist[] = {0,0,0,0,1,0,0,0,1,0,0}; vector<double> upper(upperlist,upperlist+11);// 4 and 8 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==0 && f2==1 && f3==1){
            double upperlist[] = {0,0,0,0,0,1,0,0,0,1,0}; vector<double> upper(upperlist,upperlist+11);// 5 and 9 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==0 && f3==1){
            double upperlist[] = {0,0,1,0,0,0,1,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// 2 and 6 available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else if (f0==1 && f1==1 && f2==1 && f3==1){
            double upperlist[] = {0,0,0,0,0,0,0,0,0,0,0}; vector<double> upper(upperlist,upperlist+11);// all not available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }else{
            double upperlist[] = {1,1,1,1,1,1,1,1,1,1,1}; vector<double> upper(upperlist,upperlist+11);// all available
            upperbound.insert(upperbound.end(),upper.begin(),upper.end());
        }
        lowerbound.insert(lowerbound.end(),lower.begin(),lower.end());
    }
    bounds.push_back(lowerbound);
    bounds.push_back(upperbound);
    return bounds;
}


// count unique edges
vector<vector<int> > opt::countedges()
{   
    unsigned int i, j;
    vector<int> v1, v2, v3 ,v4;
    vector<vector<int> > edges;
    for (i=0;i<triConnectivity.size();i++){
        v1.clear();v2.clear();v3.clear();
        v1.push_back(triConnectivity[i][0]);
        v1.push_back(triConnectivity[i][1]);
        v2.push_back(triConnectivity[i][1]);
        v2.push_back(triConnectivity[i][2]);
        v3.push_back(triConnectivity[i][2]);
        v3.push_back(triConnectivity[i][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        edges.push_back(v1);
        edges.push_back(v2);
        edges.push_back(v3);
    }
    for (j=0;j<quadConnectivity.size();j++){
        v1.clear();v2.clear();v3.clear();v4.clear();
        v1.push_back(quadConnectivity[j][0]);v1.push_back(quadConnectivity[j][1]);
        v2.push_back(quadConnectivity[j][1]);v2.push_back(quadConnectivity[j][2]);
        v3.push_back(quadConnectivity[j][2]);v3.push_back(quadConnectivity[j][3]);
        v4.push_back(quadConnectivity[j][3]);v4.push_back(quadConnectivity[j][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        sort(v4.begin(), v4.end());
        edges.push_back(v1);
        edges.push_back(v2);
        edges.push_back(v3);
        edges.push_back(v4);
    }
    sort(edges.begin(),edges.end());
    edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return edges;
}

// consider quadlist
vector<vector<int> > opt::countedges(vector<vector<int> > quadlist)
{   
    unsigned int i, j;
    vector<int> v1, v2, v3 ,v4;
    vector<vector<int> > edges;
    for (j=0;j<quadlist.size();j++){
        v1.clear();v2.clear();v3.clear();v4.clear();
        v1.push_back(quadlist[j][0]);v1.push_back(quadlist[j][1]);
        v2.push_back(quadlist[j][1]);v2.push_back(quadlist[j][2]);
        v3.push_back(quadlist[j][2]);v3.push_back(quadlist[j][3]);
        v4.push_back(quadlist[j][3]);v4.push_back(quadlist[j][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        sort(v4.begin(), v4.end());
        edges.push_back(v1);
        edges.push_back(v2);
        edges.push_back(v3);
        edges.push_back(v4);
    }
    sort(edges.begin(),edges.end());
    edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return edges;
}

// equality constraint for splitting optimization (the way I do it may be a little dumb)
EqualityConstraint opt::equcon()
{
    unsigned int i, j, k;
    vector<int> edgeindex, polyindex;// index
    vector<vector<int> > edges = countedges();
    // data, row and col are used to define the sparse matrix, b is the RHS of the linear constraint
    vector<int> data;
    vector<int> row;
    vector<int> col;
    vector<int> b;
    int sharededges =0;
    // check every edge if it is a shared edge
    for (i=0;i<edges.size();i++){
        // cout<<edges[i][0]<<" "<<edges[i][1]<<"\n";
        int e1 = edges[i][0];
        int e2 = edges[i][1];
        int count = 0;
        edgeindex.clear();
        polyindex.clear();
        // check if the edge belongs to a triangle
        for (j=0;j<triConnectivity.size();j++){
            if ((triConnectivity[j][0]==e1 && triConnectivity[j][1]==e2) || (triConnectivity[j][1]==e1 && triConnectivity[j][0]==e2)){
                count+=1;
                edgeindex.push_back(0);
                polyindex.push_back(j);
            }else if ((triConnectivity[j][1]==e1 && triConnectivity[j][2]==e2) || (triConnectivity[j][2]==e1 && triConnectivity[j][1]==e2)){
                count+=1;
                edgeindex.push_back(1);
                polyindex.push_back(j);
            }else if ((triConnectivity[j][2]==e1 && triConnectivity[j][0]==e2) || (triConnectivity[j][0]==e1 && triConnectivity[j][2]==e2)){
                count+=1;
                edgeindex.push_back(2);
                polyindex.push_back(j);
            }
        }
        // check if the edge belongs to a quad
        for (j=0;j<(quadConnectivity.size());j++){
            if ((quadConnectivity[j][0]==e1 && quadConnectivity[j][1]==e2) || (quadConnectivity[j][1]==e1 && quadConnectivity[j][0]==e2)){
                count+=1;
                edgeindex.push_back(0);
                polyindex.push_back(j+triConnectivity.size());
            }else if ((quadConnectivity[j][1]==e1 && quadConnectivity[j][2]==e2) || (quadConnectivity[j][2]==e1 && quadConnectivity[j][1]==e2)){
                count+=1;
                edgeindex.push_back(1);
                polyindex.push_back(j+triConnectivity.size());
            }else if ((quadConnectivity[j][2]==e1 && quadConnectivity[j][3]==e2) || (quadConnectivity[j][3]==e1 && quadConnectivity[j][2]==e2)){
                count+=1;
                edgeindex.push_back(2);
                polyindex.push_back(j+triConnectivity.size());
            }else if ((quadConnectivity[j][3]==e1 && quadConnectivity[j][0]==e2) || (quadConnectivity[j][0]==e1 && quadConnectivity[j][3]==e2)){
                count+=1;
                edgeindex.push_back(3);
                polyindex.push_back(j+triConnectivity.size());
            }
        }
        // check if the number of polygons that share the edge is 2, if it is, get the data, row, and column index
        if (count==2){
            vector<int> datai(3,1); datai.insert(datai.end(),3,-1); datai.insert(datai.end(),8,1); datai.insert(datai.end(),8,-1);datai.insert(datai.end(),22,1);
            vector<int> rowi(6,sharededges*3);rowi.insert(rowi.end(),16,sharededges*3+1);rowi.insert(rowi.end(),22,sharededges*3+2);
            vector<int> coli;
            // cout<<polyindex[1]<<"\n";
            int arrbi[] = {0,0,2};
            vector<int> bi(arrbi,arrbi+3);
            vector<int> temp1;
            vector<int> temp2;
            int arrtemp31[] = {0,1,2,3,4,5,6,7,8,9,10};
            arrayadd(arrtemp31,11,polyindex[0]*11);
            vector<int> temp31(arrtemp31,arrtemp31+11);
            int arrtemp32[] = {0,1,2,3,4,5,6,7,8,9,10};
            arrayadd(arrtemp32,11,polyindex[1]*11);
            // cout<<"arrtemp[0]: "<<arrtemp32[0]<<"\n";
            vector<int> temp32(arrtemp32,arrtemp32+11);
            for (k=0;k<2;k++){
                if (polyindex[k]<(int) triConnectivity.size()){
                    if (edgeindex[k]==0){
                        temp1.push_back(2+polyindex[k]*11);temp1.push_back(3+polyindex[k]*11);temp1.push_back(5+polyindex[k]*11);
                        temp2.push_back(0+polyindex[k]*11);temp2.push_back(1+polyindex[k]*11);temp2.push_back(4+polyindex[k]*11);temp2.push_back(6+polyindex[k]*11);
                        temp2.push_back(7+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }else if (edgeindex[k]==1){
                        temp1.push_back(0+polyindex[k]*11);temp1.push_back(3+polyindex[k]*11);temp1.push_back(4+polyindex[k]*11);
                        temp2.push_back(1+polyindex[k]*11);temp2.push_back(2+polyindex[k]*11);temp2.push_back(5+polyindex[k]*11);temp2.push_back(6+polyindex[k]*11);
                        temp2.push_back(7+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }else{
                        temp1.push_back(1+polyindex[k]*11);temp1.push_back(4+polyindex[k]*11);temp1.push_back(5+polyindex[k]*11);
                        temp2.push_back(0+polyindex[k]*11);temp2.push_back(2+polyindex[k]*11);temp2.push_back(3+polyindex[k]*11);temp2.push_back(6+polyindex[k]*11);
                        temp2.push_back(7+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }
                }else{
                    if (edgeindex[k]==0){
                        temp1.push_back(0+polyindex[k]*11);temp1.push_back(4+polyindex[k]*11);temp1.push_back(8+polyindex[k]*11);
                        temp2.push_back(1+polyindex[k]*11);temp2.push_back(2+polyindex[k]*11);temp2.push_back(3+polyindex[k]*11);temp2.push_back(5+polyindex[k]*11);
                        temp2.push_back(6+polyindex[k]*11);temp2.push_back(7+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }else if (edgeindex[k]==1){
                        temp1.push_back(1+polyindex[k]*11);temp1.push_back(5+polyindex[k]*11);temp1.push_back(9+polyindex[k]*11);
                        temp2.push_back(0+polyindex[k]*11);temp2.push_back(2+polyindex[k]*11);temp2.push_back(3+polyindex[k]*11);temp2.push_back(4+polyindex[k]*11);
                        temp2.push_back(6+polyindex[k]*11);temp2.push_back(7+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }else if (edgeindex[k]==2){
                        temp1.push_back(0+polyindex[k]*11);temp1.push_back(2+polyindex[k]*11);temp1.push_back(6+polyindex[k]*11);
                        temp2.push_back(1+polyindex[k]*11);temp2.push_back(3+polyindex[k]*11);temp2.push_back(4+polyindex[k]*11);temp2.push_back(5+polyindex[k]*11);
                        temp2.push_back(7+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }else{
                        temp1.push_back(1+polyindex[k]*11);temp1.push_back(3+polyindex[k]*11);temp1.push_back(7+polyindex[k]*11);
                        temp2.push_back(0+polyindex[k]*11);temp2.push_back(2+polyindex[k]*11);temp2.push_back(4+polyindex[k]*11);temp2.push_back(5+polyindex[k]*11);
                        temp2.push_back(6+polyindex[k]*11);temp2.push_back(8+polyindex[k]*11);temp2.push_back(9+polyindex[k]*11);temp2.push_back(10+polyindex[k]*11);
                    }
                }
            }
            coli.insert(coli.end(),temp1.begin(),temp1.end());coli.insert(coli.end(),temp2.begin(),temp2.end());
            coli.insert(coli.end(),temp31.begin(),temp31.end());coli.insert(coli.end(),temp32.begin(),temp32.end());
            data.insert(data.end(),datai.begin(),datai.end());
            row.insert(row.end(),rowi.begin(),rowi.end());
            col.insert(col.end(),coli.begin(),coli.end());
            b.insert(b.end(),bi.begin(),bi.end());
            sharededges+=1;
        }else if (count>2){
            cout<<"irregular shared edge\n";
        }
    }
    EqualityConstraint equcon(data,row,col,b,sharededges);
    return equcon;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------
// merging algorithm
// as the function name says
Indices opt::counttrisharededges()
{
    unsigned int i, j;
    vector<int> index;// index
    vector<vector<int> > triindex;// indices of two tris sharing the same edge
    vector<int> vindex;
    vector<vector<int> > verindex;
    vector<int> edgeindex;// shared edge
    vector<vector<int> > edges = countedges();
    for (i=0;i<edges.size();i++){
        int e1 = edges[i][0];
        int e2 = edges[i][1];
        int count = 0;
        index.clear();
        vindex.clear();
        // check if the edge belongs to a triangle
        for (j=0;j<triConnectivity.size();j++){
            if ((triConnectivity[j][0]==e1 && triConnectivity[j][1]==e2) || (triConnectivity[j][1]==e1 && triConnectivity[j][0]==e2)){
                count+=1;
                index.push_back(j);
                vindex.push_back(2);
            }else if ((triConnectivity[j][1]==e1 && triConnectivity[j][2]==e2) || (triConnectivity[j][2]==e1 && triConnectivity[j][1]==e2)){
                count+=1;
                index.push_back(j);
                vindex.push_back(0);
            }else if ((triConnectivity[j][2]==e1 && triConnectivity[j][0]==e2) || (triConnectivity[j][0]==e1 && triConnectivity[j][2]==e2)){
                count+=1;
                index.push_back(j);
                vindex.push_back(1);
            }
        }
        if (count==2){
            triindex.push_back(index);
            edgeindex.push_back(i);
            verindex.push_back(vindex);
        }
    }
    Indices indices(triindex, edgeindex, verindex);
    return indices;
}

// build up coefficient matrix for merging optimization
vector<double> opt::mergetriangles()
{
    unsigned int i, j;
    Indices indices = counttrisharededges();
    vector<vector<int> > edges = countedges();
    vector<int> univectice;
    vector<double> data(edges.size(),0.0);

    for (i=0;i<indices.triindex.size();i++){
        univectice.clear();
        double energybefore = 0;
        double energyafter = 0;
        double angle0=0, angle1=0;
        vector<double> angle;
        for (j=0;j<2;j++){
            vector<double> e01 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][1]],vertexCoords[triConnectivity[indices.triindex[i][j]][0]]);
            vector<double> e10 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][0]],vertexCoords[triConnectivity[indices.triindex[i][j]][1]]);
            vector<double> e12 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][2]],vertexCoords[triConnectivity[indices.triindex[i][j]][1]]);
            vector<double> e21 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][1]],vertexCoords[triConnectivity[indices.triindex[i][j]][2]]);
            vector<double> e20 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][0]],vertexCoords[triConnectivity[indices.triindex[i][j]][2]]);
            vector<double> e02 = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][2]],vertexCoords[triConnectivity[indices.triindex[i][j]][0]]);
            double triLavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(e20))/3;
            energybefore += 0.5*(w1*(pow(angleRatio(e10,e12,theta0),p) + pow(angleRatio(e21,e20,theta0),p) + pow(angleRatio(e02,e01,theta0),p))/3
                            + w2*(pow(vectNorm(e01)/triLavg,p) + pow(vectNorm(e12)/triLavg,p) + pow(vectNorm(e20)/triLavg,p))/3);
            univectice.insert(univectice.begin(),triConnectivity[indices.triindex[i][j]].begin(),triConnectivity[indices.triindex[i][j]].end());
            // sum of merged angles
            vector<double> s0v = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][indices.verindex[i][j]]],vertexCoords[edges[indices.edgeindex[i]][0]]);// vector 0v
            vector<double> vs0 = vectSubtract(vertexCoords[edges[indices.edgeindex[i]][0]],vertexCoords[triConnectivity[indices.triindex[i][j]][indices.verindex[i][j]]]);// vector v0
            vector<double> s1v = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][j]][indices.verindex[i][j]]],vertexCoords[edges[indices.edgeindex[i]][1]]);// vector 1v
            vector<double> vs1 = vectSubtract(vertexCoords[edges[indices.edgeindex[i]][1]],vertexCoords[triConnectivity[indices.triindex[i][j]][indices.verindex[i][j]]]);// vector v1
            vector<double> sf = vectSubtract(vertexCoords[edges[indices.edgeindex[i]][1]],vertexCoords[edges[indices.edgeindex[i]][0]]);// vector 01
            vector<double> sr = vectSubtract(vertexCoords[edges[indices.edgeindex[i]][0]],vertexCoords[edges[indices.edgeindex[i]][1]]);// vector 10
            angle0 += angleRatio(s0v,sf,theta1);
            angle1 += angleRatio(s1v,sr,theta1);
            angle.push_back(angleRatio(vs0,vs1,theta1));
        }
        vector<int>::iterator ip;
        sort(univectice.begin(), univectice.end()); // Sorting the array 
        unique(univectice.begin(), univectice.end());// Using std::unique 
        vector<double > e01= vectSubtract(vertexCoords[univectice[1]],vertexCoords[univectice[0]]);
        vector<double > e12= vectSubtract(vertexCoords[univectice[2]],vertexCoords[univectice[1]]);
        vector<double > e23= vectSubtract(vertexCoords[univectice[3]],vertexCoords[univectice[2]]);
        vector<double > e30= vectSubtract(vertexCoords[univectice[0]],vertexCoords[univectice[3]]);
        double quadLavg = (vectNorm(e01)+vectNorm(e12)+vectNorm(e23)+vectNorm(e30))/4;
        // third energy involving dihedral angle
        vector<double> edgevector = vectSubtract(vertexCoords[edges[indices.edgeindex[i]][1]],vertexCoords[edges[indices.edgeindex[i]][0]]);
        vector<double> diagvector = vectSubtract(vertexCoords[triConnectivity[indices.triindex[i][1]][indices.verindex[i][1]]],vertexCoords[triConnectivity[indices.triindex[i][0]][indices.verindex[i][0]]]);
        vector<double> cross = crossProduct(edgevector,diagvector);
        cross = normalize(cross);
        vector<double> e01n = normalize(e01), e12n = normalize(e12), e23n = normalize(e23), e30n = normalize(e30);
        // energy after the merging is calculated by considering both the regularity of the angles, edges and the dihedral angles
        energyafter = w1*(pow(angle0,p) + pow(angle1,p) + pow(angle[0],p) + pow(angle[1],p))/4
                            + w2*(pow(vectNorm(e01)/quadLavg,p) + pow(vectNorm(e12)/quadLavg,p) + pow(vectNorm(e23)/quadLavg,p) + pow(vectNorm(e30)/quadLavg,p))/4
                            + w3*(abs(dotProduct(cross,e01n)) + abs(dotProduct(cross,e12n)) + abs(dotProduct(cross,e23n)) + abs(dotProduct(cross,e30n)))/4;
        // the coefficent is calculated by subtracting energy before merging from the energy after merging
        double w4 = 0.7;
        double energycost = w4*energyafter - energybefore;//////////w4
        data[indices.edgeindex[i]] = energycost;

    }
    return data;
}

// find edges that belong to triangles
vector<vector<int> > opt::definetriedges()
{
    unsigned int i, j;
    vector<vector<int> > edges = countedges();
    vector<vector<int> > triedges(triConnectivity.size(), vector<int> (3,0.0));
    vector<int> v1, v2, v3;
    for (i=0;i<triConnectivity.size();i++){
        v1.clear();v2.clear();v3.clear();
        v1.push_back(triConnectivity[i][0]);
        v1.push_back(triConnectivity[i][1]);
        v2.push_back(triConnectivity[i][1]);
        v2.push_back(triConnectivity[i][2]);
        v3.push_back(triConnectivity[i][2]);
        v3.push_back(triConnectivity[i][0]);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        sort(v3.begin(), v3.end());
        for (j=0; j<edges.size();j++){
            if (v1==edges[j]){
                triedges[i][0]=j;
            }
            if (v2==edges[j]){
                triedges[i][1]=j;
            }
            if (v3==edges[j]){
                triedges[i][2]=j;
            }
        }
    }
    return triedges;
}

// create the inequality constraints
vector<vector<int> > opt::inequncon()
{
    unsigned int i, j;
    vector<vector<int> > edges = countedges();
    vector<vector<int> > A_uq(triConnectivity.size(), vector<int> (edges.size(),0.0));
    // A_uq.reserve(triConnectivity.size());
    vector<vector<int> > triedges = definetriedges();
    for (i=0;i<triedges.size();i++){
        vector<int> A_uqi(edges.size(),0.0);
        for (j=0;j<3;j++){
            A_uqi[triedges[i][j]] = 1.0;
        }
        A_uq[i]=A_uqi;
    }
    return A_uq;
}
// create the upper and lower bound of merging optimization
vector<vector<double> > opt::mergingbounds(vector<int> fixedvert) 
{
    // just set the upper bounds of the fixed edge to 0
    unsigned int i, j;
    vector<vector<int> > edges = countedges();
    vector<double> lowerbound(edges.size(),0.0);
    vector<double> upperbound(edges.size(),1.0);
    vector<vector<double> > bounds;
    Indices sharededges = counttrisharededges();
    for (i=0;i<edges.size();i++){
        if (count(sharededges.edgeindex.begin(),sharededges.edgeindex.end(),i)){
            double e0 = edges[i][0], e1 = edges[i][1];
            int f0 = 0, f1 = 0;
            if (fixedvert.size()>0){
                for (j=0;j<fixedvert.size()-1;j++){
                    if ((e0 == fixedvert[j] && e1 ==  fixedvert[j+1]) || (e1 == fixedvert[j] && e0 == fixedvert[j+1])){
                        f0 = 1;
                        f1 = 1;
                    }
                }
                if (f0==1 && f1==1){
                    upperbound[i] = 0.0;
                }
            }
        }else{
            upperbound[i] = 0.0;
        }
    }
    bounds.push_back(lowerbound);
    bounds.push_back(upperbound);
    return bounds;
}

// transformation to fully quad mesh
FullyQuad opt::fullyquadmesh(vector<int> fixed)
{   
    // splitting triangles
    vector<vector<double> > vertlist = vertexCoords;
    vector<vector<int> > quadlist;
    unsigned int i, j, k;
    for (i=0;i<triConnectivity.size();i++){
        vector<double> v0 = vertlist[triConnectivity[i][0]];
        vector<double> v1 = vertlist[triConnectivity[i][1]];
        vector<double> v2 = vertlist[triConnectivity[i][2]];
        vector<double> m01 = midPoint(v0,v1);
        vector<double> m12 = midPoint(v1,v2);
        vector<double> m20 = midPoint(v2,v0);
        vector<double> center = midPoint(v0, v1, v2);
        vertlist.push_back(center);
        int centeri = vertlist.size()-1;
        int m01i = 0, m12i =0, m20i = 0;
        for (j=0;j<vertlist.size();j++){
            if (distBetweenPoints(m01,vertlist[j])<1e-7){
                m01i = j;
            }
            if (distBetweenPoints(m12,vertlist[j])<1e-7){
                m12i = j;
            }
            if (distBetweenPoints(m20,vertlist[j])<1e-7){
                m20i = j;
            }
        }
        if (m01i == 0){
            vertlist.push_back(m01);
            m01i = vertlist.size()-1;
        }
        if (m12i == 0){
            vertlist.push_back(m12);
            m12i = vertlist.size()-1;
        }
        if (m20i == 0){
            vertlist.push_back(m20);
            m20i = vertlist.size()-1;
        }
        int quad1list[] = {triConnectivity[i][0],m01i,centeri,m20i}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int) );
        int quad2list[] = {triConnectivity[i][1],m12i,centeri,m01i}; vector<int> quad2 (quad2list, quad2list + sizeof(quad2list) / sizeof(int) );
        int quad3list[] = {triConnectivity[i][2],m20i,centeri,m12i}; vector<int> quad3 (quad3list, quad3list + sizeof(quad3list) / sizeof(int) );
        quadlist.push_back(quad1);quadlist.push_back(quad2);quadlist.push_back(quad3);
    }

    // splitting quads
    for (i=0;i<quadConnectivity.size();i++){
        vector<double> v0 = vertlist[quadConnectivity[i][0]];
        vector<double> v1 = vertlist[quadConnectivity[i][1]];
        vector<double> v2 = vertlist[quadConnectivity[i][2]];
        vector<double> v3 = vertlist[quadConnectivity[i][3]];
        vector<double> m01 = midPoint(v0,v1);
        vector<double> m12 = midPoint(v1,v2);
        vector<double> m23 = midPoint(v2,v3);
        vector<double> m30 = midPoint(v3,v0);
        vector<double> center = midPoint(v0, v1, v2, v3);
        vertlist.push_back(center);
        int centeri = vertlist.size()-1;
        int m01i = 0, m12i = 0, m23i = 0, m30i=0;
        for (j=0;j<vertlist.size();j++){
            if (distBetweenPoints(m01,vertlist[j])<1e-7){
                m01i = j;
            }
            if (distBetweenPoints(m12,vertlist[j])<1e-7){
                m12i = j;
            }
            if (distBetweenPoints(m23,vertlist[j])<1e-7){
                m23i = j;
            }
            if (distBetweenPoints(m30,vertlist[j])<1e-7){
                m30i = j;
            }
        }
        if (m01i == 0){
            vertlist.push_back(m01);
            m01i = vertlist.size()-1;
        }
        if (m12i == 0){
            vertlist.push_back(m12);
            m12i = vertlist.size()-1;
        }
        if (m23i == 0){
            vertlist.push_back(m23);
            m23i = vertlist.size()-1;
        }
        if (m30i == 0){
            vertlist.push_back(m30);
            m30i = vertlist.size()-1;
        }
        int quad1list[] = {quadConnectivity[i][0],m01i,centeri,m30i}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
        int quad2list[] = {quadConnectivity[i][1],m12i,centeri,m01i}; vector<int> quad2 (quad2list, quad2list + sizeof(quad2list) / sizeof(int));
        int quad3list[] = {quadConnectivity[i][2],m23i,centeri,m12i}; vector<int> quad3 (quad3list, quad3list + sizeof(quad3list) / sizeof(int));
        int quad4list[] = {quadConnectivity[i][3],m30i,centeri,m23i}; vector<int> quad4 (quad4list, quad4list + sizeof(quad4list) / sizeof(int));
        quadlist.push_back(quad1);quadlist.push_back(quad2);quadlist.push_back(quad3);quadlist.push_back(quad4);
    }
    vector<vector<int> > edges = countedges();
    vector<int> fixedvert = fixed;
    FullyQuad fullyquad(quadlist,vertlist,fixedvert);
                
    return fullyquad;
}

// Calculate the normal vector of each vertex and measure the curvature to determine the weight for smoothing optimization
// this function should be robust, the reason of failure may comes from invalid mesh
OptPrep opt::computenormal()
{
    unsigned int i, j, k, l;
    vector<vector<double> > vertlist= vertexCoords;
    vector<vector<int> > quadlist = quadConnectivity;
    vector<vector<int> > trilist = triConnectivity;
    vector<vector<double> > normal;
    vector<double> weight;
    for (i=0;i<vertlist.size();i++){
        vector<vector<double> > normalj;
        for (j=0;j<quadlist.size();j++){
            for (k=0;k<4;k++){
                vector<double> e1, e2, e3, e4;
                if (((int) i)==quadlist[j][k]){// see which point it is in thequad
                    if (k==0){
                        e1 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                        e2 = vectSubtract(vertlist[quadlist[j][0]],vertlist[quadlist[j][1]]);
                        e3 = vectSubtract(vertlist[quadlist[j][1]],vertlist[quadlist[j][2]]);
                        e4 = vectSubtract(vertlist[quadlist[j][2]],vertlist[quadlist[j][3]]);
                    }else if (k==3){
                        e1 = vectSubtract(vertlist[quadlist[j][2]],vertlist[quadlist[j][3]]);
                        e2 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                        e3 = vectSubtract(vertlist[quadlist[j][0]],vertlist[quadlist[j][1]]);
                        e4 = vectSubtract(vertlist[quadlist[j][1]],vertlist[quadlist[j][2]]);
                    }else{
                        e1 = vectSubtract(vertlist[quadlist[j][k-1]],vertlist[quadlist[j][k]]);
                        e2 = vectSubtract(vertlist[quadlist[j][k]],vertlist[quadlist[j][k+1]]);
                        if (k==1){
                            e3 = vectSubtract(vertlist[quadlist[j][2]],vertlist[quadlist[j][3]]);
                            e4 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                        }else {
                            e3 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                            e4 = vectSubtract(vertlist[quadlist[j][0]],vertlist[quadlist[j][1]]);
                        }
                    }
                    vector<double> normalk = crossProduct(e1,e2);
                    double area = 0.5*sqrt(dotProduct(normalk,normalk));
                    if (area==0){// area == 0 means the three points are colinear
                        normalk = crossProduct(e3,e4);
                        area = 0.5*sqrt(dotProduct(normalk,normalk));
                    }

                    if (isnan(normalk[0]) || isnan(normalk[1]) || isnan(normalk[2])){
                        cout<<"quad1: "<<j<<'\n';
                        cout<<e1[0]<<' '<<e1[1]<<' '<<e1[2]<<'\n';
                        cout<<e2[0]<<' '<<e2[1]<<' '<<e2[2]<<'\n';
                    }
                    normalk = normalize(normalk);
                    if (isnan(normalk[0]) || isnan(normalk[1]) || isnan(normalk[2])){
                        cout<<"quad2: "<<j<<'\n';
                    }
                    normalk = vectScale(normalk,area);
                    normalj.push_back(normalk);
                }
            }
        }
        for (j=0;j<trilist.size();j++){
            for (k=0;k<3;k++){
                vector<double> e1, e2;
                if (((int)i)==trilist[j][k]){
                    if (k==0){
                        e1 = vectSubtract(vertlist[trilist[j][2]],vertlist[trilist[j][0]]);
                        e2 = vectSubtract(vertlist[trilist[j][0]],vertlist[trilist[j][1]]);
                    }else if (k==1){
                        e1 = vectSubtract(vertlist[trilist[j][0]],vertlist[trilist[j][1]]);
                        e2 = vectSubtract(vertlist[trilist[j][1]],vertlist[trilist[j][2]]);
                    }else{
                        e1 = vectSubtract(vertlist[trilist[j][1]],vertlist[trilist[j][2]]);
                        e2 = vectSubtract(vertlist[trilist[j][2]],vertlist[trilist[j][0]]);
                    }
                    vector<double> normalk = crossProduct(e1,e2);
                    double area = 0.5*sqrt(dotProduct(normalk,normalk));
                    // below are a few checks to see the validity
                    if (isnan(normalk[0]) || isnan(normalk[1]) || isnan(normalk[2])){
                        cout<<"tri1: "<<j<<'\n';
                        cout<<e1[0]<<' '<<e1[1]<<' '<<e1[2]<<'\n';
                        cout<<e2[0]<<' '<<e2[1]<<' '<<e2[2]<<'\n';
                        cout<<normalk[0]<<' '<<normalk[1]<<' '<<normalk[2]<<'\n';
                    }
                    if (area==0){
                        cout<<"no\n";
                        cout<<"normalk: "<<normalk[0]<<' '<<normalk[1]<<' '<<normalk[2]<<'\n';
                        cout<<"e1: "<<e1[0]<<' '<<e1[1]<<' '<<e1[2]<<'\n';
                        cout<<"e2: "<<e2[0]<<' '<<e2[1]<<' '<<e2[2]<<'\n';
                    }
                    normalk = normalize(normalk);
                    if (isnan(normalk[0]) || isnan(normalk[1]) || isnan(normalk[2])){
                        cout<<"tri2: "<<j<<'\n';
                    }
                    normalk = vectScale(normalk,area);
                    normalj.push_back(normalk);
                }
            }
        }
        vector<double> normali = vectMean(normalj);
        normali = normalize(normali);
        // too be honest, I forgot why I did the following, but it is reasonable
        if ((isnan(normali[0]) || isnan(normali[1]) || isnan(normali[2])) && normalj.size()==2){
            vector<double> e1, e2;
            for (j=0;j<quadlist.size();j++){
                for (k=0;k<4;k++){
                    if (((int) i)==quadlist[j][k]){
                        if (k==0){
                            e1 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                            e2 = vectSubtract(vertlist[quadlist[j][0]],vertlist[quadlist[j][1]]);
                        }else if (k==3){
                            e1 = vectSubtract(vertlist[quadlist[j][2]],vertlist[quadlist[j][3]]);
                            e2 = vectSubtract(vertlist[quadlist[j][3]],vertlist[quadlist[j][0]]);
                        }else{
                            e1 = vectSubtract(vertlist[quadlist[j][k-1]],vertlist[quadlist[j][k]]);
                            e2 = vectSubtract(vertlist[quadlist[j][k]],vertlist[quadlist[j][k+1]]);
                        }
                    }
                }
            }
            normali = vectSubtract(e1,e2);
        }
        if (isnan(normali[0]) || isnan(normali[1]) || isnan(normali[2])){
            for (j=0;j<trilist.size();j++){
                for (k=0;k<3;k++){
                    if (((int) i)==trilist[j][k]){
                        cout<<i<<' '<<j<<' '<<trilist[j][0]<<' '<<trilist[j][1]<<' '<<trilist[j][2]<<'\n';
                    }
                }
            }
        }
        double weighti = 1;
        for (l=0;l<normalj.size();l++){
            vector<double> normaljl = normalj[l]; normaljl = normalize(normaljl);
            double weightl = dotProduct(normali,normaljl);
            if (weightl < weighti){
                weighti = weightl;
            }
        }
        normal.push_back(normali);
        weight.push_back(weighti);
    }
    OptPrep optprep(normal,weight);
    return optprep;
}

double opt::computeweight(vector<vector<double> > directions)
{   
    // weight is basically the smallest dot product between the two neighbouring normal vectors
    unsigned int i, j;
    double weight = -1e5;
    for (i=0;i<directions.size();i++){
        for (j=i;j<directions.size();j++){
            vector<double> direct1 = directions[i];
            vector<double> direct2 = directions[j];
            double weighti = dotProduct(direct1,direct2);
            if (weighti>weight){
                weight = weighti;
            }
        }
    }
    return weight;
}

// create a coordinate system for each vertex on its tangential plane used for smoothing optimization
vector<vector<vector<double> > > opt::computecs()
{
    unsigned int i;
    OptPrep optprep = computenormal();
    vector<vector<double> > normal = optprep.normal;
    vector<vector<double> > u(normal.size(),vector<double>(3));
    vector<vector<double> > v(normal.size(),vector<double>(3));
    for (i=0;i<normal.size();i++){
        vector<double> normali = normal[i];
        double randlist[] = {1.242345,2.242345,3.242345}; vector<double> rand (randlist, randlist + sizeof(randlist) / sizeof(double));
        vector<double> ui = crossProduct(normali,rand);
        vector<double> vi = crossProduct(normali,ui);
        ui = normalize(ui);vi=normalize(vi);
        u[i] = ui;
        v[i] = vi;
    }
    vector<vector<vector<double> > > cs;
    cs.push_back(u);cs.push_back(v);
    return cs;
}

// used in qcp and no longer necessary
vector<double> opt::computeradius()
{
    unsigned int i, j;
    OptPrep optprep = computenormal();
    vector<vector<double> > normal = optprep.normal; 
    vector<vector<int> > edges = countedges();
    vector<double> radius;
    double r, d, p, t;
    for (i=0;i<vertexCoords.size();i++){
        vector<double> n = normal[i];
        r = 1e10;
        for (j=0;j<edges.size();j++){
            if (edges[j][0]==(int)i || edges[j][1]==(int)i){
                vector<double> e0 = vertexCoords[edges[j][0]];
                vector<double> e1 = vertexCoords[edges[j][1]];
                vector<double> e = vectSubtract(e0,e1);
                d = distBetweenPoints(vertexCoords[edges[j][0]], vertexCoords[edges[j][1]]);
                p = abs(dotProduct(n,e));
                t = d*d-p*p;
                if (t<r){
                    r = t;
                }
            }
        }
        radius.push_back(r);
    }  
    return radius;
}

// based on the splitting optimizaiton results, update the mesh
SplitOutput opt::splitupdatemesh(vector<int> optres, vector<int> fixed)
{
    vector<vector<int> > trilist = triConnectivity;
    vector<vector<int> > quadlist = quadConnectivity;
    vector<int> fixedvert;
    unsigned int i, j, k;
    int trilistlen = trilist.size();
    int quadlistlen = quadlist.size();
    vector<int> trioptres(optres.begin(),optres.begin()+trilistlen);
    vector<int> quadoptres(optres.begin()+trilistlen,optres.begin()+trilistlen+quadlistlen);
    UpdateMesh trimesh = updatetrisplit(trioptres,trilist,vertexCoords);
    UpdateMesh quadmesh = updatequadsplit(quadoptres,quadlist,trimesh.vertlist);
    vector<int> tridelete = trimesh.tridelete;
    tridelete.insert(tridelete.end(),quadmesh.tridelete.begin(),quadmesh.tridelete.end());
    sort(tridelete.begin(),tridelete.end(),greater<int>());
    for (i=0;i<tridelete.size();i++){
        trilist.erase(trilist.begin()+tridelete[i]);
    }
    for (i=0;i<trimesh.triadd.size();i++){
        trilist.push_back(trimesh.triadd[i]);
    }
    for (i=0;i<quadmesh.triadd.size();i++){
        trilist.push_back(quadmesh.triadd[i]);
    }
    // delete and add quads
    vector<int> quaddelete = quadmesh.quaddelete;
    sort(quaddelete.begin(),quaddelete.end(),greater<int>());
    for (i=0;i<quaddelete.size();i++){
        quadlist.erase(quadlist.begin()+quaddelete[i]);
    }
    for (i=0;i<trimesh.quadadd.size();i++){
        quadlist.push_back(trimesh.quadadd[i]);
    }
    for (i=0;i<quadmesh.quadadd.size();i++){
        quadlist.push_back(quadmesh.quadadd[i]);
    }

    // get unique vertices
    vector<vector<double> > vertlistr = quadmesh.vertlist;
    vector<vector<double> > vertlistu = quadmesh.vertlist;
    for (i=0;i<vertlistu.size();i++){
        vector<double> vi = vertlistu[i];
        for (j=i+1;j<vertlistu.size();j++){
            vector<double> vj = vertlistu[j];
            if (distBetweenPoints(vi,vj)<1e-10){
                vertlistu.erase(vertlistu.begin()+j);
            }
        }
    }
    for (i=0;i<fixed.size();i++){
        vector<double> v0 = vertexCoords[fixed[i]];
        for (j=0;j<vertlistu.size();j++){
            vector<double> v1 = vertlistu[j];
            // cout<<"distance"<<distBetweenPoints(v0,v1)<<'\n';
            if (distBetweenPoints(v0,v1)<1e-7){
                // cout<<"nice!\n";
                fixedvert.push_back(j);
            }
        }
    }
    for (i=0;i<trilist.size();i++){
        for (j=0;j<3;j++){
            for (k=0;k<vertlistu.size();k++){
                if (distBetweenPoints(vertlistr[trilist[i][j]],vertlistu[k])<1e-7){
                    trilist[i][j]=k;
                    break;
                }
            }
        }
    }
    for (i=0;i<quadlist.size();i++){
        for (j=0;j<4;j++){
            for (k=0;k<vertlistu.size();k++){
                if (distBetweenPoints(vertlistr[quadlist[i][j]],vertlistu[k])<1e-7){
                    quadlist[i][j]=k;
                    break;
                }
            }
        }
    }
    // get splitcoord
    vector<vector<double> > splitcoord = trimesh.splitcoord;
    splitcoord.insert(splitcoord.end(),quadmesh.splitcoord.begin(),quadmesh.splitcoord.end());
    SplitOutput output(vertlistu,splitcoord,trilist,quadlist,fixedvert);
    return output;
}

// split triangles which belongs to spliting optimization
UpdateMesh opt::updatetrisplit(vector<int> optres, vector<vector<int> > trilist, vector<vector<double> > vertlist)
{
    unsigned int i;
    vector<vector<int> > triadd;
    vector<vector<int> > quadadd; 
    vector<int> tridelete;
    vector<int> quaddelete;
    vector<vector<double> > splitcoord;
    for (i=0;i<trilist.size();i++){
        int l;
        if (optres[i]==0){
            vector<double> vert = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][0],trilist[i][1],l}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int tri2list[] = {trilist[i][0],l,trilist[i][2]}; vector<int> tri2 (tri2list, tri2list + sizeof(tri2list) / sizeof(int));
            triadd.push_back(tri1);triadd.push_back(tri2);
            vector<double> split1 = vertlist[trilist[i][0]];
            vector<double> split2 = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==1){
            vector<double> vert = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][0],trilist[i][1],l}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int tri2list[] = {trilist[i][1],trilist[i][2],l}; vector<int> tri2 (tri2list, tri2list + sizeof(tri2list) / sizeof(int));
            triadd.push_back(tri1);triadd.push_back(tri2);
            vector<double> split1 = vertlist[trilist[i][1]];
            vector<double> split2 = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==2){
            vector<double> vert = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][0],l,trilist[i][2]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int tri2list[] = {trilist[i][1],trilist[i][2],l}; vector<int> tri2 (tri2list, tri2list + sizeof(tri2list) / sizeof(int));
            triadd.push_back(tri1);triadd.push_back(tri2);
            vector<double> split1 = vertlist[trilist[i][2]];
            vector<double> split2 = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==3){
            vector<double> vert1 = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            vector<double> vert2 = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            vertlist.push_back(vert1);vertlist.push_back(vert2);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][1],l,l-1}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {trilist[i][0],l-1,l,trilist[i][2]}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            vector<double> split2 = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==4){
            vector<double> vert1 = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            vector<double> vert2 = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            vertlist.push_back(vert1);vertlist.push_back(vert2);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][2],l,l-1}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {trilist[i][1],l-1,l,trilist[i][0]}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = midPoint(vertlist[trilist[i][1]],vertlist[trilist[i][2]]);
            vector<double> split2 = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==5){
            vector<double> vert1 = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            vector<double> vert2 = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            vertlist.push_back(vert1);vertlist.push_back(vert2);
            l = vertlist.size()-1;
            tridelete.push_back(i);
            int tri1list[] = {trilist[i][0],l,l-1}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {trilist[i][2],l-1,l,trilist[i][1]}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = midPoint(vertlist[trilist[i][2]],vertlist[trilist[i][0]]);
            vector<double> split2 = midPoint(vertlist[trilist[i][0]],vertlist[trilist[i][1]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }
    }
    UpdateMesh updatemesh(splitcoord,triadd,quadadd,tridelete,quaddelete,vertlist);
    return updatemesh;
}

// split quads which belongs to spliting optimization
UpdateMesh opt::updatequadsplit(vector<int> optres, vector<vector<int> > quadlist, vector<vector<double> > vertlist)
{
    unsigned int i;
    vector<vector<int> > triadd;
    vector<vector<int> > quadadd; 
    vector<int> tridelete;
    vector<int> quaddelete;
    vector<vector<double> > splitcoord;
    for (i=0;i<quadlist.size();i++){
        int l;
        if (optres[i]==0){
            vector<double> vert1 = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            vector<double> vert2 = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            vertlist.push_back(vert1);vertlist.push_back(vert2);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int quad1list[] = {quadlist[i][2],l,l-1,quadlist[i][1]}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            int quad2list[] = {quadlist[i][0],l-1,l,quadlist[i][3]}; vector<int> quad2 (quad2list, quad2list + sizeof(quad2list) / sizeof(int));
            quadadd.push_back(quad1);quadadd.push_back(quad2);
            vector<double> split1 = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            vector<double> split2 = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==1){
            vector<double> vert1 = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            vector<double> vert2 = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            vertlist.push_back(vert1);vertlist.push_back(vert2);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int quad1list[] = {quadlist[i][3],l,l-1,quadlist[i][2]}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            int quad2list[] = {quadlist[i][1],l-1,l,quadlist[i][0]}; vector<int> quad2 (quad2list, quad2list + sizeof(quad2list) / sizeof(int));
            quadadd.push_back(quad1);quadadd.push_back(quad2);
            vector<double> split1 = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            vector<double> split2 = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==2){
            vector<double> vert = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][2],l,quadlist[i][1]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][3],quadlist[i][0],quadlist[i][1],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][1]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==3){
            vector<double> vert = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][3],l,quadlist[i][2]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][0],quadlist[i][1],quadlist[i][2],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][2]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==4){
            vector<double> vert = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][0],l,quadlist[i][3]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][1],quadlist[i][2],quadlist[i][3],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][3]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==5){
            vector<double> vert = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][1],l,quadlist[i][0]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][2],quadlist[i][3],quadlist[i][0],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][0]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==6){
            vector<double> vert = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][0],l,quadlist[i][3]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][0],quadlist[i][1],quadlist[i][2],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][0]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][2]],vertlist[quadlist[i][3]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==7){
            vector<double> vert = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][1],l,quadlist[i][0]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][1],quadlist[i][2],quadlist[i][3],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][1]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][3]],vertlist[quadlist[i][0]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==8){
            vector<double> vert = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][2],l,quadlist[i][1]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][2],quadlist[i][3],quadlist[i][0],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][2]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][0]],vertlist[quadlist[i][1]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }else if (optres[i]==9){
            vector<double> vert = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            vertlist.push_back(vert);
            l = vertlist.size()-1;
            quaddelete.push_back(i);
            int tri1list[] = {quadlist[i][3],l,quadlist[i][2]}; vector<int> tri1 (tri1list, tri1list + sizeof(tri1list) / sizeof(int));
            int quad1list[] = {quadlist[i][3],quadlist[i][0],quadlist[i][1],l}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
            triadd.push_back(tri1);quadadd.push_back(quad1);
            vector<double> split1 = vertlist[quadlist[i][3]];
            vector<double> split2 = midPoint(vertlist[quadlist[i][1]],vertlist[quadlist[i][2]]);
            splitcoord.push_back(split1);splitcoord.push_back(split2);
        }
    }
    UpdateMesh updatemesh(splitcoord,triadd,quadadd,tridelete,quaddelete,vertlist);
    return updatemesh;
}

// based on the merging optimizaiton results, update the mesh
MergeOutput opt::mergeupdatemesh(vector<int> optres)
{
    unsigned int i, j, k;
    vector<vector<int> > edges = countedges();
    Indices indices = counttrisharededges();
    vector<vector<int> > quadadd;
    vector<int> tridelete;
    vector<vector<int> > trilist = triConnectivity;
    vector<vector<int> > quadlist = quadConnectivity;
    for (i=0;i<optres.size();i++){
        if (optres[i]>0.5){
            for (j=0;j<indices.triindex.size();j++){
                if (indices.edgeindex[j]==(int)i){
                    tridelete.push_back(indices.triindex[j][0]);tridelete.push_back(indices.triindex[j][1]);
                    int v0 = triConnectivity[indices.triindex[j][0]][0];int v1 = triConnectivity[indices.triindex[j][0]][1];int v2 = triConnectivity[indices.triindex[j][0]][2];
                    int e0 = edges[i][0];int e1 = edges[i][1];
                    // cout<<i<<": "<<"v0 "<<v0<<" v1 "<<v1<<" v2 "<<v2<<'\n';
                    // cout<<"e0 "<<e0<<" e1 "<<e1<<'\n';
                    if ((v0==e0 && v1==e1) || (v0==e1 && v1==e0)){
                        vector<int> tri2 = triConnectivity[indices.triindex[j][1]];
                        vector<int> edge = edges[i];
                        vector<int> diff;
                        sort(tri2.begin(), tri2.begin()+3);
                        sort(edge.begin(), edge.begin()+2);
                        set_difference(tri2.begin(), tri2.end(), edge.begin(), edge.end(),inserter(diff, diff.begin()));
                        int quad1list[] = {v2,v0,diff[0],v1}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
                        quadadd.push_back(quad1);
                        // cout<<'2';
                    }else if ((v1==e0 && v2==e1) || (v1==e1 && v2==e0)){
                        vector<int> tri2 = triConnectivity[indices.triindex[j][1]];
                        vector<int> edge = edges[i];
                        vector<int> diff;
                        sort(tri2.begin(), tri2.begin()+3);
                        sort(edge.begin(), edge.begin()+2);
                        set_difference(tri2.begin(), tri2.begin()+3, edge.begin(), edge.begin()+2,inserter(diff, diff.begin()));
                        int quad1list[] = {v0,v1,diff[0],v2}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
                        quadadd.push_back(quad1);
                        // cout<<'0';
                    }else if ((v2==e0 && v0==e1) || (v2==e1 && v0==e0)){
                        vector<int> tri2 = triConnectivity[indices.triindex[j][1]];
                        vector<int> edge = edges[i];
                        vector<int> diff;
                        sort(tri2.begin(), tri2.begin()+3);
                        sort(edge.begin(), edge.begin()+2);
                        set_difference(tri2.begin(), tri2.end(), edge.begin(), edge.end(),inserter(diff, diff.begin()));
                        int quad1list[] = {v1,v2,diff[0],v0}; vector<int> quad1 (quad1list, quad1list + sizeof(quad1list) / sizeof(int));
                        quadadd.push_back(quad1);
                        // cout<<'1';
                    }
                }
            }
        }
    }
    sort(tridelete.begin(),tridelete.end(),greater<int>());
    for (i=0;i<tridelete.size();i++){
        trilist.erase(trilist.begin()+tridelete[i]);
    }
    for (i=0;i<quadadd.size();i++){
        quadlist.push_back(quadadd[i]);
    }
    MergeOutput mergeoutput(trilist,quadlist);
    return mergeoutput;
}

// old way of creating triangulation (no longer useful)
vector<vector<int> > opt::connectTFI(vector<vector<double> > boundaries, vector<vector<double> > vertices, vector<int> num_b, vector<int> num_v, vector<int> index_b, vector<int> index_v)
{
    unsigned int i, j, k;
    vector<vector<int> > tris;
    for (i=0;i<4;i++){
        vector<vector<double> > l; l.assign(boundaries.begin()+num_b[i],boundaries.begin()+num_b[i+1]);
        vector<vector<double> > v; v.assign(vertices.begin()+num_v[i],vertices.begin()+num_v[i+1]);
        vector<int> indexl; indexl.assign(index_b.begin()+num_b[i],index_b.begin()+num_b[i+1]);
        vector<int> indexv; indexv.assign(index_v.begin()+num_v[i],index_v.begin()+num_v[i+1]);
        if (l.size()>v.size()){
            // creating tris from boundary
            for (j=0;j<l.size()-1;j++){
                double d = 1e8;
                vector<double> lj = l[j];
                int indexvk = 0;
                for (k=0;k<v.size();k++){
                    vector<double> vk = v[k];
                    double dist = distBetweenPoints(lj,vk);
                    if (dist<d){
                        d = dist;
                        indexvk = indexv[k];
                    }
                }
                if (i==0 || i==3){
                    int trigroup[] = {indexl[j],indexvk,indexl[j+1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }else{
                    int trigroup[] = {indexl[j],indexl[j+1],indexvk}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }
            }

            // creating tris form vertices
            for (j=0;j<v.size()-1;j++){
                vector<double> vj0 = v[j];
                vector<double> vj1 = v[j+1];
                int indexlk = 0;
                for (k=0;k<l.size();k++){
                    double dist0 = distBetweenPoints(vj0,l[k]);
                    double dist1 = distBetweenPoints(vj1,l[k]);
                    if (dist1<dist0){
                        indexlk = indexl[k];
                        break;
                    }
                }
                if (i==0 || i==3){
                    int trigroup[] = {indexv[j],indexv[j+1],indexlk}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }else{
                    int trigroup[] = {indexv[j],indexlk,indexv[j+1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }
                
            }
        }else if (l.size()==v.size()){
            for (j=0;j<l.size()-1;j++){
                if (i==0 || i==3){
                    int tri1group[] = {indexl[j],indexv[j+1],indexl[j+1]}; vector<int> tri1 (tri1group, tri1group + sizeof(tri1group) / sizeof(int));
                    int tri2group[] = {indexv[j],indexv[j+1],indexl[j]}; vector<int> tri2 (tri2group, tri2group + sizeof(tri2group) / sizeof(int));
                    tris.push_back(tri1); tris.push_back(tri2);
                }else{
                    int tri1group[] = {indexl[j],indexl[j+1],indexv[j+1]}; vector<int> tri1 (tri1group, tri1group + sizeof(tri1group) / sizeof(int));
                    int tri2group[] = {indexv[j],indexl[j],indexv[j+1]}; vector<int> tri2 (tri2group, tri2group + sizeof(tri2group) / sizeof(int));
                    tris.push_back(tri1); tris.push_back(tri2);
                }
            }
        }else{
            for (j=0;j<v.size()-1;j++){
                double d = 1e8;
                vector<double> vj = v[j];
                int indexlk = 0;
                for (k=0;k<l.size();k++){
                    vector<double> lk = l[k];
                    double dist = distBetweenPoints(lk,vj);
                    if (dist<d){
                        d = dist;
                        indexlk = indexl[k];
                    }
                }
                if (i==0 || i==3){
                    int trigroup[] = {indexv[j],indexv[j+1],indexlk}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }else{
                    int trigroup[] = {indexv[j],indexlk,indexv[j+1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                    tris.push_back(tri);
                }
            }
            // creating tris from boundary
            for (j=0;j<l.size()-1;j++){
                vector<double> lj0 = l[j];
                vector<double> lj1 = l[j+1];
                int indexvk = 0;
                for (k=0;k<v.size();k++){
                    double dist0 = distBetweenPoints(lj0,v[k]);
                    double dist1 = distBetweenPoints(lj1,v[k]);
                    if (dist1<dist0){
                        indexvk = indexv[k];
                        break;
                    }
                }                
                if ((((int) j)==((int) l.size()-2)) && (indexvk==0)){
                    if (i==0 || i==3){
                        int trigroup[] = {indexl[l.size()-2],indexv[v.size()-1],indexl[l.size()-1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                        tris.push_back(tri);
                    }else{
                        int trigroup[] = {indexl[l.size()-2],indexl[l.size()-1],indexv[v.size()-1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                        tris.push_back(tri);
                    }
                }else{
                    if (i==0 || i==3){
                        int trigroup[] = {indexl[j],indexvk,indexl[j+1]}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                        tris.push_back(tri);
                    }else{
                        int trigroup[] = {indexl[j],indexl[j+1],indexvk}; vector<int> tri (trigroup, trigroup + sizeof(trigroup) / sizeof(int));
                        tris.push_back(tri);
                    }
                }
            }
        }
    }
    return tris;
}

// get the average size of all edges in the mesh
double opt::averagesize()
{
    unsigned int i;
    vector<vector<int> > edges = countedges();
    double dist = 0;
    for (i=0;i<edges.size();i++){
        vector<double> e0 = vertexCoords[edges[i][0]];
        vector<double> e1 = vertexCoords[edges[i][1]];
        double disti = distBetweenPoints(e0,e1);
        dist += disti;
    }
    dist = dist/edges.size();
    return dist;
}

// determine intersection between two triangulations
TriInt opt::member_intersection(vector<vector<int> > itris, vector<vector<int> > jtris, vector<vector<double> > icoords, vector<vector<double> > jcoords, vector<int> ifixed, vector<int> jfixed)
{
    unsigned int i, j;
    vector<vector<double> > shared;
    cout<<itris.size()<<' '<<jtris.size()<<' '<<icoords.size()<<' '<<jcoords.size()<<' '<<ifixed.size()<<' '<<jfixed.size()<<'\n';
    for (i=0;i<itris.size();i++){
        // three vertices in ith tiranlge
        vector<double> itri0 = icoords[itris[i][0]];
        vector<double> itri1 = icoords[itris[i][1]];
        vector<double> itri2 = icoords[itris[i][2]];
        for (j=0;j<jtris.size();j++){
            vector<double> jtri0 = jcoords[jtris[j][0]];
            vector<double> jtri1 = jcoords[jtris[j][1]];
            vector<double> jtri2 = jcoords[jtris[j][2]];
            // intersection test
            // (jtri0, jtri1) and itri
            if ((signedvolume(jtri0,itri0,itri1,itri2)*signedvolume(jtri1,itri0,itri1,itri2)<=1e-3) &&
            (((signedvolume(jtri0,jtri1,itri0,itri1)>=1e-3) && (signedvolume(jtri0,jtri1,itri1,itri2)>=1e-3) && (signedvolume(jtri0,jtri1,itri2,itri0)>=1e-3)) ||
            ((signedvolume(jtri0,jtri1,itri0,itri1)<=1e-3) && (signedvolume(jtri0,jtri1,itri1,itri2)<=1e-3) && (signedvolume(jtri0,jtri1,itri2,itri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(jtri0,jtri1,itri0,itri1,itri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            // (jtri1, jtri2) and itri
            }else if((signedvolume(jtri1,itri0,itri1,itri2)*signedvolume(jtri2,itri0,itri1,itri2)<=1e-3) &&
            (((signedvolume(jtri1,jtri2,itri0,itri1)>=1e-3) && (signedvolume(jtri1,jtri2,itri1,itri2)>=1e-3) && (signedvolume(jtri1,jtri2,itri2,itri0)>=1e-3)) ||
            ((signedvolume(jtri1,jtri2,itri0,itri1)<=1e-3) && (signedvolume(jtri1,jtri2,itri1,itri2)<=1e-3) && (signedvolume(jtri1,jtri2,itri2,itri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(jtri1,jtri2,itri0,itri1,itri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            // (jtri2, jtri0) and itri
            }else if((signedvolume(jtri2,itri0,itri1,itri2)*signedvolume(jtri0,itri0,itri1,itri2)<=1e-3) &&
            (((signedvolume(jtri2,jtri0,itri0,itri1)>=1e-3) && (signedvolume(jtri2,jtri0,itri1,itri2)>=1e-3) && (signedvolume(jtri2,jtri0,itri2,itri0)>=1e-3)) ||
            ((signedvolume(jtri2,jtri0,itri0,itri1)<=1e-3) && (signedvolume(jtri2,jtri0,itri1,itri2)<=1e-3) && (signedvolume(jtri2,jtri0,itri2,itri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(jtri2,jtri0,itri0,itri1,itri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            // (itri0, itri1) and jtri
            }else if((signedvolume(itri0,jtri0,jtri1,jtri2)*signedvolume(itri1,jtri0,jtri1,jtri2)<=1e-3) &&
            (((signedvolume(itri0,itri1,jtri0,jtri1)>=1e-3) && (signedvolume(itri0,itri1,jtri1,jtri2)>=1e-3) && (signedvolume(itri0,itri1,jtri2,jtri0)>=1e-3)) ||
            ((signedvolume(itri0,itri1,jtri0,jtri1)<=1e-3) && (signedvolume(itri0,itri1,jtri1,jtri2)<=1e-3) && (signedvolume(itri0,itri1,jtri2,jtri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(itri0,itri1,jtri0,jtri1,jtri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            // (itri1, itri2) and jtri
            }else if((signedvolume(itri1,jtri0,jtri1,jtri2)*signedvolume(itri2,jtri0,jtri1,jtri2)<=1e-3) &&
            (((signedvolume(itri1,itri2,jtri0,jtri1)>=1e-3) && (signedvolume(itri1,itri2,jtri1,jtri2)>=1e-3) && (signedvolume(itri1,itri2,jtri2,jtri0)>=1e-3)) ||
            ((signedvolume(itri1,itri2,jtri0,jtri1)<=1e-3) && (signedvolume(itri1,itri2,jtri1,jtri2)<=1e-3) && (signedvolume(itri1,itri2,jtri2,jtri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(itri1,itri2,jtri0,jtri1,jtri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            // (itri2, itri0) and jtri
            }else if((signedvolume(itri2,jtri0,jtri1,jtri2)*signedvolume(itri0,jtri0,jtri1,jtri2)<=1e-3) &&
            (((signedvolume(itri2,itri0,jtri0,jtri1)>=1e-3) && (signedvolume(itri2,itri0,jtri1,jtri2)>=1e-3) && (signedvolume(itri2,itri0,jtri2,jtri0)>=1e-3)) ||
            ((signedvolume(itri2,itri0,jtri0,jtri1)<=1e-3) && (signedvolume(itri2,itri0,jtri1,jtri2)<=1e-3) && (signedvolume(itri2,itri0,jtri2,jtri0)<=1e-3)))){
                vector<double> interpt = intersectionpoint(itri2,itri0,jtri0,jtri1,jtri2);
                if (vectNorm(interpt)>1e-5){
                    shared.push_back(interpt);
                }
                // cout<<i<<' '<<j<<'\n';
            }             
        }
    }

    cout<<i<<' '<<j<<' '<<'\n';
    if (shared.size()>0){
        cout<<"shared:\n";
        for (i=0;i<shared.size();i++){
            cout<<shared[i][0]<<' '<<shared[i][1]<<' '<<shared[i][2]<<'\n';
        }
        sort(shared.begin(),shared.end());
        vector<double> e1 = shared[0];
        vector<double> e2 = shared[shared.size()-1];
        shared.erase(unique(shared.begin(),shared.end()),shared.end());
        vector<vector<double> > sorted = sortshared(shared);
        Reconnection itrint = reconnect(sorted,itris,icoords,ifixed);
        Reconnection jtrint = reconnect(sorted,jtris,jcoords,jfixed);
        itris = itrint.trilist;
        icoords = itrint.vertlist;
        ifixed = itrint.projlist;
        jtris = jtrint.trilist;
        jcoords = jtrint.vertlist;
        jfixed = jtrint.projlist;
    }
    TriInt triint(itris, jtris, icoords, jcoords, ifixed, jfixed);
    return triint;
}

// calculate intersection point between a line and a trianlge in 3d 
vector<double> opt::intersectionpoint(vector<double> q1, vector<double> q2, vector<double> p1, vector<double> p2, vector<double> p3)
{
    vector<double> p1p2 = vectSubtract(p2,p1);
    vector<double> p1p3 = vectSubtract(p3,p1);
    vector<double> N = crossProduct(p1p2,p1p3);
    vector<double> p1q1 = vectSubtract(q1,p1);
    vector<double> q1q2 = vectSubtract(q2,q1);
    double temp1 = - dotProduct(p1q1,N);
    double temp2 = dotProduct(q1q2,N);
    double t = temp1/temp2;
    vector<double> intersection = vectAdd(q1,q1q2,t);
    vector<double> intersectionnew(3,0.0);
    if (temp2==0){
        intersection.clear();
        intersection.assign(intersectionnew.begin(),intersectionnew.end());
    }
    return intersection;
}

// sort the shared edges to fix them easier (entire fixing algorithm should be replaced with edge_based method in the future)
vector<vector<double> > opt::sortshared(vector<vector<double> > shared)
{
    unsigned int i, j, k;
    // reverse(shared.begin(),shared.end());
    for (i=0;i<shared.size();i++){
        vector<double> veci = shared[i];
        vector<vector<double> > reduced; reduced.assign(shared.begin(),shared.end());
        reduced.erase(reduced.begin()+i);
        vector<double> v_prev, v_new;
        for (j=0;j<reduced.size();j++){
            vector<double> vecj = reduced[j];
            if (j==0){
                v_prev = vectSubtract(vecj,veci);
            }else{
                v_new = vectSubtract(vecj,veci);
                if (dotProduct(v_prev,v_new)<0){
                    break;
                }
                v_prev.assign(v_new.begin(),v_new.end());
            }
        }
        unsigned int jlimt = reduced.size();
        if (j==jlimt){
            break;
        }
    }
    vector<vector<double> > reduced; reduced.assign(shared.begin(),shared.end());
    reduced.erase(reduced.begin()+i);
    vector<double> start = shared[i];
    vector<vector<double> > sorted;
    sorted.push_back(start);
    vector<double> dists(reduced.size(),0.0);
    for (j=0;j<reduced.size();j++){
        dists[j] = distBetweenPoints(reduced[j],start);// calculate the distances between the intersection points and projected point
    }
    vector<double> reference; reference.assign(dists.begin(), dists.end());// save the distances as reference
    
    sort(dists.begin(),dists.end());// sort distance in the increasing order
    
    for (j=0;j<reduced.size();j++){// sort index by distances 
        for (k=0;k<reduced.size();k++){
            if (reference[k]==dists[j]){
                sorted.push_back(reduced[k]);
            }
        }
    }
    sorted.erase(unique(sorted.begin(),sorted.end()),sorted.end());
    sorted = uniquevet(sorted);
    cout<<"sorted size"<<sorted.size()<<'\n';
    for (i=0;i<sorted.size();i++){
        cout<<"sorted_"<<i<<' '<<sorted[i][0]<<' '<<sorted[i][1]<<' '<<sorted[i][2]<<'\n';
    }
    return sorted;
}

// get all the unique vertices
vector<vector<double> > opt::uniquevet(vector<vector<double> > sorted)
{
    unsigned int i, j;
    vector<int> index(sorted.size(),0);
    for (i=0;i<sorted.size();i++){
        vector<double> vecti = sorted[i];
        for (j=i+1;j<sorted.size();j++){
            vector<double> vectj = sorted[j];
            if (distBetweenPoints(vecti,vectj)<1e-2){
                index[j] = 1;
            }
        }
    }
    for (i=0;i<index.size();i++){
        cout<<"i="<<i<<'\n';
        cout<<"i: "<<index[i]<<'\n';
    }
    for (i=0;i<index.size();i++){
        j=index.size()-i-1;
        cout<<i<<'\n';
        cout<<j<<'\n';
        if (index[j]==1){
            sorted.erase(sorted.begin()+j);
        }
    }
    return sorted;
}

// retriangulation
Reconnection opt::reconnect(vector<vector<double> > sorted, vector<vector<int> > tris, vector<vector<double> > coords, vector<int> fixed)
{
    unsigned int i, j, k;
    vector<vector<double> > vertlist = coords;
    vector<int> projlist;
    vector<vector<int> > trilist = tris;
    for (i=0;i<sorted.size()-1;i++){// loop through all the projected points
        vector<int> projlisti;
        vector<vector<double> > inti;
        vector<int> V;
        // vector<double> mid = midPoint(projpoints[i],projpoints[i+1]);// find the midpoint of two projection points
        vector<double> midp = midPoint(sorted[i],sorted[i+1]);// find the midpoint of two projected points
        vector<double> projvector = vectSubtract(sorted[i],sorted[i+1]);
        vector<vector<int> > edges = counttriedges(trilist);
        vector<double> e0, e1;
        vector<double> p0 = sorted[i];
        vector<double> p1 = sorted[i+1];
        int count = 0;
        if (i==0){
            double distttt = 10000;
            unsigned int indexxx=0;
            for (j=0;j<vertlist.size();j++){
                if (distBetweenPoints(p0,vertlist[j])<distttt){
                    distttt = distBetweenPoints(p0,vertlist[j]);
                    indexxx=j;
                }
            }
            cout<<"================================================================================"<<distttt<<'\n';
            projlisti.push_back(indexxx);
            
        }
        vector<double> revprojdir(3,1.0);
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
                        if (distBetweenPoints(intersect,vertlist[k])<1e-2){
                            if (k!=vertlist.size()-2){
                                V.push_back(k);
                                cout<<"----------------------------------------------------------Yeah! It worked!----------------------------------------------------------\n";
                                cout<<"k:"<<k<<'\n';
                                cout<<"intersect:"<<intersect[0]<<' '<<intersect[1]<<' '<<intersect[2]<<'\n';
                                cout<<"vertlist.size"<<vertlist.size()<<'\n';
                                break;
                            }
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
        for (j=0;j<V.size();j++){
            // cout<<"V: "<<V[j]<<'\n';
        }
        if (i==sorted.size()-2){// at the second last point, add the last point
            int count = 0;
            for (k=0;k<vertlist.size();k++){
                if (distBetweenPoints(sorted[i+1],vertlist[k])<1e-2){
                    V.push_back(k);
                    break;
                }
                if (k==vertlist.size()-1){
                    count = 1;
                }
            }
            if (count==1){
                vertlist.push_back(sorted[i+1]);
                projlisti.push_back(vertlist.size()-1);
            }
        }
        projlist.insert(projlist.end(),projlisti.begin(),projlisti.end());
        // update trilist
        // cout<<"update trilist"<<'\n';
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
                    if((abs(area1-area)<1e-6) || (abs(area2-area)<1e-5) || (abs(area3-area)<1e-5)){
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
            if (tripush.size()>0){
                for (k=0;k<tripush.size();k++){
                    trilist.push_back(tripush[k]);
                }
            }
        }
    }
    fixed.insert(fixed.end(),projlist.begin(),projlist.end());
    vector<int> vertindex(vertlist.size(),-1);
    vector<double> verti, vertj;
    vector<vector<double> > univertlist;
    k=0;
    for (i=0;i<vertlist.size();i++){
        if (vertindex[i]<0){
            vertindex[i] = (int)k;
            verti = vertlist[i];
            univertlist.push_back(verti);
            for (j=i;j<vertlist.size();j++){
                vertj = vertlist[j];
                if (distBetweenPoints(verti,vertj)<1e-5){
                    vertindex[j] = (int)k;
                }
            }
            k++;
        }
    }
    for (i=0;i<trilist.size();i++){
        for (j=0;j<3;j++){
            trilist[i][j] = vertindex[trilist[i][j]];
        }
    }
    for (i=0;i<fixed.size();i++){
        fixed[i] = vertindex[fixed[i]];
    }
    int count = 0;
    for (i=0;i<fixed.size()-1;i++){
        if (count==1){
            i--;
            count = 0;
        }
        
        vector<double> vert_start = univertlist[fixed[i]];
        vector<double> vert_end = univertlist[fixed[i+1]];
        for (j=i+2;j<fixed.size();j++){
            vertj = univertlist[fixed[j]];
            if ((abs((distBetweenPoints(vert_start,vertj)+distBetweenPoints(vert_end,vertj))-distBetweenPoints(vert_start,vert_end))<1e-5) && (distBetweenPoints(vertj,vert_start)>1e-2) && (distBetweenPoints(vertj,vert_end)>1e-2)){
                fixed.insert(fixed.begin()+i+1,fixed[j]);
                count = 1;
                break;
            }
        }

    }
    Reconnection reconnection(univertlist,trilist,fixed);
    return reconnection;
}

// for a given trilist, count all the edges in it
vector<vector<int> > opt::counttriedges(vector<vector<int> > trilist)
{   
    unsigned int i;
    vector<int> v1, v2, v3;
    vector<vector<int> > edges;
    for (i=0;i<trilist.size();i++){
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

// remove duplicate trianlges if there is any
vector<int> opt::removeduplicatetris()
{
    vector<vector<int> > trilist = triConnectivity;
    unsigned int i, j;
    vector<int> output(trilist.size(),-1);
    for (i=0;i<trilist.size();i++){
        if (output[i]<0){
            output[i] = (int) i;
            int trii0 = trilist[i][0];
            int trii1 = trilist[i][1];
            int trii2 = trilist[i][2];
            // if ((int)i == 15697){
            //     cout<<trii0<<' '<<trii1<<' '<<trii2<<'\n';
            // }
            for (j=i+1;j<trilist.size();j++){
                if (output[j]<0){
                    int trij0 = trilist[j][0];
                    int trij1 = trilist[j][1];
                    int trij2 = trilist[j][2];
                    // if ((int)j == 15698){
                    //     cout<<trij0<<' '<<trij1<<' '<<trij2<<'\n';
                    // }
                    if ( ((trii0==trij0)&&(trii1==trij1)&&(trii2==trij2)) || ((trii0==trij0)&&(trii1==trij2)&&(trii2==trij1))
                    || ((trii0==trij1)&&(trii1==trij2)&&(trii2==trij0)) || ((trii0==trij1)&&(trii1==trij0)&&(trii2==trij2))
                    || ((trii0==trij2)&&(trii1==trij1)&&(trii2==trij0)) || ((trii0==trij2)&&(trii1==trij0)&&(trii2==trij1))){
                        // cout<<"fail\n";
                        output[j] = (int) i;
                    }
                }
            }
        }
    }
    return output;
}

// used for saving all the meshes into one vtk file, the same points must be detected
Uni opt::mergeduppts(double tol)
{
    unsigned int i, j;
    vector<vector<int> > quadlist = quadConnectivity;
    vector<vector<double> > vertlist = vertexCoords; 
    vector<vector<double> > univert;
    vector<int> vertindex(vertlist.size(),-1);
    int count =0;
    for (i=0;i<vertlist.size();i++){
        cout<<i<<" "<<vertlist.size()<<"\n";
        vector<double> verti = vertlist[i];
        if (vertindex[i]<0){
            univert.push_back(verti);
            vertindex[i] = count;
            for (j=i+1;j<vertlist.size();j++){
                vector<double> vertj = vertlist[j];
                if (distBetweenPoints(verti,vertj)<tol){
                    vertindex[j] = count;
                }
            }
            count++;
        }
        
    }
    for (i=0;i<quadlist.size();i++){
        for (j=0;j<4;j++){
            quadlist[i][j] = vertindex[quadlist[i][j]];
        }
    }
    Uni uni(univert, quadlist);
    return uni;
}
// ------------------------------------------------------------------------------------------------
// useful functions
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

vector<double> vectSubtract(vector<double> &a, vector<double> &b)
{
    // Vector subtraction, c=a-b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i]-b[i];
    return c;
}

vector<double> midPoint(vector<double> &a, vector<double> &b)
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> m(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        m[i] = 0.5*a[i] + 0.5*b[i];
    return m;
}

vector<double> midPoint(vector<double> &a, vector<double> &b, vector<double> &c)
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> m(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        m[i] = (a[i] + b[i] + c[i])/3;
    return m;
}

vector<double> midPoint(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d)
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> m(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        m[i] = (a[i] + b[i] + c[i] + d[i])/4;
    return m;
}

double angleOfVecs(vector<double> &a, vector<double> &b)
{
    double PI = atan(1)*4;
    double angle = acos(dotProduct(a,b)/(vectNorm(a)*vectNorm(b)))/PI*180;
    return angle;
} 

double angleRatio(vector<double> &a, vector<double> &b, double theta0)
{
    double PI = atan(1)*4;
    double angleratio = acos(dotProduct(a,b)/(vectNorm(a)*vectNorm(b)))/PI*180/theta0;
    return angleratio;
} 

double dotProduct(vector<double> &a, vector<double> &b)
{
    // Calculates dot product v1.v2
    double dp=0.0;
    for (unsigned int i=0; i<3; i++)
        dp += a[i]*b[i]; 
    return dp;
}

double vectNorm(vector<double> &a)
{
    double norm = sqrt(dotProduct(a,a));
    return norm;
}

void arrayadd(int arr[], int size, int val)
{   
    // cout<<val<<" val\n";
    for (int i=0;i<size;i++){
        arr[i] +=val;
    }

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

vector<double> normalize(vector<double> &v)
{   
    double sum=0.0;
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

double distBetweenPoints(vector<double> &p1, vector<double> &p2)
{
    // Calculate the distance between points p1 and p2, |p1-p2|
    double sum=0.0;
    for (unsigned int i=0; i<3; i++)
        sum += pow((p1[i]-p2[i]),2.0);
    sum = sqrt(sum);
    return sum;
}

vector<double> vectScale(vector<double> &p, double s)
{
    vector<double> vect(3);
    for (unsigned int i=0;i<3;i++){
        vect[i] = s*p[i];
    }
    return vect;
}

vector<double> vectMean(vector<vector<double> > p)
{
    vector<double> vect(3);
    for (unsigned int i=0;i<3;i++){
        for (unsigned int j=0;j<p.size();j++){
            vect[i] += p[j][i];
        }
    }
    // cout<<vect.size();
    // vector<double> vect1(3,0.0);
    return vect;
}

double signedvolume(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
    vector<double> ab = vectSubtract(b,a);
    vector<double> ac = vectSubtract(c,a);
    vector<double> ad = vectSubtract(d,a);
    vector<double> temp = crossProduct(ab,ac);
    double signedvolume = 1./6.*dotProduct(temp,ad);
    return signedvolume;
}