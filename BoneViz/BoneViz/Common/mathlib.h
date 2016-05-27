//////////////////////////////////////////////////////////////////////
////Header File for Math Library
////
///////////////////////////////


#if !defined(_MATHLIB_H_)
#define _MATHLIB_H_


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


double * transformVector(double * matrix, double * vector);
int * sort3(int first, int second, int third);				///Make this better
double * transpose(double * mat);
double * transpose4x4(double * mat);
double * invertMatrix(double * mat);	
double * normalizeVector(double * vector);
double * normalizeVector(double v1, double v2, double v3);
void normalizeVector(double * vector, double * output);
double vectorSize(double x, double y, double z);
double vectorSize2d(double x, double y);
double vectorSize(double * A, double * B);
double * crossProd(double * first, double * second);
double * getCentroid(double * vertices, int size);
double * matrixMult(double * first, double * second);
double determinant3x3(double * matrix);
double determinant4x4(double * m);
double determinant5x5(double * m);
double * scalarMultiplication(double * vector, double scale);
double dotProd(double * first, double * second);
char * stringCopy(char* str, int start, int end);
char * firstChar(char * str);
bool makeBool(int * i);
double getRand(double bound);
double SSStriangleAngle(double a, double b, double c);
void printMatrix(double * mat, int size);
double * transformVector3x4(double * mat, double * vector);		////transforms a 3-vector with a 4matrix, by padding a 1 for scale.
double * matrixMult4x4(double * first, double * second);		///4x4 matrix times 4x4 matrix multiplication
bool testMatrix(double * matrix);
bool testMatrix3x3(double * matrix);
double triangleArea(double * t1, double * t2, double * t3);
double triangleArea(double * points, int a, int b, int c);
double triangleArea(double * triangle);
bool intervalOverlapTest(double lo, double hi, double a, double b);	///determines if an interval a,b (a not necessarily < b) overlaps [lo,hi].


///return 2 if all even, 1 if all odd, 0 of inconsistent, -1 if size = 0;
int parityTest(int * array, int size);

//Projecting triangle t1, t2, t3 onto the unit sphere around pt, this 
//function measures the area of the spherical triangle.
double sphereicalTriangleArea(double * pt, double * t );


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///input is four points in 3d space, implemented volume on the method by Niccolò Fontana Tartaglia, which
///is a generalization of Heron's area method:
//                 [  0     1    1     1     1    ]
//                 [  1     0   d12^2 d13^2 d14^2 ]
// 288 * v^2 = det [  1   d21^2  0    d23^2 d24^2 ]
//                 [  1   d31^2 d32^2  0    d34^2 ]
//                 [  1   d41^2 d42^2 d43^2  0    ]
double computeTetrahedralVolume( double * p1, double * p2, double * p3, double * p4 );




/*----------------Macros----------------*/
///Cross Product
#ifndef CROSS
#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#endif

///Dot Product 
#ifndef DOT
#define DOT(v1,v2) ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]))
#endif



#endif


