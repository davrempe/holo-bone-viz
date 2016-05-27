////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VOLUMETRIC ANALYSIS OF SURFACE PROPERTIES (VASP)                                                           //
// BRIAN Y CHEN, HONIG LAB, 2010                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// SurfaceObject.h: interface for the SurfaceObject class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SURFACEOBJECT_H
#define SURFACEOBJECT_H

#include "defs.h"
#include "set.h"

class SurfaceObject  
{
public:
	int numPoints;
	int numTriangles;
	
	///mandatory memory allocation
	double * surfacePoints;
	double * surfaceNormals;		///these averaged normals are used for rendering
	int * triangles;
	double * triangleNormals;	///these normals are used for flat shading and interior checking
	double * centroid;

	///optional memory allocation
	set_t highlights;
	set_t edges; 				///edge highlights
	double * colors;			///optional colors, initialized only through separate request.  One vector for each pt.

	///////////////////////////////////////////////////////
	// Constructors / Destructors
	SurfaceObject();
	SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris);
	SurfaceObject(set_t pts, set_t tris);
	virtual ~SurfaceObject();
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///returns the coords of the requested triangle
	///returns null if out of range.
	///notice that the triangle must be formatted in this way in the array (3 then 3 then 3)
	double * getTriangle(int num); 
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///Computes the centroid based on the point positions
	double * getCentroid(); 
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///This function flips the normals of the surfaceObject backwards 
	///so that we can treat it as a negative volume
	void flipNormals();
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///Old Constructors / Destructors
	/// SurfaceObject();
	/// SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris, bool elimIntCavities);
	/// SurfaceObject(int numPts, double * pts, double * norms, int numTris, int * tris, double * triangleNorms);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///Adds the geometry of another SurfaceObject to this SurfaceObject.
	///ASSUMES THAT THIS SURFACEOBJECT DOES NOT COLLIDE WITH THE OTHER SURFACE OBJECT.
	///DOES NOT CHECK FOR COLLISION, SO MUST BE DONE BEFORE THIS FUNCTION IS CALLED.
	/// void addObject(SurfaceObject * obj);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///Copies a surfaceObject
	SurfaceObject * copy(void);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///add colors
	/// void addColors(double * c);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///print out details
	/// void toString();
	/// void printSummary();
	///////////////////////////////////////////////////////
	
	//int countNumCommon(set_t set1, set_t set2);
};

#endif 



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VOLUMETRIC ANALYSIS OF SURFACE PROPERTIES (VASP)                                                           //
// BRIAN Y CHEN, HONIG LAB, 2010                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


