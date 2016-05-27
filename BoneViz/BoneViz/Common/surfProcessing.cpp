////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VOLUMETRIC ANALYSIS OF SURFACE PROPERTIES (VASP)                                                           //
// BRIAN Y CHEN, HONIG LAB, 2010                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// surfProcessing.cpp : Source file for surfProcessing.h
// Reads files in the SURF format and generates a SurfaceObject,
// Writes files in the SURF format, using a SurfaceObject as input.
//

#include "pch.h"
#include "surfProcessing.h"


// #################################################################################################################
// ###    ##  ##  #     ##     ###     #    #  ###     ###     ####   ###     ###    ##    #  ##  ##    ############
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  ##     ##  ##  #  ##  ##  ##   #  #  ##  ###########
// ##  #####  ##  #  ##  #  ######  #####  ##  ###  ######  ##  ##  #  ##  ##  #  ######  ##      #  ###############
// ###  ####  ##  #    ###    ####    ###  ##  ###    ####  ##  #  ###  #    ####  #####  ##  #   #  ###############
// #####  ##  ##  #  #  ##  ######  #####  ##  ###  ######     ##       #  #  #####  ###  ##  ##  #  #   ###########
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  #####  ###  #  ##  #  ##  ##  ##  ##  #  ##  ###########
// ##  ##  #      #  ##  #  ######  #####  ##  ###  ######  #####  ###  #  ##  #  ##  ##  ##  ##  #  ##  ###########
// ###    ###    ##  ##  #  ######  ####    #    #     ###  #####  ###  #  ##  ##    ##    #  ##  ##    ############
// #################################################################################################################

/////////////////////////////////////////////////////////
// Top level parsing control.  Run this function to parse a given file name into 
// a SurfaceObject.  It calls the functions below.
/////////////////////////////////////////////////////////

SurfaceObject * parseGeometryFile(char * fileName)
{
	FILE * surfFile = fopen(fileName, "r");
	if(surfFile == NULL){
		printf("ERROR: Surface file [%s] not found.  Exitting.\n", fileName);
		exit(1);
	}
	else{
		printf("Parsing Surface Geometry File [%s] ...", fileName); fflush(stdout);
	}

	int numPoints = getNumberOfPoints(surfFile);
	double * VectorsAndNormals = getSurfaceGeometryAndNormals(surfFile, numPoints);
	int numTriangles = getNumberOfTriangles(surfFile);
	int * Topology = getTopology(surfFile, numTriangles);
	double * colors = getColors(surfFile, numPoints);
	
//	This is removed for VASP (it remains in surfaceExtractor) because we do not do this kind
//	of filtering in VASP.  (It should be dealt with in SurfaceExtractor before giving to VASP)
//
//	bool eliminateInteriorCavities = false;
//	bool eliminateInteriorCavities = true;
//	SurfaceObject * result = new SurfaceObject(numPoints, VectorsAndNormals, numTriangles, Topology, eliminateInteriorCavities);
	SurfaceObject * result = new SurfaceObject(numPoints, VectorsAndNormals, numTriangles, Topology);
	
	//result->addColors(colors);
	
	delete[](VectorsAndNormals);
	delete[](Topology);
	delete[](colors);

	fclose(surfFile);

	printf(" done.\n");
	return result;
}


/////////////////////////////////////////////////////////
// find out how many points on the surface
/////////////////////////////////////////////////////////
int getNumberOfPoints(FILE * surfaceFile)
{
	char * line = new char[1000];
	char * ptr;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "GEOMETRY:") == 0){
			break;
		}
		else{
			printf("ERROR: Expected GEOMETRY: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}

	////if the loop has ended, then we are at the GEOMETRY: line.
	////find out how many points there are.
	ptr = strtok(NULL, " \n\t\r");

	int result = atoi(ptr);

	delete[](line);

	return result;
}


/////////////////////////////////////////////////////////
// get the vectors and surface normals on the surface
/////////////////////////////////////////////////////////
double * getSurfaceGeometryAndNormals(FILE * surfaceFile, int numPts)
{
	int i = 0;
	char * line = new char[1000];
	char * ptr;
	double * result = new double[6*numPts];
	for(i = 0; i<numPts; i++){
		fgets(line, 1000, surfaceFile);

/*
		///if the line is a comment, warn and continue;
		bool iscomment = false;
		if(line[0] == '#'){ iscomment = true; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ iscomment = true; }
		if(iscomment){
			printf("WARNING: COMMENTS DETECTED INSIDE GEOMETRY: DEFINITION.  IGNORING\n");
			i = i - 1;
			continue;
		}
*/

		///otherwise its not a comment, so parse it.  Use linecopy.
		ptr = strtok(line, " \n\t\r");
		result[6*i+0] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+1] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+2] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+3] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+4] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+5] = atof(ptr);
	}

	delete[](line);

	return result;
}


/////////////////////////////////////////////////////////
// get the number of triangles
/////////////////////////////////////////////////////////
int getNumberOfTriangles(FILE * surfaceFile)
{
	char * line = new char[1000];
	char * ptr;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "TOPOLOGY:") == 0){
			break;
		}
		else{
			printf("ERROR: Expected GEOMETRY: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}

	////if the loop has ended, then we are at the GEOMETRY: line.
	////find out how many points there are.
	ptr = strtok(NULL, " \n\t\r");

	int result = atoi(ptr);

	delete[](line);

	return result;
}


/////////////////////////////////////////////////////////
// get the triangles
/////////////////////////////////////////////////////////
int * getTopology(FILE * surfaceFile, int numTris)
{
	int i = 0;
	char * line = new char[1000];
	char * ptr;
	int * result = new int[3*numTris];
	for(i = 0; i<numTris; i++){
		fgets(line, 1000, surfaceFile);

/*
		///if the line is a comment, warn and continue;
		bool iscomment = false;
		if(line[0] == '#'){ iscomment = true; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ iscomment = true; }
		if(iscomment){
			printf("WARNING: COMMENTS DETECTED INSIDE GEOMETRY: DEFINITION.  IGNORING\n");
			i = i - 1;
			continue;
		}
*/

		///otherwise its not a comment, so parse it.  Use linecopy.
		ptr = strtok(line, " \n\t\r");
		result[3*i+0] = atoi(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+1] = atoi(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+2] = atoi(ptr);
	}

	delete[](line);

	return result;	


}


/////////////////////////////////////////////////////////
// Gets the colors out of the SURF file.
/////////////////////////////////////////////////////////
double * getColors(FILE * surfaceFile, int numColors)
{
	double * result = new double[3*numColors];
	int i = 0;

	char * line = new char[1000];
	char * ptr;

	bool foundColors = false;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "COLORS:") == 0){
			foundColors = true;
			break;
		}
		else{
			printf("ERROR: Expected COLORS: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}
	
	if(!foundColors){
		delete[](result);
		delete[](line);
		return NULL;
	}
	

	///parse the colors
	for(i = 0; i<numColors; i++){
		fgets(line, 1000, surfaceFile);
		ptr = strtok(line, " \n\t\r");
		result[3*i+0] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+1] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+2] = atof(ptr);
	}
	
	delete[](line);
	return result;
}





// #################################################################################################################
// ###    ##  ##  #     ##     ###     #    #  ###     ####    ##  ##  #      #     ##  ##  #      #################
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  ##  #  ##  ###  ###################
// ##  #####  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  ##  #  ##  ###  ###################
// ###  ####  ##  #    ###    ####    ###  ##  ###    ####  ##  #  ##  ###  ###     ##  ##  ###  ###################
// #####  ##  ##  #  #  ##  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ##  ##  #      #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ###    ###    ##  ##  #  ######  ####    #    #     ####    ###    ####  ###  ######    ####  ###################
// #################################################################################################################

////generateSURF:  This function generates a SURF file based on a surfaceObject and a set of
////                  point indices that reference a patch on the surface.  The output is a
////                  SURF file that contains only the points relevant to render the patch.
////                  If the set of points is NULL, the entire surface is returned.
////                  
////                  Inputs:
////                     a SurfaceObject *			describes the surface
////                     set_t points				describes which vertices we want (NULL for all)
////                     char * outputFileName         filename of the output surface file
////
////                  Outputs:
////                     Void - only a file is generated
////                  
////                  =========================================================
////                    OUTPUT FORMAT:
////                       -Anywhere in the file, empty lines and line starting with '#' are
////                       considered comments, and are ignored.
////                       -File has 2 sections: geometry and topology
////                  
////                       The geometry section is always first, and starts with:
////                       GEOMETRY: <int>
////                       where <int> is the number of points to be specified.
////                       then there are <int> lines as follows:
////                       <float> <float> <float> <float> <float> <float>
////                       which stand for x,y,z, and xnormal, ynormal, znormal.
////                  
////                       The topology section is always second, and starts with:
////                       TOPOLOGY: <int>
////                       where <int> specifies the number of triangles on the surface
////                       Then there are <int> lines as follows:
////                       <int> <int> <int>
////                       Each int stands for the 3 corners of the triangle
////                       and is an index into the array of points provided in the
////                       Geometry section.  I.e. the geometry section is indexed
////                       starting at zero, and ending at (size-1).
////                  =========================================================
////                  
void generateSURF(SurfaceObject * surf, set_t points, char * outputFileName)
{
	int i = 0;
	char * tempString = new char[300];
	FILE * currentOutputFile = fopen(outputFileName, "w");
	
	int numVerts = surf->numPoints;
	if(points != NULL){
		numVerts = size_set(points);
	}

	if(points != NULL){
		printf("Generating SURF file [%s]: numVerts overall: %i, numVerts selected for output: %i\n", outputFileName, surf->numPoints, size_set(points));
	}
	else{
		printf("Generating SURF file [%s]: numVerts overall: %i, numVerts selected for output: %i\n", outputFileName, surf->numPoints, surf->numPoints);
	}

	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "## Surface File Generated by Brian Chen and modified Trollbase Lib   ##\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "\n");
	///output original info
	sprintf(tempString, "#Original PDB file: none (Surface computations)\n");
	fprintf( currentOutputFile, tempString);

	///generate the Geometry Header
	fprintf( currentOutputFile, "########################################\n");
	fprintf( currentOutputFile, "###  ##   ##  ## ### #   #   #   ## # ##\n");
	fprintf( currentOutputFile, "## #### ### ## #  #  # #### ## ## # # ##\n");
	fprintf( currentOutputFile, "## #  #  ## ## # # # #  ### ##   ### ###\n");
	fprintf( currentOutputFile, "## ## # ### ## # ### # #### ## # ### ###\n");
	fprintf( currentOutputFile, "###  ##   ##  ## ### #   ## ## ## ## ###\n");
	fprintf( currentOutputFile, "########################################\n");
	sprintf(tempString, "GEOMETRY: %i\n", numVerts);
	fprintf( currentOutputFile, tempString);

	///generate the Geometry
	int * map = NULL;
	int * rmap = NULL;
	///if points is null, jsut output the geometry
	if(points==NULL){
		for(i = 0; i<numVerts; i++){
			sprintf(tempString, "%f %f %f %f %f %f\n", 
				surf->surfacePoints[3*i+0], surf->surfacePoints[3*i+1], surf->surfacePoints[3*i+2], 
				surf->surfaceNormals[3*i+0], surf->surfaceNormals[3*i+1], surf->surfaceNormals[3*i+2] );
			fprintf( currentOutputFile, tempString);
		}
	}
	//if points is not null, prepare forward and backward maps so that we cna quickly process the triangles.
	//then generate only the points used by this pocket.
	else{
		map = new int[size_set(points)];
		rmap = new int[surf->numPoints];
		for(i = 0; i<surf->numPoints; i++){ rmap[i] = -1; }
		for(i = 0; i<size_set(points); i++){
			map[i] = points[i];
			rmap[points[i]] = i;
		}
		for(i = 0; i<numVerts; i++){
			int tempVal = points[i];
			sprintf(tempString, "%f %f %f %f %f %f\n", 
				surf->surfacePoints[3*tempVal+0], surf->surfacePoints[3*tempVal+1], surf->surfacePoints[3*tempVal+2], 
				surf->surfaceNormals[3*tempVal+0], surf->surfaceNormals[3*tempVal+1], surf->surfaceNormals[3*tempVal+2] );
			fprintf( currentOutputFile, tempString);
		}
	}

	///generate the Topology Header
	fprintf( currentOutputFile, "########################################\n");
	fprintf( currentOutputFile, "##   ##  ##   ###  ## ####  ###  ## # ##\n");
	fprintf( currentOutputFile, "### ## ## # ## # ## # ### ## # #### # ##\n");
	fprintf( currentOutputFile, "### ## ## #   ## ## # ### ## # #  ## ###\n");
	fprintf( currentOutputFile, "### ## ## # #### ## # ### ## # ## ## ###\n");
	fprintf( currentOutputFile, "### ###  ## #####  ##   ##  ###  ### ###\n");
	fprintf( currentOutputFile, "########################################\n");
	
	///count the number of triangles:
	int numTriangles = 0;
	///If points is null, generate all triangles
	if(points==NULL){
		numTriangles = surf->numTriangles;
	}
	///If points is not null, count the triangles list for the triangles which are
	///entirely contained within points, and only output those.
	else{
		for(i = 0; i<surf->numTriangles; i++){
			bool test1 = contains_set(points, surf->triangles[3*i+0]);
			bool test2 = contains_set(points, surf->triangles[3*i+1]);
			bool test3 = contains_set(points, surf->triangles[3*i+2]);
			if(test1 && test2 && test3){
				numTriangles++;
			}
		}		
	}	
	///output the number of triangles.	
	sprintf(tempString, "TOPOLOGY: %i\n", numTriangles);
	fprintf( currentOutputFile, tempString);

	//now actually output the triangles that we need.Use the reverse map.
	int point1Idx, point2Idx, point3Idx;
	///if points is null, just output everything.
	if(points==NULL){
		for(i = 0; i<numTriangles; i++){
			///get the triangle indices from surfTriangles
			point1Idx = surf->triangles[3*i+0];
			point2Idx = surf->triangles[3*i+1];
			point3Idx = surf->triangles[3*i+2];
			sprintf(tempString, "%i %i %i\n", point1Idx, point2Idx, point3Idx);
			fprintf( currentOutputFile, tempString);
		}
	}
	//If poitns is non-null, output only the ones we care about.  Use the reverse map.
	else{
		for(i = 0; i<surf->numTriangles; i++){
			bool test1 = contains_set(points, surf->triangles[3*i+0]);
			bool test2 = contains_set(points, surf->triangles[3*i+1]);
			bool test3 = contains_set(points, surf->triangles[3*i+2]);
			if(test1 && test2 && test3){
				sprintf(tempString, "%i %i %i\n", rmap[surf->triangles[3*i+0]], rmap[surf->triangles[3*i+1]], rmap[surf->triangles[3*i+2]]);
				fprintf( currentOutputFile, tempString);
			}
		}
	}

	if(surf->colors != NULL){	
		///generate the Colors Header
		fprintf( currentOutputFile, "\n");
		fprintf( currentOutputFile, "########################################\n");
		fprintf( currentOutputFile, "###  ##  ## ####  ##  ###   ############\n");
		fprintf( currentOutputFile, "## ### ## # ### ## # # # ###############\n");
		fprintf( currentOutputFile, "## ### ## # ### ## #  ###  #############\n");
		fprintf( currentOutputFile, "## ### ## # ### ## # # #### ############\n");
		fprintf( currentOutputFile, "###  ##  ##   ##  ## # #   #############\n");
		fprintf( currentOutputFile, "########################################\n");
		fprintf( currentOutputFile, "COLORS: \n");
		for(i = 0; i<surf->numPoints; i++){
			sprintf(tempString, "%f %f %f\n", surf->colors[3*i+0], surf->colors[3*i+1], surf->colors[3*i+2] );
			fprintf( currentOutputFile, tempString);
		}
	}



	///now output the footer.
	fprintf( currentOutputFile, "\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "## Surface File Generated by Brian Chen and modified Trollbase Lib   ##\n");
	fprintf( currentOutputFile, "#######################################################################\n");
	fprintf( currentOutputFile, "#######################################################################\n");

	///done outputting.  Close file.
	fclose( currentOutputFile );
	
	printf("Output SURF File [%s] complete.\n", outputFileName);

	delete[](tempString);
	
//	no need to delete a closed file pointer
//	delete(currentOutputFile);
	if(map != NULL){
		delete[](map);
	}
	if(rmap != NULL){
		delete[](rmap);
	}
}








////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VOLUMETRIC ANALYSIS OF SURFACE PROPERTIES (VASP)                                                           //
// BRIAN Y CHEN, HONIG LAB, 2010                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////






















