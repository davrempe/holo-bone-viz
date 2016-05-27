//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//Stand in Math Library by Brian Chen/////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/*********************************************************

	Contents:

	3x3 matrix Inversion
	3x3 matrix X 3d vector transformation
	3 element sort
	3 dimensional Cross Product
	3x3 matrix determinant
	3x3 times 3x3 matrix multiplication
	3d vector normalization
	3d vector magnitude
	3 3d vertex centroid

	3 dimensional scalar multiplication
	3 dimensional dot product

*/

#include "pch.h"
#include "mathlib.h"




double * transformVector(double * matrix, double * vector)
{
	//Matrix Arrangement			Vector Arrangement
	//[0, 3. 6]						[0]
	//[1, 4. 7]				X		[1]
	//[2, 5. 8]						[2]

	double * result = new double[3];

	////apply matrix multiplication to the vector
		result[0] = (matrix[0]*vector[0]) + (matrix[3]*vector[1]) + (matrix[6]*vector[2]);
		result[1] = (matrix[1]*vector[0]) + (matrix[4]*vector[1]) + (matrix[7]*vector[2]);
		result[2] = (matrix[2]*vector[0]) + (matrix[5]*vector[1]) + (matrix[8]*vector[2]);

	return result;

}



////////////////////////////
///3 element sorter
////FIXME: make me more efficient.
int * sort3(int first, int second, int third)
{

	int * result = new int[3];
	int state = 0;

	/////find the first element
		if		( ( first >= second) && (first >= third) )	{result[0] = first; state = 0;}
		else if ( ( second >= first) && (second >= third) )	{result[0] = second; state = 1;}
		else if ( ( third >= first) && (third >= second) )	{result[0] = third; state = 2;}


	/////find the second element, and set the third element
		switch (state){
		case 0 :	if (second <= third)	{result[1] = third; result[2] = second; break;}
					else					{result[1] = second; result[2] = third; break;}
		
		case 1 :	if (first <= third)		{result[1] = third; result[2] = first; break;}
					else					{result[1] = first; result[2] = third; break;}

		case 2 :	if (second <= first)	{result[1] = first; result[2] = second; break;}
					else					{result[1] = second; result[2] = first; break;}
		}


	return result;

}

////simple determinant for 3x3 matrices
double determinant3x3(double * matrix){

	double det = 0;

	////calculate determinant
		det = matrix[0]*( (matrix[4]*matrix[8]) - (matrix[7]*matrix[5]) );
		det = det - matrix[3]*( (matrix[1]*matrix[8]) - (matrix[7]*matrix[2]) );
		det = det + matrix[6]*( (matrix[1]*matrix[5]) - (matrix[4]*matrix[2]) );


	return det;

}


////simple determinant for 4x4 matrices
double determinant4x4(double * m)
{
	double det = 0;

	// 0  4  8  12
	// 1  5  9  13
	// 2  6  10 14
	// 3  7  11 15

	double * d0 = new double[9];
	d0[0]=m[5]; d0[1]=m[6]; d0[2]=m[7]; d0[3]=m[9]; d0[4]=m[10]; d0[5]=m[11]; d0[6]=m[13]; d0[7]=m[14]; d0[8]=m[15];
	double * d1 = new double[9];
	d1[0]=m[1]; d1[1]=m[2]; d1[2]=m[3]; d1[3]=m[9]; d1[4]=m[10]; d1[5]=m[11]; d1[6]=m[13]; d1[7]=m[14]; d1[8]=m[15];
	double * d2 = new double[9];
	d2[0]=m[1]; d2[1]=m[2]; d2[2]=m[3]; d2[3]=m[5]; d2[4]=m[6]; d2[5]=m[7]; d2[6]=m[13]; d2[7]=m[14]; d2[8]=m[15];
	double * d3 = new double[9];
	d3[0]=m[1]; d3[1]=m[2]; d3[2]=m[3]; d3[3]=m[5]; d3[4]=m[6]; d3[5]=m[7]; d3[6]=m[9]; d3[7]=m[10]; d3[8]=m[11];

	////calculate determinant
	if(m[ 0] != 0){ det += m[ 0] * determinant3x3(d0); }
	if(m[ 4] != 0){ det -= m[ 4] * determinant3x3(d1); }
	if(m[ 8] != 0){ det += m[ 8] * determinant3x3(d2); }
	if(m[12] != 0){ det -= m[12] * determinant3x3(d3); }

	delete[](d0);
	delete[](d1);
	delete[](d2);
	delete[](d3);

	return det;
}

////simple determinant for 5x5 matrices
double determinant5x5(double * m)
{
	double det = 0;

	// 0  5  10 15 20
	// 1  6  11 16 21
	// 2  7  12 17 22
	// 3  8  13 18 23
	// 4  9  14 19 24

	double * d0 = new double[16];
	d0[ 0]=m[ 6]; d0[ 1]=m[ 7]; d0[ 2]=m[ 8]; d0[ 3]=m[ 9];	d0[ 4]=m[11]; d0[ 5]=m[12]; d0[ 6]=m[13]; d0[ 7]=m[14]; 
	d0[ 8]=m[16]; d0[ 9]=m[17]; d0[10]=m[18]; d0[11]=m[19];	d0[12]=m[21]; d0[13]=m[22]; d0[14]=m[23]; d0[15]=m[24]; 
	double * d1 = new double[16];
	d1[ 0]=m[ 1]; d1[ 1]=m[ 2]; d1[ 2]=m[ 3]; d1[ 3]=m[ 4];	d1[ 4]=m[11]; d1[ 5]=m[12]; d1[ 6]=m[13]; d1[ 7]=m[14]; 
	d1[ 8]=m[16]; d1[ 9]=m[17]; d1[10]=m[18]; d1[11]=m[19];	d1[12]=m[21]; d1[13]=m[22]; d1[14]=m[23]; d1[15]=m[24]; 
	double * d2 = new double[16];
	d2[ 0]=m[ 1]; d2[ 1]=m[ 2]; d2[ 2]=m[ 3]; d2[ 3]=m[ 4];	d2[ 4]=m[ 6]; d2[ 5]=m[ 7]; d2[ 6]=m[ 8]; d2[ 7]=m[ 9]; 
	d2[ 8]=m[16]; d2[ 9]=m[17]; d2[10]=m[18]; d2[11]=m[19];	d2[12]=m[21]; d2[13]=m[22]; d2[14]=m[23]; d2[15]=m[24]; 
	double * d3 = new double[16];
	d3[ 0]=m[ 1]; d3[ 1]=m[ 2]; d3[ 2]=m[ 3]; d3[ 3]=m[ 4];	d3[ 4]=m[ 6]; d3[ 5]=m[ 7]; d3[ 6]=m[ 8]; d3[ 7]=m[ 9];
	d3[ 8]=m[11]; d3[ 9]=m[12]; d3[10]=m[13]; d3[11]=m[14];	d3[12]=m[21]; d3[13]=m[22]; d3[14]=m[23]; d3[15]=m[24]; 
	double * d4 = new double[16];
	d4[ 0]=m[ 1]; d4[ 1]=m[ 2]; d4[ 2]=m[ 3]; d4[ 3]=m[ 4];	d4[ 4]=m[ 6]; d4[ 5]=m[ 7]; d4[ 6]=m[ 8]; d4[ 7]=m[ 9]; 
	d4[ 8]=m[11]; d4[ 9]=m[12]; d4[10]=m[13]; d4[11]=m[14];	d4[12]=m[16]; d4[13]=m[17]; d4[14]=m[18]; d4[15]=m[19]; 

	////calculate determinant
	if(m[ 0] != 0){ det += m[ 0] * determinant4x4(d0); }
	if(m[ 5] != 0){ det -= m[ 5] * determinant4x4(d1); }
	if(m[10] != 0){ det += m[10] * determinant4x4(d2); }
	if(m[15] != 0){ det -= m[15] * determinant4x4(d3); }
	if(m[20] != 0){ det += m[20] * determinant4x4(d4); }

	delete[](d0);
	delete[](d1);
	delete[](d2);
	delete[](d3);
	delete[](d4);

	return det;
}



double * transpose(double * mat)
{
	//Matrix Arrangement
	//[0, 3. 6]
	//[1, 4. 7]
	//[2, 5. 8]

	double * result = new double[9];

	result[0] = mat[0];
	result[1] = mat[3];
	result[2] = mat[6];
	result[3] = mat[1];
	result[4] = mat[4];
	result[5] = mat[7];
	result[6] = mat[2];
	result[7] = mat[5];
	result[8] = mat[8];

	return result;

}

double * transpose4x4(double * mat)
{
	//Matrix Arrangement
	//[0, 4. 8  12]
	//[1, 5. 9  13]
	//[2, 6. 10 14]
	//[3, 7, 11,15] 

	double * result = new double[16];

	result[0]  = mat[0];
	result[1]  = mat[4];
	result[2]  = mat[8];
	result[3]  = mat[12];
	result[4]  = mat[1];
	result[5]  = mat[5];
	result[6]  = mat[9];
	result[7]  = mat[13];
	result[8]  = mat[2];
	result[9]  = mat[6];
	result[10] = mat[10];
	result[11] = mat[14];
	result[12] = mat[3];
	result[13] = mat[7];
	result[14] = mat[11];
	result[15] = mat[15];

	return result;

}

double * invertMatrix(double * mat)
{

	//Matrix Arrangement
	//[0, 3. 6]
	//[1, 4. 7]
	//[2, 5. 8]

	double * matrix = new double[9];
	double det;
	double * adjoint = new double[9];
	double * inverse = new double[9];

	////copy out the values
		matrix[0] = mat[0];	
		matrix[1] = mat[1];	
		matrix[2] = mat[2];	
		matrix[3] = mat[3];	
		matrix[4] = mat[4];	
		matrix[5] = mat[5];	
		matrix[6] = mat[6];	
		matrix[7] = mat[7];	
		matrix[8] = mat[8];	

	///calculate determinant
		det = determinant3x3(matrix);

	///calculate adjoint - it's the Transpose of the matrix of component determinants.
		adjoint[0] =		( (matrix[4]*matrix[8]) - (matrix[7]*matrix[5]) );
		adjoint[3] =	-	( (matrix[3]*matrix[8]) - (matrix[6]*matrix[5]) );
		adjoint[6] =		( (matrix[3]*matrix[7]) - (matrix[6]*matrix[4]) );
		adjoint[1] =	-	( (matrix[1]*matrix[8]) - (matrix[7]*matrix[2]) );
		adjoint[4] =		( (matrix[0]*matrix[8]) - (matrix[6]*matrix[2]) );
		adjoint[7] =	-	( (matrix[0]*matrix[7]) - (matrix[6]*matrix[1]) );
		adjoint[2] =		( (matrix[1]*matrix[5]) - (matrix[4]*matrix[2]) );
		adjoint[5] =	-	( (matrix[0]*matrix[5]) - (matrix[3]*matrix[2]) );
		adjoint[8] =		( (matrix[0]*matrix[4]) - (matrix[3]*matrix[1]) );

	///calculate inverse
		inverse[0] = adjoint[0]/det;
		inverse[1] = adjoint[1]/det;
		inverse[2] = adjoint[2]/det;
		inverse[3] = adjoint[3]/det;
		inverse[4] = adjoint[4]/det;
		inverse[5] = adjoint[5]/det;
		inverse[6] = adjoint[6]/det;
		inverse[7] = adjoint[7]/det;
		inverse[8] = adjoint[8]/det;

	delete[](matrix);
	delete[](adjoint);

	return inverse;
}



////////////////////
//Returns normalized vector in input direction
double * normalizeVector(double * vector)
{
	double * result = new double[3];

	////calculate the magnitude
		double mag = sqrt(	(vector[0]*vector[0]) + 
							(vector[1]*vector[1]) + 
							(vector[2]*vector[2]) );

	////set the normalized vector
		result[0] = vector[0]/mag;
		result[1] = vector[1]/mag;
		result[2] = vector[2]/mag;

	return result;

}

double * normalizeVector(double v0, double v1, double v2)
{
	double * result = new double[3];
	double mag = sqrt( (v0*v0) + (v1*v1) + (v2*v2) );
	result[0] = v0/mag;
	result[1] = v1/mag;
	result[2] = v2/mag;

	return result;	
}

void normalizeVector(double * vector, double * output)
{
	double v0 = vector[0];
	double v1 = vector[1];
	double v2 = vector[2];

	double mag = sqrt( (v0*v0) + (v1*v1) + (v2*v2) );

	output[0] = v0/mag;
	output[1] = v1/mag;
	output[2] = v2/mag;

}



///finds the matrix magnitude
double vectorSize(double x, double y, double z)
{
	return sqrt( (x*x) + (y*y) + (z*z) );
}

double vectorSize2d(double x, double y)
{
	return sqrt( (x*x) + (y*y) );

}

///finds the matrix magnitude
double vectorSize(double * A, double * B)
{
	double x = B[0]-A[0];
	double y = B[1]-A[1];
	double z = B[2]-A[2];

	return sqrt( (x*x) + (y*y) + (z*z) );
}



double * crossProd(double * first, double * second)
{
	///Fix size of array to 3 doubles
	double * result = new double[3];	

	/////calculate cross product
		result[0] = (first[1] * second[2]) - (first[2] * second[1]);
		result[1] = -( (first[0] * second[2]) - (first[2] * second[0]) );
		result[2] = (first[0] * second[1]) - (first[1] * second[0]);

	return result;

}

double * getCentroid(double * vertices, int size)
{
	double * result = new double[3];
	double tempSum = 0.0;
	int i = 0;
	
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+0];}
	result[0] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+1];}
	result[1] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+2];}
	result[2] = tempSum/size;

	return result;
}



///////////////////////////
///3x3 matrix times 3x3 matrix multiplication
double * matrixMult(double * first, double * second)
{
	double * result = new double[9];
	
		result[0] = (first[0]*second[0]) + (first[3]*second[1]) + (first[6]*second[2]);
		result[1] = (first[1]*second[0]) + (first[4]*second[1]) + (first[7]*second[2]);
		result[2] = (first[2]*second[0]) + (first[5]*second[1]) + (first[8]*second[2]);
		result[3] = (first[0]*second[3]) + (first[3]*second[4]) + (first[6]*second[5]);
		result[4] = (first[1]*second[3]) + (first[4]*second[4]) + (first[7]*second[5]);
		result[5] = (first[2]*second[3]) + (first[5]*second[4]) + (first[8]*second[5]);
		result[6] = (first[0]*second[6]) + (first[3]*second[7]) + (first[6]*second[8]);
		result[7] = (first[1]*second[6]) + (first[4]*second[7]) + (first[7]*second[8]);
		result[8] = (first[2]*second[6]) + (first[5]*second[7]) + (first[8]*second[8]);

	return result;
}


double * scalarMultiplication(double * vector, double scale)
{
	double * result = new double[3];
	
	/////calculate the scaled vector
		result[0] = vector[0]*scale;
		result[1] = vector[1]*scale;
		result[2] = vector[2]*scale;

	return result;

}



double dotProd(double * first, double * second)
{
	double result = 0;

	/////calculate the dot product
		result =	(first[0]*second[0]) +
					(first[1]*second[1]) +
					(first[2]*second[2]);

	return result;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting values out of the file when the position is fixed//////////////
/////Unlike strcpy which copies everything; I want to just copy specifics////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
char * stringCopy(char* str, int start, int end)
{
	int size = end - start + 1;
	char * result = new char[size + 1];

	for(int i = 0; i<size; i++){
		result[i] = str[start+i];	
	}
	result[size] = '\0';

	return result;

}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting the first character out of a string
/////////////////////////////////////////////////////////////////////////////////////////////////
char * firstChar(char * str)
{
	char * result = new char[2];
	if (str != NULL){
		result[0] = str[0];
		result[1] = '\0';
	}
	return result;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting a bool fast
/////////////////////////////////////////////////////////////////////////////////////////////////
bool makeBool(int * i)
{
	if(i[0] == 1){
		return true;
	}
	else{
		return false;
	}


}



///returns a random double between 0 and the bound
///MAKE SURE TO SEED THE RANDOM NUMBER GENERATOR FIRST
double getRand(double bound)
{
	if(bound == 0.0){
		return 0.0;
	}
	double stuff = (double) rand();
	double stuff2 = (double)(stuff/RAND_MAX);
	double stuff3 = stuff2 * bound;
	return stuff3;

}


///Given side lengths a, b, c, returns the angle opposite A, in degrees
double SSStriangleAngle(double a, double b, double c)
{
	if( (a+b-c > 0) && (a+b-c < .0001 )){			//process extreme triangles to prevent strange vals
		return 0.0;
	}
	if( (c+a-b > 0) && (c+a-b < .0001 )){			//process extreme triangles to prevent strange vals
		return 0.0;
	}
	if( (c+b-a > 0) && (c+b-a < .0001 )){			//process extreme triangles to prevent strange vals
		return 180.0;
	}

	double param, result;
	param = ( ( (a*a)-(b*b)-(c*c) ) / (2*b*c) );
	result = 180 - (acos(param) * 180 / PI);
	return result;
}


////prints out a matrix for debugging
void printMatrix(double * mat, int size)
{
	switch(size){
	case 9:
		printf("%06f %06f %06f \n", mat[0], mat[3], mat[6]);
		printf("%06f %06f %06f \n", mat[1], mat[4], mat[7]);
		printf("%06f %06f %06f \n", mat[2], mat[5], mat[8]);
		break;
	
	case 16:
		printf("%06f %06f %06f %06f\n", mat[0], mat[4], mat[8], mat[12]);
		printf("%06f %06f %06f %06f\n", mat[1], mat[5], mat[9], mat[13]);
		printf("%06f %06f %06f %06f\n", mat[2], mat[6], mat[10], mat[14]);
		printf("%06f %06f %06f %06f\n\n", mat[3], mat[7], mat[11], mat[15]);
		break;
	
	}
}

////transforms a 3-vector with a 4matrix, by padding a 1 for scale.
double * transformVector3x4(double * mat, double * vector)
{
	double * result = new double[3];

	result[0] = mat[0]*vector[0] + mat[4]*vector[1] + mat[8]*vector[2] + mat[12];
	result[1] = mat[1]*vector[0] + mat[5]*vector[1] + mat[9]*vector[2] + mat[13];
	result[2] = mat[2]*vector[0] + mat[6]*vector[1] + mat[10]*vector[2] + mat[14];

	return result;
}


///////////////////////////
///4x4 matrix times 4x4 matrix multiplication
double * matrixMult4x4(double * first, double * second)
{
	double * result = new double[16];

	result[0]  = first[ 0]*second[ 0] + first[ 4]*second[ 1] + first[ 8]*second[ 2] + first[12]*second[ 3];
	result[1]  = first[ 1]*second[ 0] + first[ 5]*second[ 1] + first[ 9]*second[ 2] + first[13]*second[ 3];
	result[2]  = first[ 2]*second[ 0] + first[ 6]*second[ 1] + first[10]*second[ 2] + first[14]*second[ 3];
	result[3]  = first[ 3]*second[ 0] + first[ 7]*second[ 1] + first[11]*second[ 2] + first[15]*second[ 3];

	result[4]  = first[ 0]*second[ 4] + first[ 4]*second[ 5] + first[ 8]*second[ 6] + first[12]*second[ 7];
	result[5]  = first[ 1]*second[ 4] + first[ 5]*second[ 5] + first[ 9]*second[ 6] + first[13]*second[ 7];
	result[6]  = first[ 2]*second[ 4] + first[ 6]*second[ 5] + first[10]*second[ 6] + first[14]*second[ 7];
	result[7]  = first[ 3]*second[ 4] + first[ 7]*second[ 5] + first[11]*second[ 6] + first[15]*second[ 7];

	result[8]  = first[ 0]*second[ 8] + first[ 4]*second[ 9] + first[ 8]*second[10] + first[12]*second[11];
	result[9]  = first[ 1]*second[ 8] + first[ 5]*second[ 9] + first[ 9]*second[10] + first[13]*second[11];
	result[10] = first[ 2]*second[ 8] + first[ 6]*second[ 9] + first[10]*second[10] + first[14]*second[11];
	result[11] = first[ 3]*second[ 8] + first[ 7]*second[ 9] + first[11]*second[10] + first[15]*second[11];

	result[12] = first[ 0]*second[12] + first[ 4]*second[13] + first[ 8]*second[14] + first[12]*second[15];
	result[13] = first[ 1]*second[12] + first[ 5]*second[13] + first[ 9]*second[14] + first[13]*second[15];
	result[14] = first[ 2]*second[12] + first[ 6]*second[13] + first[10]*second[14] + first[14]*second[15];
	result[15] = first[ 3]*second[12] + first[ 7]*second[13] + first[11]*second[14] + first[15]*second[15];

	return result;
}

bool testMatrix(double * matrix)
{
	double * matrixa = new double[16];
	matrixa[0] = matrix[0];
	matrixa[1] = matrix[1];
	matrixa[2] = matrix[2];
	matrixa[3] = matrix[3];
	matrixa[4] = matrix[4];
	matrixa[5] = matrix[5];
	matrixa[6] = matrix[6];
	matrixa[7] = matrix[7];
	matrixa[8] = matrix[8];
	matrixa[9] = matrix[9];
	matrixa[10] = matrix[10];
	matrixa[11] = matrix[11];
	matrixa[12] = 0;
	matrixa[13] = 0;
	matrixa[14] = 0;
	matrixa[15] = matrix[15];
//	printMatrix(matrixa, 16);
	double * matrix1 = transpose4x4(matrixa);
//	printMatrix(matrix1, 16);

	double * matrix2 = matrixMult4x4(matrix1, matrixa);
//	printMatrix(matrix2, 16);
//
	double * matrixI = new double[16];
	matrixI[0]  = 1 - matrix2[0];
	matrixI[1]  = 0 - matrix2[1];
	matrixI[2]  = 0 - matrix2[2];
	matrixI[3]  = 0 - matrix2[3];
	matrixI[4]  = 0 - matrix2[4];
	matrixI[5]  = 1 - matrix2[5];
	matrixI[6]  = 0 - matrix2[6];
	matrixI[7]  = 0 - matrix2[7];
	matrixI[8]  = 0 - matrix2[8];
	matrixI[9]  = 0 - matrix2[9];
	matrixI[10] = 1 - matrix2[10];
	matrixI[11] = 0 - matrix2[11];
	matrixI[12] = 0 - matrix2[12];
	matrixI[13] = 0 - matrix2[13];
	matrixI[14] = 0 - matrix2[14];
	matrixI[15] = 1 - matrix2[15];

//	printMatrix(matrixI, 16);

	int i = 0;
	for(i = 0; i<16; i++){
		if(matrixI[i] > 0.01 || matrixI[i] < -0.01 ){
			printf("WHOA\n");

			delete[](matrixa);
			delete[](matrix1);
			delete[](matrix2);
			delete[](matrixI);
			return false;
		}
	}

	delete[](matrixa);
	delete[](matrix1);
	delete[](matrix2);
	delete[](matrixI);
	return true;

}



bool testMatrix3x3(double * matrix)
{
	double * matrixa = new double[9];
	matrixa[0] = matrix[0];
	matrixa[1] = matrix[1];
	matrixa[2] = matrix[2];
	matrixa[3] = matrix[3];
	matrixa[4] = matrix[4];
	matrixa[5] = matrix[5];
	matrixa[6] = matrix[6];
	matrixa[7] = matrix[7];
	matrixa[8] = matrix[8];

	printMatrix(matrixa, 9);
	double * matrix1 = transpose(matrixa);
	printMatrix(matrix1, 9);

	double * matrix2 = matrixMult(matrix1, matrixa);
	printMatrix(matrix2, 9);

	double * matrixI = new double[9];
	matrixI[0]  = 1 - matrix2[0];
	matrixI[1]  = 0 - matrix2[1];
	matrixI[2]  = 0 - matrix2[2];
	matrixI[3]  = 0 - matrix2[3];
	matrixI[4]  = 1 - matrix2[4];
	matrixI[5]  = 0 - matrix2[5];
	matrixI[6]  = 0 - matrix2[6];
	matrixI[7]  = 0 - matrix2[7];
	matrixI[8]  = 1 - matrix2[8];

	printMatrix(matrixI, 9);

	int i = 0;
	for(i = 0; i<9; i++){
		if(matrixI[i] > 0.01 || matrixI[i] < -0.01 ){
			printf("WHOA\n");

			delete[](matrixa);
			delete[](matrix1);
			delete[](matrix2);
			delete[](matrixI);
			return false;
		}
	}

	delete[](matrixa);
	delete[](matrix1);
	delete[](matrix2);
	delete[](matrixI);

	return true;

}

///From Wikipedia:
///The area of a parallelogram can be calculated using vectors. 
///Let vectors AB and AC point respectively from A to B and from A to C. 
///The area of parallelogram ABDC is then |AB × AC|, which is the magnitude 
///of the cross product of vectors AB and AC. |AB × AC| is equal to |h × AC|, 
///where h represents the altitude h as a vector.
double triangleArea(double * t1, double * t2, double * t3)
{
	double t10 = t1[0];	double t11 = t1[1];	double t12 = t1[2];
	double t20 = t2[0];	double t21 = t2[1];	double t22 = t2[2];
	double t30 = t3[0];	double t31 = t3[1];	double t32 = t3[2];
	
	double cross0 = ( (t21-t11)*(t32-t12) ) - ( (t22-t12)*(t31-t11) );
	double cross1 = -( ((t20-t10)*(t32-t12)) - ((t22-t12)*(t30-t10)) );
	double cross2 = ( (t20-t10)*(t31-t11) ) - ( (t21-t11)*(t30-t10) );

	double magABDC = vectorSize(cross0, cross1, cross2);
	double result = magABDC/2;

	return result;
}


double triangleArea(double * points, int a, int b, int c)
{
	double t10 = points[3*a+0];	double t11 = points[3*a+1];	double t12 = points[3*a+2];
	double t20 = points[3*b+0];	double t21 = points[3*b+1];	double t22 = points[3*b+2];
	double t30 = points[3*c+0];	double t31 = points[3*c+1];	double t32 = points[3*c+2];
	double cross0 = ( (t21-t11)*(t32-t12) ) - ( (t22-t12)*(t31-t11) );
	double cross1 = -( ((t20-t10)*(t32-t12)) - ((t22-t12)*(t30-t10)) );
	double cross2 = ( (t20-t10)*(t31-t11) ) - ( (t21-t11)*(t30-t10) );

	double magABDC = vectorSize(cross0, cross1, cross2);
	double result = magABDC/2;

	return result;
}


double triangleArea(double * triangle)
{
	double t10 = triangle[0];	double t11 = triangle[1];	double t12 = triangle[2];
	double t20 = triangle[3];	double t21 = triangle[4];	double t22 = triangle[5];
	double t30 = triangle[6];	double t31 = triangle[7];	double t32 = triangle[8];
	double cross0 = ( (t21-t11)*(t32-t12) ) - ( (t22-t12)*(t31-t11) );
	double cross1 = -( ((t20-t10)*(t32-t12)) - ((t22-t12)*(t30-t10)) );
	double cross2 = ( (t20-t10)*(t31-t11) ) - ( (t21-t11)*(t30-t10) );

	double magABDC = vectorSize(cross0, cross1, cross2);
	double result = magABDC/2;

	return result;
}



double sphereicalTriangleArea(double * pt, double * t )
{
	//t is organized as {0,1,2}, {3,4,5}, {6,7,8}
	
	///get the R vectors
	double * R1 = new double[3];   R1[0] = t[0]-pt[0];   R1[1] = t[1]-pt[1];   R1[2] = t[2]-pt[2]; 
	double * R2 = new double[3];   R2[0] = t[3]-pt[0];   R2[1] = t[4]-pt[1];   R2[2] = t[5]-pt[2]; 
	double * R3 = new double[3];   R3[0] = t[6]-pt[0];   R3[1] = t[7]-pt[1];   R3[2] = t[8]-pt[2]; 

	///get the r magnitudes;
	double r1 = vectorSize( R1[0], R1[1], R1[2] );
	double r2 = vectorSize( R2[0], R2[1], R2[2] );
	double r3 = vectorSize( R3[0], R3[1], R3[2] );
	
	double * cross = new double[3];
	CROSS(cross, R2, R3);
	
	//evaluate the expression from Oosterom and Strackee, except for arctan.
	double topVal = DOT(R1,cross);
	double botVal = (r1*r2*r3) + (DOT(R1,R2)*r3) + (DOT(R1,R3)*r2) + (DOT(R2,R3)*r1);
	
	double fracResult = topVal/botVal;
	
	if(fracResult != fracResult){ printf("WTF fracresult: %f\n", fracResult); }
	
	double tanResult  = atan(fracResult);
	
	double result = 2*tanResult;
	
	delete[](R1);
	delete[](R2);
	delete[](R3);
	delete[](cross);
		
	return result;
}


bool intervalOverlapTest(double a1, double a2, double b1, double b2)
{
	double bota = a1;	double topa = a2;
	if(a2<a1){ topa = a1; bota = a2; }
	double botb = b1;	double topb = b2;
	if(b2<b1){ topb = b1; botb = b2; }

	bool result;
	if( (topb < bota) || (botb > topa) ){ result = false; }
	else{ result = true; }
	
	return result;
}



///return 2 if all even, 1 if all odd, 0 of inconsistent, -1 if size = 0;
int parityTest(int * array, int size)
{
	if(size == 0){ return -1; }
	
	int numEven = 0;
	int numOdd = 0;
	int i = 0;
	for(i = 0; i<size; i++){
		if( array[i] % 2 == 0 ){
			numEven++;
		}
		else{
			numOdd++;
		}
		if(numEven>0 && numOdd>0){
			return 0;
		}
	}
	
	int result = -1;

	if( numEven>0 && numOdd==0 ){
		result = 2;
	}
	if( numEven==0 && numOdd>0 ){
		result = 1;
	}
	
	return result;
}




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
//
//OPTIMIZATION: implement this so that it runs the 4x4 matrices at a more optimal line than the first.
//              (not coded yet)
//

double computeTetrahedralVolume( double * p1, double * p2, double * p3, double * p4 )
{
	/////////////////////////////////
	double d12 = vectorSize( p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2] );
	double d13 = vectorSize( p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2] );
	double d14 = vectorSize( p1[0]-p4[0], p1[1]-p4[1], p1[2]-p4[2] );
	
	double d23 = vectorSize( p2[0]-p3[0], p2[1]-p3[1], p2[2]-p3[2] );
	double d24 = vectorSize( p2[0]-p4[0], p2[1]-p4[1], p2[2]-p4[2] );
	
	double d34 = vectorSize( p3[0]-p4[0], p3[1]-p4[1], p3[2]-p4[2] );
	/////////////////////////////////
	double s12 = d12*d12;
	double s13 = d13*d13;
	double s14 = d14*d14;

	double s23 = d23*d23;
	double s24 = d24*d24;

	double s34 = d34*d34;
	/////////////////////////////////
	double * m = new double[25];
	m[ 0] = 0;		m[ 5] = 1;		m[10] = 1;		m[15] = 1;		m[20] = 1;
	m[ 1] = 1;		m[ 6] = 0;		m[11] = s12;	m[16] = s13;	m[21] = s14;
	m[ 2] = 1;		m[ 7] = s12;	m[12] = 0;		m[17] = s23;	m[22] = s24;
	m[ 3] = 1;		m[ 8] = s13;	m[13] = s23;	m[18] = 0;		m[23] = s34;
	m[ 4] = 1;		m[ 9] = s14;	m[14] = s24;	m[19] = s34;	m[24] = 0;
	/////////////////////////////////
	double det = determinant5x5(m);

	delete[](m);

	double result = sqrt(det/288);
	
	if (result != result){
		printf("WARNING: ComputeTetrahedralVolume generated a NAN tetrahedron!\n");
		result = 0;
	}
	
	
	return result;
}



