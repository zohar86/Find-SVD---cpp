#pragma once

#ifndef SVD_H
#define SVD_H

using namespace std;

class SVD{

	size_t row;
	size_t col;
	string fileInput;
	double eps = 0.000001;
	int maxIter;

public: //array
	double* A = NULL; //A matrix 
	double* ATA = NULL; // ATA matrix
	double* segma = NULL; //output array with eigenvaluesMatrix
	double* segmaTemp = NULL; //output array with eigenvalues only
	double* VT = NULL; //output array VT matrix
	double* U = NULL; // output array U matrix

public:
	/*constructor*/
	SVD(double *A, double* mat, size_t sizeATA, size_t sizeU, size_t sizeSegma, size_t row, size_t col, int maxIter, string fileInput);
	/*Destructor */
	~SVD();
	/*Calculation of matrix U*/
	void find_U_Matrix();
	/**/
	void mat_identity(double v[]);
	/*SVD calculation with the help of a Jacobian algorithm*/
	void find_svd();

	
private:

};

#endif

/**********Function*************/
bool compare(const pair<int, int>& i, const pair<int, int>& j);
void find_size(string inputfile, size_t &row, size_t &col);
void read_csv(double* mat, string inputfile);
void write_csv(double* mat, string inputfile, size_t ind, int type);
void printMartix(double VT[], double segma[], double U[], size_t row, size_t col);
/**********************************************************/

