#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>   // for file-reading
#include <sstream>   // for file-reading
#include <algorithm>
#include <vector>
#include <chrono>

#include "SVD.h"

using namespace std;

SVD::SVD(double *A, double* ATA, size_t sizeATA, size_t sizeU, size_t sizeSegma, size_t row, size_t col, int maxItr, string fileInput){
    this->row = row;
    this->col = col;
    this->maxIter = maxItr;
    this->fileInput = fileInput;

    this->A = A;
    this->ATA = ATA;

    segma = new double[sizeATA] {0};
    segmaTemp = new double [sizeSegma] {0};
    VT = new double [sizeATA] {0};
    U = new double [sizeU] {0};
}

SVD::~SVD() {
    delete[] segma;
    delete[] segmaTemp;
    delete[] VT;
    delete[] U;
}

void SVD::find_U_Matrix() { 
    double* tempA = new double[row * col] {0};
    size_t j = 0;
    size_t z = 0;
    size_t g = 0;
    size_t ind1 = 1, ind2 = col;
    size_t n = 0;
    size_t m = 0;
    size_t min = row;
    size_t max = row;

    if (row > col) {
        min = col;
    }

    cout << "---------------------------------------U-----------------------------------------------" << endl;
    for (size_t i = 0; i < min; i++) {
        if (segmaTemp[i] != 0) {
            double temp = 1 / sqrt(fabs(segmaTemp[i]));
            for(size_t j=0; j<row*col; j++){
                tempA[j] = temp * A[j];
                j++;
            }
            j = 0;
            while (n < max) {
                while (m < min) {
                    U[z] = U[z] + tempA[j] * VT[g];
                    //cout << "U:" << U[z] << "  VT: " << VT[g] << "  tempA: " << tempA[j] <<" g: "<<g<< endl;
                    j++;
                    g++;
                    m++;
                }
                g = ind2 - col;
                z = z + max;
                n++;
                m = 0;
            }
            z = ind1;
            g = ind2;
            ind2 = ind2 + row;
            ind1++;
            j = 0;
            m = 0;
            n = 0;
        }
        else {
            z = ind1;
            g = ind2;
            ind2 = ind2 + col;
            ind1++;
            j = 0;
            m = 0;
            n = 0;
        }
    }
    cout << "-------------------------------------End U-------------------------------------------" << endl;
    cout << endl << endl;
}

void SVD:: mat_identity(double v[]) {

    size_t i = 0;

    for(i; i < (col * col); ){
        v[i] = 1;
        i++;
        i = i + col;
    }
    return;
}
/*SVD finding*/
void SVD::find_svd() {

    size_t i, j, p1, p2, q1, q2, k = 0;
    double maxL_abs, maxL = -1;
    size_t b;
    size_t indI, indJ, indMax;
    double* arrI1 = new double [col] {0};
    double theta;
    int flag = 1;

    mat_identity(VT);

    cout << "The file is: " << fileInput << endl;
    cout << "The number of iterations: " << maxIter << endl;

    cout << endl << endl;
    cout << "------------------------------------VT Matrix------------------------------------------" << endl;

    for (size_t o=0; o < maxIter; o++) {
        if (flag) {
            i = 1;
            j = 2;
            b = col;
            maxL_abs = fabs(ATA[i]);
            maxL = ATA[i];

            for (i; i < col * col;) {
                while (i < b) {
                    if (fabs(ATA[i]) >= maxL_abs) {
                        maxL_abs = fabs(ATA[i]);
                        maxL = ATA[i];
                        indI = i / col;
                        indJ = i % col;
                        indMax = i;
                    }
                    i++;
                }
                i = i + j;
                j++;
                b = b + col;
            }
            if (fabs(maxL) <= eps) {
                flag = 0;
            }

            size_t ii = indI * col + indI;
            size_t jj = indJ * col + indJ;
            double c;
            double s;

            if (fabs(ATA[ii] - ATA[jj]) < eps) {
                theta = (M_PI_4);
            }
            else {
                double m = (ATA[jj] - ATA[ii]);
                double res = (2 * ATA[indMax]) / m;
                theta = atan(res) * 0.5;
            }

            c = cos(theta);
            s = sin(theta);

            size_t tempindI = indJ;
            size_t tempindJ = indI;
            size_t tempindMax = tempindI * col + tempindJ;

            double tempI = ATA[ii];//Value in diagonal ii
            double tempJ = ATA[jj];//Value in diagonal jj
            double tempMaxIJ = ATA[indMax];//upper max
            double tempMaxJI = ATA[tempindMax];//lower max

            double m = tempI - tempJ;
            if (fabs(m) < eps) {
                m = 0;
            }

            ATA[ii] = (c * c * tempI) - (2 * s * c * tempMaxIJ) + (s * s * tempJ);
            ATA[jj] = (s * s * tempI) + (2 * s * c * tempMaxIJ) + (c * c * tempJ);
            //ATA[indMax] = ATA[tempindMax] = (c * c - s * s) * (tempMaxIJ) + (s * c * (tempI - tempJ));
            ATA[indMax] = ATA[tempindMax] = (c * c - s * s) * (tempMaxIJ)+(s * c * (m));

            p1 = col * indI;
            p2 = indI;
            q1 = indJ;
            q2 = col * indJ;

            size_t x = 0, y = 0;

            size_t itr = col * (indI + 1);

            while (p1 < itr) {

                if (p1 != ii && p1 != indMax) {
                    arrI1[x++] = ATA[p1];
                    ATA[p1] = ATA[p2] = c * ATA[p1] - s * ATA[q1];
                }
                p1++;
                p2 = p2 + row;
                q1 = q1 + row;
            }

            p1 = col * indI;
            p2 = indI;
            q1 = indJ;
            q2 = col * indJ;

            x = 0;

            while (p1 < itr) {
                if (q2 != jj && q2 != tempindMax) {
                    ATA[q1] = ATA[q2] = s * arrI1[x] + c * ATA[q1];
                }
                p1++;
                x++;
                q1 = q1 + col;
                q2++;
            }

            size_t VindI = indI;
            size_t VindJ = indJ;
            size_t z = 0;

            for(size_t z=0; z<row; z++){
                double vip = VT[VindI];
                double viq = VT[VindJ];
                VT[VindI] = c * vip - s * viq;
                VT[VindJ] = c * viq + s * vip;
                VindI = VindI + col;
                VindJ = VindJ + col;
            }

            k = 0;
            b = col + 1;
            i = 0;
            size_t n = col;

            while (k < col * col) {
                while (k < n) {
                    if (fabs(ATA[k]) < eps && k != b && k != 0) {
                        if (ATA[k] != 0)
                            ATA[k] = 0;
                    }
                    k++;
                }
                i++;
                b = b + col + i;
                n = n + col;
            }
        }
    }

    int r = 0;
    size_t stop = col;

    vector<pair<int, double>> tempValsort;
    vector<int> index;
    vector<double> tempVecsort;
    int* arrindex = new int[col];
    int* arrindexSort = new int[col];

    i = 0;
    j = 0;
    for(i; i < col*col; ){
        segmaTemp[j] = ATA[i];
        arrindex[j] = j;
        tempValsort.push_back(make_pair(arrindex[j], (segmaTemp[j])));
        i++;
        i = i + col;
        j++;
    }

    sort(tempValsort.begin(), tempValsort.end(), compare);

    for (size_t i = 0; i < col; i++) {
        arrindexSort[i] = tempValsort[i].first;
    }

    size_t g = 0;
    for (size_t i = 0; i < col; i++) {
        j = arrindexSort[i];
        for(size_t g=0; g<col;g++){
            tempVecsort.push_back(VT[j]);
            j = j + col;
        }
    }

    j = 0;
    size_t z = 0, ind = 1;
    for (size_t i = 0; i < col; i++) {
        for (size_t g = 0; g < col; g++) {
            VT[z] = tempVecsort[z];
            z = z + col;
        }
        z = ind;
        ind++;
    }

    cout << "------------------------------------end VT-----------------------------------------" << endl;
    cout << endl << endl;

    cout << "-----------------------------------Segma----------------------------------------------" << endl;
    z = 0;
    i = 0;
    while (i < col * col) {
        segmaTemp[z] = tempValsort[z].second;
        segma[i] = tempValsort[z].second;
        z++;
        i++;
        i = i + col;
    }
    cout << "----------------------------------End Segma--------------------------------------------" << endl;
    cout << endl << endl;

    find_U_Matrix();//A, VT, segmaTemp, U

    delete[] arrindex;
    delete[] arrindexSort;
    delete[] arrI1;
}


/*****************************************************************/
/*Sorting eigenvalues*/
bool compare(const pair<int, int>& i, const pair<int, int>& j) {
    return i.second > j.second;
}
/*Find the size of the matrix rows and columns*/
void find_size(string inputfile, size_t &row, size_t &col) {

    string line;
    ifstream file(inputfile.c_str());// creates objects of ifstream class to make operations with the input file, i.e. reading from the input file
    size_t size = 0;
    double x;

    int flag = 1;
    
    if (file.is_open()) {
        while (getline(file, line)) {
            stringstream lineStream(line);
            string bit;
            if (flag) {
                while (getline(lineStream, bit, ',')) {
                    x = stof(bit);
                    col++;
                }
                flag = 0;
            }
            row++;
        }
    }
    cout << "rows: " << row << endl;
    cout << "cols: " << col << endl;
}
/*Convert text file to a matrix in C++*/
void read_csv(double* mat, string inputfile) {

    string line;
    ifstream file(inputfile.c_str());
    double x;
    size_t i = 0;
    while (getline(file, line)) {
        stringstream lineStream(line);
        string bit;
        while (getline(lineStream, bit, ',')) {
            x = stof(bit);
            mat[i++] = x;
        }
    }
    file.close();
}
/*write the result to file*/
void write_csv(double* mat, string inputfile, size_t ind, int type) {
    size_t i = 0, s = 0, b = 0;
    ofstream file(inputfile.c_str());
    string str = " ";

    if (type == 1) {
        b = ind * ind - ind;
        while (s < ind) {
            while (i < ind * ind) {
                str = to_string(mat[i]);

                if (i == b)
                    file << str << " ";
                else
                    file << str << ", ";
                i = i + ind;
            }
            s++;
            i = s;
            file << endl;
            b++;
        }
    }
    else if (type == 2) {
        b = ind;
        while (i < b) {
             str = to_string(sqrt(fabs(mat[i++])));
             if (i == b)
                 file << str << " ";

             else
                 file << str << ", ";
        }
    }
    else {
        b = ind;
        while (s < ind) {
            while (i < b) {
                str = to_string(mat[i]);

                if (i == b)
                    file << str << " ";
                else
                    file << str << ", ";
                i++;
            }
            s++;
            file << endl;
            b = b + ind;
        }
    }
    file.close();
}
/*print the results on screen*/
void printMartix(double VT[], double segma[], double U[], size_t row ,size_t col) {

    size_t j = 1, r = 0, stop = col;
    size_t index = 1;

    cout << "VT temp matrix: " << endl;
    for (int t = 0; t < col; t++) {
        while (r < col * col) {
            cout << VT[r] << " ";
            r = r + col;
        }
        r = j;
        j++;
        stop = stop + col;
        cout << endl;
    }
    cout << endl << endl;

    r = 0;
    stop = col;
    cout << "segma Matrix: " << endl;
    for (int t = 0; t < col; t++) {
        while (r < stop) {
            cout << sqrt(fabs(segma[r++])) << " ";
        }
        stop = stop + col;
        cout << endl;
    }
    r = 0;
    stop = col;
    cout << "segma Matrix: " << endl;

    for (int t = 0; t < col; t++) {
        cout << sqrt(fabs(segma[t])) << " ";
    }
    cout << endl << endl;

    r = 0;
    j = 0;
    stop = row;
    
    cout << "U Matrix: " << endl;
    for (int t = 0; t < row; t++) {
        while (j < stop) {
            cout << U[r++] << " ";
            j++;
        }
        j = 0;
        cout << endl;
    }
}