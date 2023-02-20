#include <iostream>
#include <cmath>
#include <string>
#include <fstream>   // for file-reading
#include <sstream>   // for file-reading
#include <algorithm>
#include <vector>
#include <chrono>

#include "SVD.h"

#define MATRIX_A "matrixA.txt"
#define DATA "data.txt"
#define EXAMPLE "example7.txt"

/***********main*****************/
int main() {

    auto t1 = std::chrono::high_resolution_clock::now();
    size_t ROW=0, COL=0;
    int max_Iter = 1000;

    find_size(DATA, ROW, COL);//

    size_t row1 = ROW, col1 = COL;// size of matrix A
    size_t row2 = COL, col2 = ROW;//size of matrix AT

    double* A = new double[row1 * col1];

    read_csv(A, DATA);//

    double* AT = new double[row2 * col2];

    size_t k = 0;
    size_t i = 0;
    size_t j = 1;
    size_t d = 0;

    while (i < ROW * COL) {
        while (k < ROW * COL) {
            AT[i++] = A[k];
            k = k + COL;
        }
        k = j;
        j++;
    }

    int r = 0;
    size_t stop = col2;

    size_t sizeATA = row2 * col1;
    size_t sizeU = row1 * col2;
    size_t sizeSegma = sizeATA / col1;

    double* ATA = new double[sizeATA] {0};

    i = 0;
    k = 0;
    j = 0;
    d = 0;

    size_t stop1 = col1, stop2 = col2;
    size_t indJ = 1, indK = 0;

    while (i < sizeATA) {
        while (i < stop1) {
            while (d < stop2) {
                ATA[i] = ATA[i] + AT[k] * A[j];
                k++;
                j = j + col1;
                d++;
            }
            i++;
            k = indK;
            j = indJ;
            indJ++;
            d = 0;

        }
        stop1 = stop1 + col1;
        indK = indK + col2;
        k = indK;
        j = 0;
        indJ = 1;
    }

    size_t sizeS = col1;
    if (row1 > col1)
        sizeS = row1;

    SVD* obj_svd = new SVD(A, ATA, sizeATA, sizeU, sizeSegma, ROW, COL, max_Iter, DATA);
  
    obj_svd->find_svd();

    write_csv(obj_svd->VT, "VT_res.txt", COL, 1);
    write_csv(obj_svd->segmaTemp, "SEGMA_res.txt", COL, 2);
    write_csv(obj_svd->U, "U_res.txt", ROW, 3);

    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

    cout << "The total time: " << duration;
    cout << endl << endl;

    //printMartix(obj_svd->VT, obj_svd->segmaTemp, obj_svd->U, ROW, COL);//print result

    delete[] A;
    delete[] AT;
    delete[] ATA;
 
    return 0;
}

/**************************/
/*print matrxi AT and ATA*/
/*cout << "AT: " << endl;
for (int t = 0; t < row2; t++) {
    while (r < stop) {
        cout << AT[r++] << " ";
    }
    stop = stop + col2;
    cout << endl;
}*/

/*r = 0;
   stop = COL;
   cout << "matrix ATA: " << endl;
   for (int t = 0; t < COL; t++) {
       while (r < stop) {
           cout << ATA[r++] << " ";
       }
       stop = stop + COL;
       cout << endl;
 }*/