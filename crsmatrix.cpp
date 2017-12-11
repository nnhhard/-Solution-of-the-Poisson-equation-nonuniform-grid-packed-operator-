#include <crsmatrix.h>
#include <iostream>



CRS_Matrix AllocCRS(int m, int n, int countNonZero)
{
    if(countNonZero > m*n)
        throw std::runtime_error("nz > m*n");
    CRS_Matrix mat;
    mat.values = new double[countNonZero];
    mat.cols   = new int[countNonZero];
    mat.srows  = new int[m+1];
    mat.m = m;
    mat.n = n;
    mat.nz = countNonZero;
    return mat;
}
void DeallocCRS(CRS_Matrix& m)
{
    delete[] m.values;
    delete[] m.cols;
    delete[] m.srows;
}
CRS_Matrix ConvertToCRS(const double** denseMatrix, int m, int n) {
    int nzLoc = 0;
    for(int i = 0;i < m;i++)
        for(int j = 0;j < n;j++)
            if(denseMatrix[i][j] != 0)
                ++nzLoc;


    CRS_Matrix mat = AllocCRS(m, n, nzLoc);

    mat.cols[0] = 1;
    mat.cols[1] = 1;

    int valueIndex = 0;
    mat.srows[0] = 0;
    for(int i = 0;i < m;i++) {
        int countNZinCol = 0;
        for(int j = 0;j < n;j++) {
            if(denseMatrix[i][j] != 0) {
                mat.values[valueIndex] = denseMatrix[i][j];
                mat.cols[valueIndex] = j;
                ++valueIndex;
                ++countNZinCol;
            }
        }
        mat.srows[i+1] = mat.srows[i] + countNZinCol;
    }
    return mat;
}
