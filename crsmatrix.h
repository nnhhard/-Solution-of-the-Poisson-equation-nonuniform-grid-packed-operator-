#ifndef SRSMATRIX_H
#define SRSMATRIX_H

struct CRS_Matrix
{
    double* values;
    int* cols;
    int* srows;
    int m;
    int n;
    int nz;
};

void DeallocCRS(CRS_Matrix& m);
CRS_Matrix AllocCRS(int m, int n, int countNonZero);
CRS_Matrix ConvertToCRS(const double** denseMatrix, int m, int n);

#endif // SRSMATRIX_H
