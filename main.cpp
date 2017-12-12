#include <QCoreApplication>
#include <math.h>
#include <crsmatrix.h>

#define N 20
#define M 20

double a_x = 0, b_x = 1;
double a_y = 0, b_y = 1;

const int flag = 1;

/*
 * if flag == 0 then solve by min result (Operator)
 * if flag == 1 then solve by Gaus (Matrix)
 * if flag == 2 then solve by min result (CSR Matrix)
 */

double Norm(double *vector, int ln) {
    double max = fabs(vector[0]);

    for(int i = 0; i < ln; i++) {
        if(fabs(vector[i])>max) {
            max=fabs(vector[i]);
        }
    }
    return max;
}
double Sc(double *vector_a, double *vector_b, int ln) {
    double s = 0;
    for(int i = 0; i < ln; i++) {
            s += vector_a[i] * vector_b[i];
    }
    return s;
}

double Norm(double **vector, int ln, int lm) {
    double max = fabs(vector[0][0]);

    for(int i = 0; i <= ln; i++) {
        for(int j = 0; j <= lm; j++) {
            if(fabs(vector[i][j])>max) {
                max=fabs(vector[i][j]);
            }
        }
    }
    return max;
}
double Sc(double **vector_a, double **vector_b, int ln, int lm) {
    double s = 0;
    for(int i = 0; i <= ln; i++) {
        for(int j = 0; j <= lm; j++) {
            s += vector_a[i][j] * vector_b[i][j];
        }
    }
    return s;
}
void null(double **vector, int ln, int lm) {
    for(int i = 0; i <= ln; i++)
        for(int j = 0; j <= lm; j++)
            vector[i][j]=0;
}
double Function(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
    //return x * x + y * y;
}
double RightPart(double x, double y) {
    return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    //return 4;
}
double Operator(double **vector, int i, int j, double *hx, double *hy) {
    return ( ( ( ( vector[i+1][j] - vector[i][j] ) / (hx[i]) ) - ( ( vector[i][j] - vector[i-1][j] ) / (hx[i-1] ) ) ) / ( ( hx[i] + hx[i-1] ) / 2.0 ) )
         + ( ( ( ( vector[i][j+1] - vector[i][j] ) / (hy[j]) ) - ( ( vector[i][j] - vector[i][j-1] ) / (hy[j-1] ) ) ) / ( ( hy[j] + hy[j-1] ) / 2.0 ) );
}
void SolveByMinResultOperator(double **u, double **mask_u, int ln, int lm, double *x, double *y, double *hx, double *hy) {
    for(int i = 0; i <= ln; i++){
        for(int j = 0; j <= lm; j++){
            if(mask_u[i][j] == 0){
                u[i][j] = Function(x[i], y[j]);
            }
            else{
                u[i][j] = 10;
            }
        }
    }

    double **u1 = new double*[ln+1];
    for(int i = 0; i <= ln; i++)
        u1[i] = new double[lm+1];

    double **r = new double*[ln+1];
    for(int i = 0; i <= ln; i++)
        r[i] = new double[lm+1];

    double **ar = new double*[ln+1];
    for(int i = 0; i <= ln; i++)
        ar[i] = new double[lm+1];

    null(r,ln,lm);
    null(ar,ln,lm);


    for(int i = 0; i <= ln; i++)
        for(int j = 0; j <= lm; j++)
            if(mask_u[i][j] != 0)
                r[i][j] = Operator(u,i,j,hx,hy) - RightPart(x[i],y[j]);


    for(int i = 0; i <= ln; i++)
        for(int j = 0; j <= lm; j++)
            if(mask_u[i][j] != 0)
                ar[i][j] = Operator(r,i,j,hx,hy);

    double tau = Sc(r,ar,ln,lm)/Sc(ar,ar,ln,lm);

    double eps = 1e-6;

    while(Norm(r,ln,lm)>eps)
    {
        for(int i = 0; i <= ln; i++)
            for(int j = 0; j <= lm; j++)
                u1[i][j]=u[i][j]-tau*r[i][j];

        for(int i = 0; i <= ln; i++)
            for(int j = 0; j <= lm; j++)
                if(mask_u[i][j] != 0)
                    r[i][j] = Operator(u1,i,j,hx,hy) - RightPart(x[i],y[j]);

        for(int i = 0; i <= ln; i++)
            for(int j = 0; j <= lm; j++)
                if(mask_u[i][j] != 0)
                    ar[i][j]=Operator(r,i,j,hx,hy);

        tau = Sc(r,ar,ln,lm)/Sc(ar,ar,ln,lm);

        for(int i = 0; i <= ln; i++)
            for(int j = 0; j <= lm; j++)
                if(mask_u[i][j] != 0)
                    u[i][j]=u1[i][j];
    }

    for(int i = 0; i <= ln; i++){
        delete []u1[i];
        delete []r[i];
        delete []ar[i];
    }

    delete []u1;
    delete []r;
    delete []ar;
}
void SolveByMinResultMatrix(CRS_Matrix &mat, double *b, double *x) {

    for(int i = 0; i < mat.m; i++) {
        x[i] = 0;
    }

    double *r = new double[mat.m];
    double *ar = new double[mat.m];

    for(int i = 0;i < mat.m;i++) {
            r[i] = 0;
            for(int j = mat.srows[i];j < mat.srows[i+1];j++)
                r[i] += mat.values[j] * x[mat.cols[j]];
            r[i] = r[i] - b[i];
        }

    for(int i = 0;i < mat.m;i++) {
            ar[i] = 0;
            for(int j = mat.srows[i];j < mat.srows[i+1];j++)
                ar[i] += mat.values[j] * r[mat.cols[j]];
        }

    double tau = Sc(r,ar,mat.m)/Sc(ar,ar,mat.m);

    double eps = 1e-9;


    int n = 0;
    while(Norm(r,mat.m) > eps)
    {
        n++;
        for(int i = 0; i < mat.m; i++) {
            x[i] = x[i] - tau * r[i];
        }

        for(int i = 0;i < mat.m;i++) {
                r[i] = 0;
                for(int j = mat.srows[i];j < mat.srows[i+1];j++)
                    r[i] += mat.values[j] * x[mat.cols[j]];
                r[i] = r[i] - b[i];
            }

        for(int i = 0;i < mat.m;i++) {
                ar[i] = 0;
                for(int j = mat.srows[i];j < mat.srows[i+1];j++)
                    ar[i] += mat.values[j] * r[mat.cols[j]];
            }

        tau = Sc(r,ar,mat.m)/Sc(ar,ar,mat.m);

        //printf("Norm = %lf\n",Norm(r,mat.m));
    }


}
void SolveByGaus(double **A,double *B,double *x,int ln) {

    for ( int k=0; k<ln-1; k++)
    {
        for ( int i=k+1; i<ln; i++)
        {
            double t=A[i][k]/A[k][k];
            B[i]=B[i]-t*B[k];
            for (int j=k; j<ln; j++)
                A[i][j]=A[i][j]-t*A[k][j];
        }
    }



    for(int i=ln-1;i>=0;i--)
    {
        double sum=0;
        for(int j=i+1;j<ln;j++)
            sum+=A[i][j]*x[j];
        x[i]=(B[i]-sum)/A[i][i];
    }

}
double** getMatrixOfOperator(int L, double **mask, int ln, int lm,double *hx, double *hy) {
        double  **A=new double*[L];
        for(int i = 0; i < L; i++)
            A[i]=new double[L];

        double **e = new double*[ln + 1];
        for(int i = 0; i <= ln; i++) {
            e[i] = new double[lm + 1];
        }

        double **a = new double*[ln + 1];
        for(int i = 0; i <= ln; i++) {
            a[i] = new double[lm + 1];
        }

        double **mask_temp = new double*[ln + 1];
        for(int i = 0; i <= ln; i++) {
            mask_temp[i] = new double[lm + 1];
        }

        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                mask_temp[i][j] = mask[i][j];
            }
        }

        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                if(mask[i][j] == 0)
                    e[i][j] = 0;
            }
        }

        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                e[i][j] = 0;
            }
        }

        bool flag_break = 0;
        int i_fixed = 0;
        int j_fixed = 0;

        for(int k = 0; k < L; k++)
        {

            for(int i = 0; i <= ln; i++) {
                for(int j = 0; j <= lm; j++) {
                    if(mask_temp[i][j] != 0) {
                        e[i][j] = 1;
                        mask_temp[i][j] = 0;
                        i_fixed = i;
                        j_fixed = j;
                        flag_break = 1;
                        break;
                    }
                }
                if(flag_break == 1)
                    break;
            }


            for(int i = 0; i <= ln; i++) {
                for(int j = 0; j <= lm; j++) {
                    a[i][j] = 0;
                    if(mask[i][j] != 0){
                        a[i][j] = Operator(e,i,j,hx,hy);

                    }
                }
            }

            int q = 0;
            for(int i = 0; i <= ln; i++) {
                for(int j = 0; j <= lm; j++) {
                    if(mask[i][j] != 0) {
                        A[q][k] = a[i][j];
                        q++;
                    }
                }
            }

            flag_break = 0;
            e[i_fixed][j_fixed] = 0;
        }

        for(int i = 0; i <= lm; i++) {
            delete []e[i];
            delete []a[i];
            delete []mask_temp[i];
        }
        delete []e;
        delete []a;
        delete []mask_temp;
        e = NULL;
        a = NULL;
        mask_temp = NULL;
        return A;
}
double* getVectorOfRightPart(int L,double **mask, int ln, int lm, double *x, double *y,double *hx,double *hy) {
    double *b=new double[L];

    double **f = new double*[ln+1];
    for(int i=0;i<=ln;i++) {
        f[i] = new double[lm+1];
    }

    for(int i=0;i<=ln;i++) {
        for(int j=0;j<=lm;j++) {
            if(mask[i][j] != 0)
                f[i][j]=RightPart(x[i],y[j]);
            else
                f[i][j]=0;
        }
    }

    double **z = new double*[ln+1];
    for(int i=0;i<=ln;i++) {
        z[i] = new double[lm+1];
    }
    for(int i=0;i<=ln;i++) {
        for(int j=0;j<=lm;j++) {
            if(mask[i][j] == 0)
                z[i][j] = Function(x[i],y[j]);
            else
                z[i][j] = 0;
        }
    }

    for(int i=0;i<=ln;i++) {
        for(int j=0;j<=lm;j++) {
            if(mask[i][j] != 0)
                f[i][j] = - Operator(z,i,j,hx,hy) + f[i][j];
        }
    }

    int k=0;
    for(int i=0;i<=ln;i++) {
        for(int j=0;j<=lm;j++) {
            if(mask[i][j] != 0) {
                b[k] = f[i][j];
                k++;
            }
        }
    }
    for(int i = 0; i <= lm; i++) {
        delete []f[i];
        delete []z[i];
    }
    delete []f;
    delete []z;
    f = NULL;
    z = NULL;
    return b;
}
void Solve(double **u, double **mask_u, int ln, int lm, double *x, double *y, double *hx, double *hy) {
    if(flag == 0) {
        printf("SolveByMinResultOperator: start\n");
            SolveByMinResultOperator(u,mask_u,ln,lm,x,y,hx,hy);
        printf("SolveByMinResultOperator: end\n");
    }
    else if(flag == 1) {
        printf("Get Matrix and Vector: start\n");
        int L = (ln + 1 - 2) * (lm + 1 - 2);

        double *b = getVectorOfRightPart(L,mask_u,ln,lm,x,y,hx,hy);


        double  **A = getMatrixOfOperator(L,mask_u,ln,lm,hx,hy);

        double *vec_u=new double[L];
        printf("Get Matrix and Vector: end\n");

        printf("SolveByGaus: start\n");
        SolveByGaus(A,b,vec_u,L);
        printf("SolveByGaus: end\n");
        int k=0;
        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                if(mask_u[i][j] != 0) {
                    u[i][j] = vec_u[k];
                    k++;
                }
                else {
                    u[i][j] = Function(x[i],y[j]);
                }
            }
        }
        for(int i = 0; i < L; i++) {
            delete []A[i];
        }
        delete []A;
        delete []vec_u;
        delete []b;
        A = NULL;
        b = NULL;
        vec_u = NULL;
    }
    else if(flag == 2) { 
        printf("Get Matrix and Vector: start\n");
        int L = (ln + 1 - 2) * (lm + 1 - 2);
        double *b = getVectorOfRightPart(L,mask_u,ln,lm,x,y,hx,hy);



        double  **A = getMatrixOfOperator(L,mask_u,ln,lm,hx,hy);

        double *vec_u=new double[L];
        printf("Get Matrix and Vector: end\n");

        printf("Convert Matrix to CRS: start\n");
        CRS_Matrix mat = ConvertToCRS((const double **)A,L,L);
        printf("Convert Matrix to CRS: end\n");

        printf("SolveByMinResultMatrixCRS: start\n");
        SolveByMinResultMatrix(mat,b,vec_u);
        printf("SolveByMinResultMatrixCRS: end\n");

        int k=0;
        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                if(mask_u[i][j] != 0) {
                    u[i][j] = vec_u[k];
                    k++;
                }
                else {
                    u[i][j] = Function(x[i],y[j]);
                }
            }
        }
        for(int i = 0; i < L; i++) {
            delete []A[i];
        }
        delete []A;
        delete []vec_u;
        delete []b;
        A = NULL;
        b = NULL;
        vec_u = NULL;
    }
    else {
        printf("Error solvers\n");
    }

}



int main()
{
    double *hx = new double[N];
    double *hy = new double[M];

    double h = (b_x - a_x) / (N);

    if(0 == 1)
    {
        hx[0] = h/4.0;
        hx[N-1] = h/4.0;

        for(int i = 1; i <= N-2; i++)
        {
            hx[i] = ( b_x - h/2.0 - a_x )/(N-2);
        }
    }
    else
    {
        for(int i = 0; i <= N-1; i++)
        {
            hx[i] = h;
        }
    }

    h = (b_y - a_y) / (M);
    for(int i = 0; i <= M-1; i++)
    {
        hy[i] = h;
    }

    double *x = new double[N+1];

    x[0]=0;
    for(int i = 1; i <= N; i++)
    {
        x[i] = x[i-1] + hx[i-1];
    }

    double *y = new double[M+1];

    y[0]=0;
    for(int i = 1; i <= M; i++)
    {
        y[i] = y[i-1] + hy[i-1];
    }

    double **u = new double*[N+1];
    for(int i = 0; i <= N; i++)
        u[i] = new double[M+1];

    double **ut = new double*[N+1];
    for(int i = 0; i <= N; i++)
        ut[i] = new double[M+1];

    double **w = new double*[N+1];
    for(int i = 0; i <= N; i++)
        w[i] = new double[M+1];

    double **mask_u = new double*[N+1];
    for(int i = 0; i <= N; i++)
        mask_u[i] = new double[M+1];

    null(mask_u,N,M);
    null(w,N,M);

    for(int i = 1; i <= N-1; i++){
        for(int j = 1; j <= M-1; j++){
            mask_u[i][j] = 1;
        }
    }

    for(int i = 0; i <= N; i++){
        for(int j = 0; j <= M; j++){
            ut[i][j] = Function(x[i], y[j]);
        }
    }

    Solve(u,mask_u,N,M,x,y,hx,hy);

    for(int i = 0; i <= N; i++) {
        for(int j = 0; j <= M; j++) {
            if(mask_u[i][j] != 0)
                w[i][j] = u[i][j] - ut[i][j];
        }
    }

    printf("Norm(W) = %lf \n",Norm(w,N,M));

    printf("\nFinalisation work\n");

    for(int i = 0; i <= N; i++){
        delete []u[i];
        delete []ut[i];
        delete []w[i];
        delete []mask_u[i];
    }

    delete []u;
    delete []ut;
    delete []w;
    delete []mask_u;
    return 0;
}
