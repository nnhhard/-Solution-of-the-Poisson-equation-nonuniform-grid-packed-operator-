#include <QCoreApplication>
#include <math.h>
#include <crsmatrix.h>

#define N 10
#define M 10

double a_x = 0, b_x = 1;
double a_y = 0, b_y = 1;

const int flag = 1;

/*
 * if flag == 0 then solve by min result (Operator)
 * if flag == 1 then solve by min result (CSR Matrix)
 * if flag == 2 then solve by Gaus (Matrix)
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
    //return sin(M_PI * x) * sin(M_PI * y);
    return x * x + y * y;
}
double RightPart(int ln, int lm, double x, double y, int i, int j) {
    if( (i == 0 && j == 0) || (i == 0 && j == lm) || (i == ln && j == 0) || (i == ln && j == lm) )
           return Function(x, y);
       else if(i==0)
           return Function(x, y);
       else if(i == ln)
           return Function(x, y);
       else if(j == 0)
           return Function(x, y);
       else if(j == lm)
           return Function(x, y);
       else
           //return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
           return 4;
}


double Operator(double **vector, int ln, int lm, int i, int j, double *hx, double *hy) {
    if( (i == 0 && j == 0) || (i == 0 && j == lm) || (i == ln && j == 0) || (i == ln && j == lm) )
            return vector[i][j];
        else if(i==0)
            return vector[i][j];
        else if(i==ln)
            return vector[i][j];
        else if(j==0)
            return vector[i][j];
        else if(j==lm)
            return vector[i][j];
        else
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
                r[i][j] = Operator(u,ln,lm,i,j,hx,hy) - RightPart(ln,lm,x[i],y[j],i,j);


    for(int i = 0; i <= ln; i++)
        for(int j = 0; j <= lm; j++)
            if(mask_u[i][j] != 0)
                ar[i][j] = Operator(r,ln,lm,i,j,hx,hy);

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
                    r[i][j] = Operator(u1,ln,lm,i,j,hx,hy) - RightPart(ln,lm,x[i],y[j],i,j);

        for(int i = 0; i <= ln; i++)
            for(int j = 0; j <= lm; j++)
                if(mask_u[i][j] != 0)
                    ar[i][j]=Operator(r,ln,lm,i,j,hx,hy);

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

        printf("Norm = %lf\n",Norm(r,mat.m));
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
void getVectorOfRightPart(double *b, int ln, int lm, double *x, double *y) {
    double  **f = new double*[ln+1];
    for(int i = 0 ;i <= ln; i++)
        f[i] = new double[lm+1];

    for(int i = 0; i <= ln; i++) {
        for(int j = 0; j <= lm ; j++) {
            f[i][j] = RightPart(ln,lm,x[i],y[j],i,j);
        }
    }
    for(int i = 0; i <= ln; i++) {
        for(int j = 0; j <= lm; j++) {
            b[(lm+1)*i+j]=f[i][j];
        }
    }


    for(int i=0;i<=ln;i++) {
        delete []f[i];
    }
    delete []f;
}
void getMatrixOfOperator(double **A, int ln, int lm,double *hx, double *hy) {
    double *F1=new double[(ln + 1) * (lm + 1)];

    double  **e=new double*[ln + 1];
    for(int i = 0; i <= ln; i++)
        e[i] = new double[lm + 1];

    double  **u1=new double*[ln + 1];
    for(int i = 0; i<= ln; i++)
        u1[i] = new double[lm+1];

    int q = 0;
    int l = 0;
    for(int k = 0; k < ( (ln + 1) * (lm + 1)); k++)
    {

        for(int i = 0; i <= ln; i++)
        {
            for(int j = 0; j <= lm; j++)
            {
                e[i][j] = 0;
            }
        }
        e[q][l] = 1;
        l++;
        if(l > lm)
        {
            l = 0;
            q++;
        }


        for(int i = 0; i <= ln; i++)
        {
            for(int j = 0; j <= lm; j++)
            {
                u1[i][j] = Operator(e,ln,lm,i,j,hx,hy);
            }
        }


        for(int i = 0; i <= ln; i++)
        {
            for(int j = 0; j <= lm; j++)
            {
                F1[(lm + 1) * i + j] = u1[i][j];
            }
        }

        for(int i = 0; i < ((ln + 1) * (lm + 1)); i++)
        {
            A[i][k] = F1[i];
        }
    }
    delete []F1;

    for(int i = 0; i <= ln; i++)
    {
        delete []e[i];
        delete []u1[i];
    }
    delete []e;
    delete []u1;
}
void Solve(double **u, double **mask_u, int ln, int lm, double *x, double *y, double *hx, double *hy) {
    if(flag == 0) {
        printf("SolveByMinResultOperator: start\n");
            SolveByMinResultOperator(u,mask_u,ln,lm,x,y,hx,hy);
        printf("SolveByMinResultOperator: end\n");
    }
    else if(flag == 1) {
        printf("Get Matrix and Vector: start\n");
        int L = (ln + 1) * (lm + 1);
        double *b=new double[L];
        getVectorOfRightPart(b,ln,lm,x,y);

        double *vec_u=new double[L];

        double  **A=new double*[L];
        for(int i = 0; i < L; i++)
            A[i]=new double[L];

        getMatrixOfOperator(A,ln,lm,hx,hy);
        printf("Get Matrix and Vector: end\n");

        printf("Convert Matrix to CRS: start\n");
        CRS_Matrix mat = ConvertToCRS((const double**)A,L,L);
        printf("Convert Matrix to CRS: end\n");
        printf("SolveByMinResultMatrix: start\n");
        SolveByMinResultMatrix(mat,b,vec_u);
        printf("SolveByMinResultMatrix: end\n");
        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                u[i][j]=vec_u[(lm+1)*i+j];
            }
        }
    }
    else if(flag == 2) {
        printf("Get Matrix and Vector: start\n");
        int L = (ln + 1) * (lm + 1);
        double *b=new double[L];
        getVectorOfRightPart(b,ln,lm,x,y);

        double *vec_u=new double[L];

        double  **A=new double*[L];
        for(int i = 0; i < L; i++)
            A[i]=new double[L];

        getMatrixOfOperator(A,ln,lm,hx,hy);
        printf("Get Matrix and Vector: end\n");

        printf("SolveByGaus: start\n");
        SolveByGaus(A,b,vec_u,L);
        printf("SolveByGaus: end\n");
        for(int i=0;i<=ln;i++) {
            for(int j=0;j<=lm;j++) {
                u[i][j]=vec_u[(lm+1)*i+j];
            }
        }
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
            printf("%lf ",u[i][j]);
        }
        printf("\n");
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
