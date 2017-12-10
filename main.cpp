#include <QCoreApplication>
#include <math.h>



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
}
double RightPart(double x, double y) {
    return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
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

    for(int i = 0; i <= n; i++){
        delete []u1[i];
        delete []r[i];
        delete []ar[i];
    }

    delete []u1;
    delete []r;
    delete []ar;
}
void Solve(double **u, double **mask_u, int ln, int lm, double *x, double *y, double *hx, double *hy) {
    SolveByMinResultOperator(u,mask_u,n,m,x,y,hx,hy);
}

#define n 14
#define m 10

double a_x = 0, b_x = 1;
double a_y = 0, b_y = 1;

int main()
{
    double *hx = new double[n];
    double *hy = new double[m];

    double h = (b_x - a_x) / (n);

    if(0 == 0)
    {
        hx[0] = h/4.0;
        hx[n-1] = h/4.0;

        for(int i = 1; i <= n-2; i++)
        {
            hx[i] = ( b_x - h/2.0 - a_x )/(n-4);
        }
    }
    else
    {
        for(int i = 0; i <= n-1; i++)
        {
            hx[i] = h;
        }
    }

    h = (b_y - a_y) / (m);
    for(int i = 0; i <= m-1; i++)
    {
        hy[i] = h;
    }

    double *x = new double[n+1];

    x[0]=0;
    for(int i = 1; i <= n; i++)
    {
        x[i] = x[i-1] + hx[i-1];
    }

    double *y = new double[m+1];

    y[0]=0;
    for(int i = 1; i <= m; i++)
    {
        y[i] = y[i-1] + hy[i-1];
    }

    double **u = new double*[n+1];
    for(int i = 0; i <= n; i++)
        u[i] = new double[m+1];

    double **ut = new double*[n+1];
    for(int i = 0; i <= n; i++)
        ut[i] = new double[m+1];

    double **w = new double*[n+1];
    for(int i = 0; i <= n; i++)
        w[i] = new double[m+1];

    double **mask_u = new double*[n+1];
    for(int i = 0; i <= n; i++)
        mask_u[i] = new double[m+1];

    null(mask_u,n,m);
    null(w,n,m);

    for(int i = 1; i <= n-1; i++){
        for(int j = 1; j <= m-1; j++){
            mask_u[i][j] = 1;
        }
    }

    for(int i = 0; i <= n; i++){
        for(int j = 0; j <= m; j++){
            ut[i][j] = Function(x[i], y[j]);
        }
    }

    Solve(u,mask_u,n,m,x,y,hx,hy);

    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= m; j++)
            if(mask_u[i][j] != 0)
                w[i][j] = u[i][j] - ut[i][j];

    printf("Norm minR = %lf \n",Norm(w,n,m));

    printf("\nFinalisation work\n");

    for(int i = 0; i <= n; i++){
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
