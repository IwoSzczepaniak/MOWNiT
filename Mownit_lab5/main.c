#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>


void naive(double** A, double** B, double** C,  int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

void better(double** A, double** B, double** C,  int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int k = 0; k < size; k++)
        {
            for (int j = 0; j < size; j++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

double** calloc_2d_array(int dim1, int dim2)
{
    double **C;
    C = (double**)calloc(dim1, sizeof(double));
    for (int i = 0; i < dim2; i++)
    {
        C[i] = (double*)calloc(dim2, sizeof(double));
    }
    return C;
}

void free2d(double** C, int dim1)
{
    for(int i = 0; i<dim1; i++)
    {
        free(C[i]);
    }
    free(C);
}

void blas(double *a, double *b, double *c, int size)
{
    gsl_matrix_view A = gsl_matrix_view_array(a, size, size);
    gsl_matrix_view B = gsl_matrix_view_array(b, size, size);
    gsl_matrix_view C = gsl_matrix_view_array(c, size, size);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);
}

double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

int main()
{
    double **A, **B, **C, *a, *b, *c;
    double time1, time2, time3;
    struct timespec start, end;
    clock_t resolution = sysconf(_SC_CLK_TCK);
    srand(time(NULL));   
    FILE *file = fopen("C_results.csv","a");
    fprintf(file, "size,naive,better,blas\n");
    
    for (int n = 5; n<510; n += 10)
    {
        C = calloc_2d_array(n,n);
        A = calloc_2d_array(n,n);
        B = calloc_2d_array(n,n);

        for(int j = 0; j < 10; j++)
        {  
            for (int a = 0; a < n; a++)
            {
                for(int b = 0; b < n; b++)
                {
                    A[a][b] = drand(0,1);
                    B[a][b] = drand(0,1);
                }
            }

            clock_gettime(CLOCK_MONOTONIC, &start);
            naive(A,B,C,n); 
            clock_gettime(CLOCK_MONOTONIC, &end);
            time1 = (double) (end.tv_sec-start.tv_sec) + (double)(end.tv_nsec-start.tv_nsec)/10e9;

            clock_gettime(CLOCK_MONOTONIC, &start);
            better(A,B,C,n); 
            clock_gettime(CLOCK_MONOTONIC, &end);
            time2 = (double) (end.tv_sec-start.tv_sec) + (double)(end.tv_nsec-start.tv_nsec)/10e9;

            a = calloc(n*n, sizeof(double));
            b = calloc(n*n, sizeof(double));
            c = calloc(n*n, sizeof(double));

            for (int k = 0; k < n*n; k++){
                a[k] = A[k/n][k%n];
                b[k] = B[k/n][k%n];
            }   


            clock_gettime(CLOCK_MONOTONIC, &start);
            blas(a,b,c,n); 
            clock_gettime(CLOCK_MONOTONIC, &end);

            time3 = (double) (end.tv_sec-start.tv_sec) + (double)(end.tv_nsec-start.tv_nsec)/10e9;

            fprintf(file, "%d,%f, %f, %f\n", n, time1, time2, time3);
        }

        free2d(A, n);
        free2d(B, n);
        free2d(C, n);
        free(a);
        free(b);
        free(c);
    }

    return 0;
}