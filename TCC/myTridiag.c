#include "myTridiag.h"

/*
double* myTridiag(double a[], double b[], double c[], double d[], const size_t n) {
    // FUNCIONA TAMBEM
    double *w = (double*) calloc(n-1, sizeof(double));
    double *g = (double*) calloc(n, sizeof(double));
    double *p = (double*) calloc(n, sizeof(double));

    w[0] = c[0]/b[0];
    g[0] = d[0]/b[0];

    int i;
    for(i = 0; i < n-1; i++) {
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1]);
    }
    for(i = 0; i < n; i++) {
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1]);
    }
    p[n-1] = g[n-1];
    for(i = n-1; i > 0; i--) {
        p[i-1] = g[i-1] - w[i-1]*p[i];
    }
    free(w);
    free(g);
    return p;
}
*/

/*
double* myTridiag(double a[], double b[], double c[], double d[], int n) 
{
    double *y = (double*) calloc(n, sizeof(double));
    n--;
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i] * d[i + 1];
    }
    
    for(int i = 0; i <= n; i++) {
        y[i] = d[i];
    }
    return y;
}
*/


double* myTridiag(double b[], double a[], double c[], double f[], const size_t n) {
    // FUNCIONA !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double *v = (double *) calloc(n, sizeof(double)); // inicializando com 0s
    double *y = (double *) calloc(n, sizeof(double)); // inicializando com 0s
    
    double w = a[0];
    y[0] = f[0]/w;
    
    for(int i = 1; i < n; i++) {
        v[i-1] = c[i-1]/w;
        w = a[i] - b[i]*v[i-1];
        y[i] = (f[i] - b[i]*y[i-1])/w;
    }
    for(int j = n-2; j >= 0; j--) {
        y[j] = y[j] - v[j]*y[j+1];
    }
    free(v);
    return y;
    // source: MATLAB file tridiag.m
}


//  double * myTridiag(double a[], double b[], double c[], double x[], int N) {    // versao do iury
//    
//      /*
//      resolve Ax = v onde A é uma matriz tridiagonal composta pelos veores a, b, c
//      note que o conteúdo do vetor de entrada c será modificado, tornando esta uma função de um /único tempo de uso
//      x[] - inicialmente contém o vector de entrada v e retorna a solução x. indexados por[0, ..., N - 1]
//      N - número de equações
//      a[] - subdiagonal (diagonal abaixo da diagonal principal) -- indexados por [1, ..., N - 1]
//      b[] - matriz principal, indexados por [0, ..., N - 1]
//      c[] - superdiagonal (diagonal acima da diagonal principal) -- indexedos por [0, ..., N - 2]
//      */
//      double * y = (double *) malloc(N*sizeof(double ));   // Don't waste time initializing
//  
//      double m;
//      int n;
//  
//      c[0] = c[0] / b[0];
//      x[0] = x[0] / b[0];
//  
//      /* loop de 1 a N - 1 inclusive */
//      for (n = 1; n < N; n++) {
//        m = 1.0/(b[n] - a[n] * c[n - 1]);
//        c[n] = c[n] * m;
//        x[n] = (x[n] - a[n] * x[n - 1]) * m;
//      }
//  
//      /* loop de N - 2 a 0 inclusive */
//      for (n = N - 2; n >= 0; n-- )
//        x[n] = x[n] - c[n] * x[n + 1];
//  
//  
//      for (n = 0; n < N; n++)
//        y[n] = x[n];
//  
//      return y;
//  
//        // https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
//  }
/*
double* myTridiag(double e[], double f[], double g[], double r[], int n) {
    // input:
    // e = subdiagonal vector
    // f = diagonal vector
    // g = superdiagonal vector
    // r = right hand side vector
    // output:
    // x = solution vector
    
    // Memory allocation:
    //double* x = calloc(n, sizeof(double));  // Initialize with zeroes
    double* x = (double*) malloc(n*sizeof(double));   // Don't waste time initializing
    // Counter:
    int k;
    // *** Thomas algorithm ***
    // forward elimination
    double factor;
    for(k = 1; k < n; k++) {
      factor = e[k]/f[k-1];
      f[k] = f[k] - factor*g[k-1];
      r[k] = r[k] - factor*r[k-1];
    }
    // back substitution
    x[n-1] = r[n-1]/f[n-1];
    for(k = n-2; k >= 0; k--) {
      x[k] = (r[k] - g[k]*x[k+1])/f[k];
    }
    return x;
}
*/