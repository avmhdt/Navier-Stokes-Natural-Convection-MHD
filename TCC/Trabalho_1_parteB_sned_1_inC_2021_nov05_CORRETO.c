#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 0.0
#define B 1.0
#define Nel 256
#define h (((double)B - (double)A)/(double)Nel)

#define I Nel+1
#define J I

#define Re (double) 2500.0

#define which_h 1
/***************************************/
/******** Function Declarations ********/

double max(double num1, double num2);

double min(double num1, double num2);

// double* seq(double a, double b, double by);

double* seq(double start, double finish, int size);

double max_vec(double vec[], int size);

double min_vec(double vec[], int size);

int which_max(double vec[], int size);

int which_min(double vec[], int size);

void print_results(double **u, int rows, int cols, char filename[]);

// void plot(double x[], double y[]);

// double* myTridiag(double e[], double f[], double g[], double r[], int n);

double** createArray(int n, int m);

void cavity();

/***************************************/
/**************** Utils ****************/

double max(double num1, double num2) {
    return (num1 > num2 ) ? num1 : num2;
}

double min(double num1, double num2) {
    return (num1 > num2 ) ? num2 : num1;
}

/*
double* seq(double a, double b, double by) {
    int n = ceil((b-a)/by);
    double *x = calloc(n+1, sizeof(double));
    x[0] = (double)a; x[n] = (double)b;
    for(int i = 1; i < n; i++) {
        x[i] = a + (double)i*by;
    }
    return x;
}
*/

double* seq(double start, double finish, int size) {
    double by = ((double)finish - (double)start)/(double)(size-1);
    // double *x = calloc(size, sizeof(double));
    double *x = (double*) malloc(size*sizeof(double));
    x[0] = (double)start; x[size-1] = (double)finish;
    for(int i = 1; i < size-1; i++) {
        x[i] = start + i*by;
    }
    return x;
}

double max_vec(double vec[], int size) {
    double m = vec[0];
    for(int i = 1; i < size; i++) {
        if(vec[i] > m)
            m = vec[i];
    }
    return m;
}

double min_vec(double vec[], int size) {
    double m = vec[0];
    for(int i = 1; i < size; i++) {
        if(vec[i] < m)
            m = vec[i];
    }
    return m;
}

double max_mat(double mat[][J-2], int n, int m) {
    double thisMax = mat[0][0];
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(mat[i][j] > thisMax) {
                thisMax = mat[i][j];
            }
        }
    }
    return thisMax;
}

double min_mat(double mat[][J-2], int n, int m) {
    double thisMin = mat[0][0];
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(mat[i][j] < thisMin) {
                thisMin = mat[i][j];
            }
        }
    }
    return thisMin;
}

int which_max(double vec[], int size) {
    int m = 0;
    for(int i = 1; i < size; i++) {
        if(vec[i] > vec[m])
            m = i;
    }
    return m;
}

int which_min(double vec[], int size) {
    int m = 0;
    for(int i = 1; i < size; i++) {
        if(vec[i] < vec[m])
            m = i;
    }
    return m;
}


void print_results(double **u, int rows, int cols, char filename[]) {
    FILE* f = fopen(filename, "w");
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols-1; j++) {
            fprintf(f, "%.20f,", u[i][j]);
        }
        fprintf(f, "%.20f", u[i][cols-1]);
        fprintf(f, "\n");
    }
    fclose(f);
    return;
}

// void plot(double* x, double* y);
/*
double* myTridiag(double a[], double b[], double c[], double d[], int n) {
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
double** createArray(int n, int m) {
    double* values = (double*) calloc(n*m, sizeof(double)); // Initialize with zeroes
    double** rows = malloc(n*sizeof(double*));
    for(int i = 0; i < n; i++) {
        rows[i] = values + i*m;
    }
    return rows;
}

void destroyArray(double** arr) {
    free(*arr);
    free(arr);
}

void zeroArray(double** arr, int n, int m) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            arr[i][j] = 0.0;
        }
    }
}

/**********************************************************/
/***** Solving quasi-transient cavity problem via ADI *****/

void cavity() {
    double *x, *y;
    x = seq(A, B, I);
    y = seq(A, B, J);
    
    double dt;
    switch(which_h) {
        case 1:
            dt = pow(h, 2.0); // Nel = 32, Re = 1000 Final n = 221
            break;
        case 2:
            dt = pow(h, 2.0)/2.0;
            break;
        case 3:
            dt = pow(h, 2.0)/4.0;
            break;
        case 4:
            dt = pow(h, 2.0)/10.0; // Nel = 32, Re = 1000, Final n = 2205
            break;
        case 5:
            dt = pow(h, 3.0); // FUNCIONA COM Nel = 32, Re = 1000, Final n = 7055!!!
            break;
        case 6:
            dt = pow(h, 4.0); // FUNCIONA COM Nel = 32, Re = 1000, Final n = 225748!!!
            break;
        default:
            printf("ERROR: Invalid which_h\n");
    }

    
    double sig;
    sig = dt/(2*pow(h, 2));
    
    double **w[2], **psi[2];
    double **what, **psihat;
    int n; // time counter
    int i, j; // x and y dimension counters
    
    int k; // other counter
    
    for(k = 0; k < 2; k++) {
        w[k] = createArray(I, J);
        psi[k] = createArray(I, J);
    }
    what = createArray(I, J);
    psihat = createArray(I, J);
    
    /* Vectors from AU = F linear system of equations from ADI method */
    double F[J-2] = {0};           
    double Di[J-2] = {0}, Dp[J-2] = {0}, Ds[J-2] = {0};   // Diagonals from A matrix
    double *U = NULL;
    // RhoX and RhoY:
    double RhoX[J-2] = {0}, RhoY[J-2] = {0};
    
    // Residuals:
    double Rpsi[I-2][J-2] = {{0}}, Rw[I-2][J-2] = {{0}};
    double R, maxRpsi, maxRw;
    
    // Initial conditions:
    zeroArray(w[0], I, J);
    zeroArray(psi[0], I, J);
    zeroArray(what, I, J);
    zeroArray(psihat, I, J);
    
    /* Boundary conditions for vorticity */
    /* // Tentando impor o contorno na condicao inicial
    for(j = 1; j < J-1; j++) {
        w[0][0][j] = -2*psi[0][1][j]/pow(h, 2);             // Left
        w[0][I-1][j] = -2*psi[0][I-2][j]/pow(h, 2);         // Right
    }
    for(i = 0; i < I; i++) {
        w[0][i][0] = -2*psi[0][i][1]/pow(h, 2);             // Inferior
        w[0][i][J-1] = -2*psi[0][i][J-2]/pow(h, 2) - 2/h;   // Superior
    }
    */
    
    int stop = 0; // boolean to check stop criteria
    n = 1;
    
    while(!stop) {
        printf("**************** n = %d\n", n);
        
        /*************** ADI method ***************/
        
        /************ Stream function: ************/
        /* Boundary conditions for stream function */
        for(i = 0; i < I; i++) {
            psi[1][i][0] = 0.0;
            psihat[i][0] = 0.0;
            psi[1][i][J-1] = 0.0;
            psihat[i][J-1] = 0.0;
        }
        for(j = 0; j < J; j++) {
            psi[1][0][j] = 0.0;
            psihat[0][j] = 0.0;
            psi[1][I-1][j] = 0.0;
            psihat[I-1][j] = 0.0;
        }
        
        //printf("line %d\n", 262);
        
        // Filling diagonals from A matrix first
        for(k = 0; k < J-2; k++) {
            Di[k] = -sig;
            Dp[k] = 1 + 2*sig;
            Ds[k] = -sig;  
        }
        
        //printf("line %d\n", 271);
        
        /******* Step 1: Implicit x *******/
        for(j = 1; j < J-1; j++) {
            // Filling source term F:
            for(i = 1; i < I-1; i++) {
                F[i-1] = (dt/2)*w[0][i][j] + (1-2*sig)*psi[0][i][j] + sig*(psi[0][i][j+1] + psi[0][i][j-1]);
            }
            F[0] = F[0] + sig*psihat[0][j];
            F[I-3] = F[I-3] + sig*psihat[I-1][j];
            
            U = myTridiag(Di, Dp, Ds, F, I-2);
            for(i = 1; i < I-1; i++) {
                psihat[i][j] = U[i-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 290);
        
        /******* Step 2: Implicit y *******/
        for(i = 1; i < I-1; i++) {
            // Filling source term F:
            for(j = 1; j < J-1; j++) {
                F[j-1] = (dt/2)*w[0][i][j] + (1-2*sig)*psihat[i][j] + sig*(psihat[i+1][j] + psihat[i-1][j]);
            }
            F[0] = F[0] + sig*psi[1][i][0];
            F[J-3] = F[J-3] + sig*psi[1][i][J-1];
            
            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                psi[1][i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 309);
        
        /*************** Vorticity: ***************/
        /* Boundary conditions for vorticity */
        for(j = 1; j < J-1; j++) { // j = 1,...J-2 para nao preencher o contorno superior nem inferior, loop abaixo     // ESSE LOOP TEM QUE SER PRIMEIRO 
            w[1][0][j] = -2*psi[1][1][j]/pow(h, 2);             // Left
            what[0][j] = -2*psi[1][1][j]/pow(h, 2);             // Left
            w[1][I-1][j] = -2*psi[1][I-2][j]/pow(h, 2);         // Right
            what[I-1][j] = -2*psi[1][I-2][j]/pow(h, 2);         // Right
        }
        for(i = 0; i < I; i++) {      /// ESSE LOOP TEM QUE SER POR ULTIMO PRA QUE O CONTORNO SUPERIOR SEJA PREENCHIDO.
            w[1][i][0] = -2*psi[1][i][1]/pow(h, 2);             // Inferior
            what[i][0] = -2*psi[1][i][1]/pow(h, 2);             // Inferior
            w[1][i][J-1] = -2*psi[1][i][J-2]/pow(h, 2) - 2.0/h;   // Superior
            what[i][J-1] = -2*psi[1][i][J-2]/pow(h, 2) - 2.0/h;   // Superior
        }
        
        //printf("line %d\n", 326);
        
        /******* Step 1: Implicit x *******/
        for(j = 1; j < J-1; j++) {
            for(i = 1; i < I-1; i++) {
                RhoX[i-1] = Re*(psi[1][i+1][j] - psi[1][i-1][j])/4;
                RhoY[i-1] = Re*(psi[1][i][j+1] - psi[1][i][j-1])/4;
                
                Di[i-1] = -sig*(RhoY[i-1] + 1);
                Dp[i-1] = 1 + 2*sig;
                Ds[i-1] = sig*(RhoY[i-1] - 1);
                
                F[i-1] = (1-2*sig)*w[0][i][j] + sig*(RhoX[i-1] + 1)*w[0][i][j+1] - sig*(RhoX[i-1] - 1)*w[0][i][j-1];
            }
            
            F[0] = F[0] + sig*(RhoY[0] + 1)*what[0][j];
            F[I-3] = F[I-3] - sig*(RhoY[I-3] - 1)*what[I-1][j];
            
            U = myTridiag(Di, Dp, Ds, F, I-2);
            for(i = 1; i < I-1; i++) {
                what[i][j] = U[i-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 352);
        
        /******* Step 2: Implicit y *******/
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                RhoX[j-1] = Re*(psi[1][i+1][j] - psi[1][i-1][j])/4;
                RhoY[j-1] = Re*(psi[1][i][j+1] - psi[1][i][j-1])/4;
                
                Di[j-1] = sig*(RhoX[j-1] - 1);
                Dp[j-1] = 1 + 2*sig;
                Ds[j-1] = -sig*(RhoX[j-1] + 1);
                
                F[j-1] = (1-2*sig)*what[i][j] - sig*(RhoY[j-1] - 1)*what[i+1][j] + sig*(RhoY[j-1] + 1)*what[i-1][j];
            }
            
            F[0] = F[0] - sig*(RhoX[0] - 1)*w[1][i][0];
            F[J-3] = F[J-3] + sig*(RhoX[J-3] + 1)*w[1][i][J-1];
            
            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                w[1][i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 378);
        
        /* Calculating approximation residuals */
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                Rpsi[i-1][j-1] = fabs((psi[1][i+1][j] + psi[1][i-1][j] + psi[1][i][j+1] + psi[1][i][j-1] - 4*psi[1][i][j])/pow(h, 2) + w[1][i][j]);
                
                Rw[i-1][j-1] = fabs((w[1][i+1][j] + w[1][i-1][j] + w[1][i][j+1] + w[1][i][j-1] - 4*w[1][i][j])/(Re*pow(h, 2)) - ((psi[1][i][j+1] - psi[1][i][j-1])*(w[1][i+1][j] - w[1][i-1][j]) - (psi[1][i+1][j] - psi[1][i-1][j])*(w[1][i][j+1] - w[1][i][j-1]))/(4*pow(h, 2)));
                
                /* Rpsi[i-1][j-1] = fabs((psi[1][i+1][j] - 2*psi[1][i][j] + psi[1][i-1][j])/pow(h, 2) +
                                 (psi[1][i][j+1] - 2*psi[1][i][j] + psi[1][i][j-1])/pow(h, 2) +
                                 w[1][i][j]);
                */
                
                // printf("\t%f", Rw[i-1][j-1]);
            }
            //printf("\n");
            //printf("max(Rpsi[%d][j]) = %f\nmax(Rw[%d][j]) = %f\n", i-1, i-1, max_vec(Rpsi[i-1], J-2), max_vec(Rw[i-1], J-2));
        }
        
        //printf("line %d\n", 389);
        
        maxRpsi = max_mat(Rpsi, I-2, J-2);
        maxRw = max_mat(Rw, I-2, J-2);
        R = max(maxRpsi, maxRw);
        //R = min(maxRpsi, maxRw);
        printf("max(Rpsi) = %.20f\nmax(Rw) = %.20f\n", maxRpsi, maxRw);
        
        //printf("line %d\n", 393);
                
        printf("max(Rw, Rpsi) = %f\n", R);
        //printf("min(Rw, Rpsi) = %f\n", R);
        char buffer[100];
        int cx;
        if(R <= 1.0e-6) {
            stop = 1;
            cx = snprintf(buffer, 100,  "./Results/omega_%d_%0.0f_dt%d.csv", Nel, Re, which_h);
            print_results(w[1], I, J, buffer); //"./Results/omega_64_1000_h2s10.csv");
            cx = snprintf(buffer, 100, "./Results/psi_%d_%0.0f_dt%d.csv", Nel, Re, which_h);
            print_results(psi[1], I, J, buffer); //"./Results/psi_64_1000_h2s10.csv");
        }
        else {
            for(i = 0; i < I; i++) {
                for(j = 0; j < J; j++) {
                    w[0][i][j] = w[1][i][j];
                    psi[0][i][j] = psi[1][i][j];
                    // printf("\t%f", w[0][i][j]);
                }
                // printf("\n");
            }
            zeroArray(w[1], I, J);
            zeroArray(psi[1], I, J);
            zeroArray(what, I, J);
            zeroArray(psihat, I, J);
            n++; // Number of iterations
        }
        
        //printf("line %d\n", 413);
        
    }
    
    printf("\n\n\n******************\n Final n = %d\n", n);
    
    /* Freeing memory */
    free(x);
    free(y);
    for(k = 0; k < 2; k++) {
        destroyArray(w[k]);
        destroyArray(psi[k]);
    }
    destroyArray(what);
    destroyArray(psihat);
    // free(U);
}

int main() {
    cavity();
    printf("\nNel = %d\n", Nel);
    printf("\nReynolds = %f\n\n", Re);
    return 0;
}






