#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 0.0
#define B 1.0
#define Nel 256
#define h (((double)B - (double)A)/(double)Nel)

#define I Nel+1
#define J I

#define Re (double) 1000.0

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
    sig = dt/(2.*pow(h, 2));

    double **w, **psi;
    double **what, **psihat;
    int n; // time counter
    int i, j; // x and y dimension counters
    
    int k; // other counter
    
    w = createArray(I, J);
    psi = createArray(I, J);
    what = createArray(I, J);
    psihat = createArray(I, J);
    
    /* Vectors from AU = F linear system of equations from ADI method */
    double F[J-2] = {0};           
    double Di[J-2] = {0}, Dp[J-2] = {0}, Ds[J-2] = {0};   // Diagonals from A matrix
    double *U = NULL;
    // RhoX and RhoY:
    double RhoX[I-2][J-2] = {{0}}, RhoY[I-2][J-2] = {{0}};
    
    // Residuals:
    double Rpsi[I-2][J-2] = {{0}}, Rw[I-2][J-2] = {{0}};
    double R, maxRpsi, maxRw;
    
    // Initial conditions:
    zeroArray(w, I, J);
    zeroArray(psi, I, J);
    zeroArray(what, I, J);
    zeroArray(psihat, I, J);
    
    int stop = 0; // boolean to check stop criteria
    n = 1;
    
    while(!stop) {
        printf("**************** n = %d\n", n);
        
        /*************** ADI method ***************/
        
        //printf("line %d\n", 262);
        
        // Filling diagonals from A matrix first
        for(k = 0; k < J-2; k++) {
            Di[k] = -sig;
            Dp[k] = 1. + 2.*sig;
            Ds[k] = -sig;  
        }
        
        //printf("line %d\n", 271);
        
        /******* Step 1: Implicit x *******/
        for(j = 1; j < J-1; j++) {
            // Filling source term F:
            for(i = 1; i < I-1; i++) {
                F[i-1] = (dt/2.)*w[i][j] + (1.-2.*sig)*psi[i][j]
                       + sig*(psi[i][j+1] + psi[i][j-1]);
            }
            
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
                F[j-1] = (dt/2.)*w[i][j] + (1.-2.*sig)*psihat[i][j]
                       + sig*(psihat[i+1][j] + psihat[i-1][j]);
            }
            
            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                psi[i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 309);
        
        for(j = 0; j < J; j++) {
            what[0][j] = -2.*psi[1][j]/pow(h, 2);             // Left
            what[I-1][j] = -2.*psi[I-2][j]/pow(h, 2);         // Right
            what[j][0] = -2.*psi[j][1]/pow(h, 2);             // Inferior
            what[j][J-1] = -2.*psi[j][J-2]/pow(h, 2) - 2.0/h;   // Superior
        }

        //printf("line %d\n", 326);
        for(j = 1; j < J-1; j++) {
            for(i = 1; i < I-1; i++) {
                RhoX[i-1][j-1] = Re*(psi[i+1][j] - psi[i-1][j])/4.;
                RhoY[i-1][j-1] = Re*(psi[i][j+1] - psi[i][j-1])/4.;
            }
        }
        /******* Step 1: Implicit x *******/
        for(j = 1; j < J-1; j++) {
            for(i = 1; i < I-1; i++) {
                
                Di[i-1] = -sig*(RhoY[i-1][j-1] + 1.);
                Dp[i-1] = 1. + 2.*sig;
                Ds[i-1] = sig*(RhoY[i-1][j-1] - 1.);
                
                F[i-1] = (1.-2.*sig)*w[i][j] + sig*(RhoX[i-1][j-1] + 1.)*w[i][j+1]
                                                - sig*(RhoX[i-1][j-1] - 1.)*w[i][j-1];
            }
            
            F[0] = F[0] + sig*(RhoY[0][j-1] + 1.)*what[0][j];
            F[I-3] = F[I-3] - sig*(RhoY[I-3][j-1] - 1.)*what[I-1][j];
            
            U = myTridiag(Di, Dp, Ds, F, I-2);
            for(i = 1; i < I-1; i++) {
                what[i][j] = U[i-1];

            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 352);
        
        for(j = 0; j < J; j++) {
            w[0][j] = -2.*psi[1][j]/pow(h, 2);             // Left
            w[I-1][j] = -2.*psi[I-2][j]/pow(h, 2);         // Right
            w[j][0] = -2.*psi[j][1]/pow(h, 2);             // Inferior
            w[j][J-1] = -2.*psi[j][J-2]/pow(h, 2) - 2.0/h;   // Superior
        }

        /******* Step 2: Implicit y *******/
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                
                Di[j-1] = sig*(RhoX[i-1][j-1] - 1.);
                Dp[j-1] = 1. + 2.*sig;
                Ds[j-1] = -sig*(RhoX[i-1][j-1] + 1.);
                
                F[j-1] = (1.-2.*sig)*what[i][j] - sig*(RhoY[i-1][j-1] - 1.)*what[i+1][j]
                                                + sig*(RhoY[i-1][j-1] + 1.)*what[i-1][j];
            }
            
            F[0] = F[0] - sig*(RhoX[i-1][0] - 1.)*w[i][0];
            F[J-3] = F[J-3] + sig*(RhoX[i-1][J-3] + 1.)*w[i][J-1];
            
            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                w[i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }
        
        //printf("line %d\n", 378);
        
        /* Calculating approximation residuals */
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                Rpsi[i-1][j-1] = fabs((psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - 4.*psi[i][j])/pow(h, 2) + w[i][j]);
                
                Rw[i-1][j-1] = fabs((w[i+1][j] + w[i-1][j] + w[i][j+1] + w[i][j-1] - 4.*w[i][j])/(Re*pow(h, 2)) - ((psi[i][j+1] - psi[i][j-1])*(w[i+1][j] - w[i-1][j]) - (psi[i+1][j] - psi[i-1][j])*(w[i][j+1] - w[i][j-1]))/(4.*pow(h, 2)));
                
            }
        }
        
        //printf("line %d\n", 389);
        
        maxRpsi = max_mat(Rpsi, I-2, J-2);
        maxRw = max_mat(Rw, I-2, J-2);
        R = fmax(maxRpsi, maxRw);
        //R = min(maxRpsi, maxRw);
//        printf("max(Rpsi) = %.20f\nmax(Rw) = %.20f\n", maxRpsi, maxRw);
        
        //printf("line %d\n", 393);
                
        printf("max(Rw, Rpsi) = %.20f\n", R);
        //printf("min(Rw, Rpsi) = %f\n", R);
        char buffer[100];
        int cx;
        if(R <= 1.0e-6) {
            stop = 1;
            cx = snprintf(buffer, 100,  "./Results/omega_%d_%0.0f_dt%d.txt", Nel, Re, which_h);
            print_results(w, I, J, buffer); //"./Results/omega_64_1000_h2s10.csv");
            cx = snprintf(buffer, 100, "./Results/psi_%d_%0.0f_dt%d.txt", Nel, Re, which_h);
            print_results(psi, I, J, buffer); //"./Results/psi_64_1000_h2s10.csv");
        }
        else {
            n++; // Number of iterations
        }
        
        //printf("line %d\n", 413);
        
    }
    
    printf("\n\n\n******************\n Final n = %d\n", n);
    
    /* Freeing memory */
    free(x);
    free(y);
    destroyArray(w);
    destroyArray(psi);
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






