#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils_mhd.h"
#include "myTridiag.h"
#include "parameters_mhd.h"

#define PI 3.14159265358979323846

void mhd() {

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

    double h2 = pow(h, 2.);
    double sig;
    sig = dt/(2.*h2);
    
    double **w, **psi, **T;
    double **what, **psihat, **That;
    int n; // time counter
    int i, j; // x and y dimension counters
    
    int k; // other counter

    w = createArray(I, J);
    psi = createArray(I, J);
    T = createArray(I, J);
    what = createArray(I, J);
    psihat = createArray(I, J);
    That = createArray(I, J);
    
    /* Vectors from AU = F linear system of equations from ADI method */
    double F[J-2] = {0.};           
    double Di[J-2] = {0.}, Dp[J-2] = {0.}, Ds[J-2] = {0.};   // Diagonals from A matrix
    double *U = NULL;

    // RhoX and RhoY:
    double RhoX[I][J] = {{0.}}, RhoY[I][J] = {{0.}};
    
    // For Temperature:
    double F1[I] = {0.};
    double Di1[I] = {0.}, Dp1[I] = {0.}, Ds1[I] = {0.};   // Diagonals from A matrix
    double *U1 = NULL;

    // Theta for vorticity:
    double theta[I-2][J-2] = {{0.}};

    // Psi_x^2 for vorticity:
    double PSIx2[I-2][J-2] = {{0.}};
    
    // Residuals:
    double Rpsi[I-2][J-2] = {{0.}}, Rw[I-2][J-2] = {{0.}}, RT[I-2][J-2] = {{0.}};
    double R, maxRpsi, maxRw, maxRT;
    
    // Initial conditions:
    zeroArray(w, I, J);
    zeroArray(psi, I, J);
    zeroArray(T, I, J);
    zeroArray(what, I, J);
    zeroArray(psihat, I, J);
    zeroArray(That, I, J);

    int stop = 0; // boolean to check stop criteria
    n = 1;
    while(!stop) {
        printf("**************** n = %d\n", n);
        // Stream function:
        // Implicit x:
        for(j = 1; j < J-1; j++) {
            for(i = 1; i < I-1; i++) {
                Di[i-1] = -sig;
                Dp[i-1] = (1./alpha_psi + 2.*sig);
                Ds[i-1] = -sig;

                F[i-1] = dt*w[i][j]/2. + (1./alpha_psi - 2.*sig)*psi[i][j] + sig*(psi[i][j+1] + psi[i][j-1]);
            }

            U = myTridiag(Di, Dp, Ds, F, I-2);
            for(i = 1; i < I-1; i++) {
                psihat[i][j] = U[i-1];
            }
            free(U);
            U = NULL;
        }
        //Implicit y:
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                Di[j-1] = -sig;
                Dp[j-1] = (1./alpha_psi + 2.*sig);
                Ds[j-1] = -sig;

                F[j-1] = dt*w[i][j]/2. + (1./alpha_psi - 2.*sig)*psihat[i][j] + sig*(psihat[i+1][j] + psihat[i-1][j]);
            }

            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                psi[i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }

        // Temperature:
        // Boundary conditions:
        for(i = 0; i < I; i++) {
            That[i][0] = N*sin(2.*PI*x[i]); // bottom
            That[i][J-1] = 0.;              // top
            T[i][0] = N*sin(2.*PI*x[i]);    // bottom
            T[i][J-1] = 0.;                 // top
        }

        for(i = 1; i < I-1; i++) {
            for(j = 0; j < J; j++) {
                RhoX[i][j] = (psi[i+1][j] - psi[i-1][j])/4.;
                RhoY[j][i] = (psi[j][i+1] - psi[j][i-1])/4.;
            }
        }
        for(j = 0; j < J; j++) {
            RhoX[0][j] = RhoX[1][j];
            RhoX[I-1][j] = RhoX[I-2][j];
        }
        for(i = 0; i < I; i++) {
            RhoY[i][0] = RhoY[i][1];
            RhoY[i][J-1] = RhoY[i][J-2];
        }
        // Implicit x:
        for(j = 1; j < J-1; j++) {
            for(i = 0; i < I; i++) {
                Di1[i] = -sig*(RhoY[i][j] + 1.);
                Dp1[i] = (1./alpha_T + 2.*sig);
                Ds1[i] = sig*(RhoY[i][j] - 1.);

                F1[i] = (1./alpha_T - 2.*sig)*T[i][j] + sig*(RhoX[i][j] + 1.)*T[i][j+1] + sig*(1. - RhoX[i][j])*T[i][j-1];
            }
            
            Ds1[0] = Ds1[0] - sig*(RhoY[0][j] + 1.);
            Di1[I-1] = Di1[I-1] + sig*(RhoY[I-1][j] - 1.);

            U1 = myTridiag(Di1, Dp1, Ds1, F1, I);
            for(i = 0; i < I; i++) {
                That[i][j] = U1[i];
            }
            free(U1);
            U1 = NULL;
        }
        // Implicit y:
        for(i = 0; i < I; i++) {
            for(j = 1; j < J-1; j++) {
                Di[j-1] = sig*(RhoX[i][j] - 1.);
                Dp[j-1] = (1./alpha_T + 2.*sig);
                Ds[j-1] = -sig*(RhoX[i][j] + 1.);

                if (i == 0) {
                    F[j-1] = (1./alpha_T - 2.*sig)*That[i][j] + sig*(1. - RhoY[i][j])*That[i+1][j] + sig*(1. + RhoY[i][j])*That[i+1][j];
                } else if (i == I-1) {
                    F[j-1] = (1./alpha_T - 2.*sig)*That[i][j] + sig*(1. - RhoY[i][j])*That[i-1][j] + sig*(1. + RhoY[i][j])*That[i-1][j];
                } else {
                    F[j-1] = (1./alpha_T - 2.*sig)*That[i][j] + sig*(1. - RhoY[i][j])*That[i+1][j] + sig*(1. + RhoY[i][j])*That[i-1][j];
                }
            }

            F[0] = F[0] - sig*(RhoX[i][1] - 1.)*T[i][0];
            F[J-3] = F[J-3] + sig*(RhoX[i][J-2] + 1.)*T[i][J-1];

            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                T[i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }

        // Vorticity:
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                theta[i-1][j-1] = dt*Ra*Pr*(T[i+1][j] - T[i-1][j])/(4.*h);
                PSIx2[i-1][j-1] = dt*pow(Ha, 2.)*Pr*(psi[i+1][j] - 2.*psi[i][j] + psi[i-1][j])/(2.*h2);
            }
        }
        // Boundary conditions:
        for(j = 0; j < J; j++) {
            what[0][j] = -2.*psi[1][j]/h2; // left
            what[I-1][j] = -2.*psi[I-2][j]/h2; // right
        }
        for(i = 0; i < I; i++) {
            what[i][0] = -2.*psi[i][1]/h2; // bottom
            what[i][J-1] = -2.*psi[i][J-2]/h2; // top
        }
        // Implicit x:
        for(j = 1; j < J-1; j++) {
            for(i = 1; i < I-1; i++) {
                Di[i-1] = -sig*(RhoY[i][j] + Pr);
                Dp[i-1] = (1./alpha_w + 2.*Pr*sig);
                Ds[i-1] = sig*(RhoY[i][j] - Pr);

                F[i-1] = theta[i-1][j-1] + PSIx2[i-1][j-1] + (1./alpha_w - 2.*Pr*sig)*w[i][j] + sig*(RhoX[i][j] + Pr)*w[i][j+1] + sig*(Pr - RhoX[i][j])*w[i][j-1];
            }

            F[0] = F[0] + sig*(RhoY[1][j] + Pr)*what[0][j];
            F[I-3] = F[I-3] - sig*(RhoY[I-2][j] - Pr)*what[I-1][j];

            U = myTridiag(Di, Dp, Ds, F, I-2);
            for(i = 1; i < I-1; i++) {
                what[i][j] = U[i-1];
            }
            free(U);
            U = NULL;
        }
        // Boundary conditions:
        for(j = 0; j < J; j++) {
            w[0][j] = what[0][j]; // left
            w[I-1][j] = what[I-1][j]; // right
        }
        for(i = 0; i < I; i++) {
            w[i][0] = what[i][0]; // bottom
            w[i][J-1] = what[i][J-1]; // top
        }
        // Implicit y:
        for(i = 1; i < I-1; i++) { 
            for(j = 1; j < J-1; j++) {
                Di[j-1] = sig*(RhoX[i][j] - Pr);
                Dp[j-1] = (1./alpha_w + 2.*Pr*sig);
                Ds[j-1] = -sig*(RhoX[i][j] + Pr);

                F[j-1] = theta[i-1][j-1] + PSIx2[i-1][j-1] + (1./alpha_w - 2.*Pr*sig)*what[i][j] + sig*(Pr - RhoY[i][j])*what[i+1][j] + sig*(Pr + RhoY[i][j])*what[i-1][j];
            }

            F[0] = F[0] - sig*(RhoX[i][1] - Pr)*w[i][0];
            F[J-3] = F[J-3] + sig*(RhoX[i][J-2] + Pr)*w[i][J-1];

            U = myTridiag(Di, Dp, Ds, F, J-2);
            for(j = 1; j < J-1; j++) {
                w[i][j] = U[j-1];
            }
            free(U);
            U = NULL;
        }

        /* Calculating approximation residuals */
        for(i = 1; i < I-1; i++) {
            for(j = 1; j < J-1; j++) {
                /*
                Rpsi[i-1][j-1] = fabs((psi[i+1][j] - 2.*psi[i][j] + psi[i-1][j] + psi[i][j+1] - 2.*psi[i][j] + psi[i][j-1])/h2 + w[i][j]);
                Rw[i-1][j-1] = fabs(Pr*(w[i+1][j] - 2.*w[i][j] + w[i-1][j] + w[i][j+1] - 2.*w[i][j] + w[i][j-1])/(Ra*h2)
                                    - (psi[i][j+1] - psi[i][j-1])*(w[i+1][j] - w[i-1][j])/(Ra*4.*h2)
                                    + (psi[i+1][j] - psi[i-1][j])*(w[i][j+1] - w[i][j-1])/(Ra*4.*h2)
                                    + Pr*(T[i+1][j] - T[i-1][j])/(2.*h));
                RT[i-1][j-1] = fabs((T[i+1][j] - 2.*T[i][j] + T[i-1][j] + T[i][j+1] - 2.*T[i][j] + T[i][j-1])/h2
                                     - (psi[i][j+1] - psi[i][j-1])*(T[i+1][j] - T[i-1][j])/(4.*h2)
                                      + (psi[i+1][j] - psi[i-1][j])*(T[i][j+1] - T[i][j-1])/(4.*h2));
                */
                Rpsi[i-1][j-1] = fabs((psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - 4.*psi[i][j])/h2 + w[i][j]);
                Rw[i-1][j-1] = fabs(Pr*(w[i+1][j] + w[i-1][j] + w[i][j+1] + w[i][j-1] - 4.*w[i][j])/(Ra*h2)
                                    - ((psi[i][j+1] - psi[i][j-1])*(w[i+1][j] - w[i-1][j])
                                    - (psi[i+1][j] - psi[i-1][j])*(w[i][j+1] - w[i][j-1]))/(Ra*4.*h2)
                                    + Pr*(T[i+1][j] - T[i-1][j])/(2.*h)
                                    + pow(Ha, 2.)*Pr*(psi[i+1][j] - 2.*psi[i][j] + psi[i-1][j])/(Ra*h2));
                RT[i-1][j-1] = fabs((T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4.*T[i][j])/h2
                                     - ((psi[i][j+1] - psi[i][j-1])*(T[i+1][j] - T[i-1][j])
                                      - (psi[i+1][j] - psi[i-1][j])*(T[i][j+1] - T[i][j-1]))/(4.*h2));
            }
        }
        maxRpsi = max_mat(Rpsi, I-2, J-2);
        maxRw = max_mat(Rw, I-2, J-2);
        maxRT = max_mat(RT, I-2, J-2);

        R = max(maxRpsi, maxRw);
        R = max(R, maxRT);

        printf("max(Rpsi) = %.20f\nmax(Rw) = %.20f\nmax(RT) = %.20f\n", maxRpsi, maxRw, maxRT);
        printf("max(Rw, Rpsi) = %f\n", R);

        char buffer[100];
        int cx;
        if(R <= 1.0e-6) {
            stop = 1;
            cx = snprintf(buffer, 100,  "./Results/MHD/omega_%d_%0.2f_%0.0f_%0.0f_%0.1f_%0.2f_%0.2f_%0.2f_dt%d_%d.csv", Nel, Pr, Ra, Ha, N, alpha_psi, alpha_w, alpha_T, which_h, n);
            print_results(w, I, J, buffer);
            cx = snprintf(buffer, 100, "./Results/MHD/psi_%d_%0.2f_%0.0f_%0.0f_%0.1f_%0.2f_%0.2f_%0.2f_dt%d_%d.csv", Nel, Pr, Ra, Ha, N, alpha_psi, alpha_w, alpha_T, which_h, n);
            print_results(psi, I, J, buffer);
            cx = snprintf(buffer, 100, "./Results/MHD/T_%d_%0.2f_%0.0f_%0.0f_%0.1f_%0.2f_%0.2f_%0.2f_dt%d_%d.csv", Nel, Pr, Ra, Ha, N, alpha_psi, alpha_w, alpha_T, which_h, n);
            print_results(T, I, J, buffer);
        }
        else {
            n++; // Number of iterations
        }
    }

    printf("\n\n\n******************\n Final n = %d\n", n);
    
    /* Freeing memory */
    free(x);
    free(y);
    destroyArray(w);
    destroyArray(psi);
    destroyArray(T);
    destroyArray(what);
    destroyArray(psihat);
    destroyArray(That);
}

int main() {
    mhd();
    printf("\nNel = %d\n", Nel);
    printf("\nPrandtl = %f\n\n", Pr);
    printf("\nRayleigh = %f\n\n", Ra);
    printf("\nHartmann = %f\n\n", Ha);
    return 0;
}
