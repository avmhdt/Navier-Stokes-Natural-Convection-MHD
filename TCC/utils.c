#include "utils.h"

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