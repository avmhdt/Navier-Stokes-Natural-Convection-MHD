#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include "parameters_boussinesq.h"

double max(double num1, double num2);

double min(double num1, double num2);

// double* seq(double a, double b, double by);

double* seq(double start, double finish, int size);

double max_vec(double vec[], int size);

double min_vec(double vec[], int size);

double max_mat(double mat[][J-2], int n, int m);

double min_mat(double mat[][J-2], int n, int m);

int which_max(double vec[], int size);

int which_min(double vec[], int size);

void print_results(double **u, int rows, int cols, char filename[]);

double** createArray(int n, int m);

void destroyArray(double** arr);

void zeroArray(double** arr, int n, int m);

#endif