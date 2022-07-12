#ifndef PARAMETERS_MHD_H
#define PARAMETERS_MHD_H

#define A 0.0
#define B 1.0

#define Nel 128
#define I (Nel+1)
#define J I

#define h (((double)B - (double)A)/(double)Nel)

#define Pr (double) 1.0
#define Ra (double) 1.e6
#define Ha (double) 0.0//100.0

#define alpha_psi (double) 0.1//1.0
#define alpha_w (double) 0.01//1.0
#define alpha_T (double) 0.5//1.0

#define N 0.5

#define which_h 1

#endif