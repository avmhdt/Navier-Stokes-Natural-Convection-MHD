#ifndef PARAMETERS_BOUSSINESQ_H
#define PARAMETERS_BOUSSINESQ_H

#define A 0.0
#define B 1.0

#define Nel 80
#define I (Nel+1)
#define J I

#define h (((double)B - (double)A)/(double)Nel)

#define Pr (double) 0.71
#define Ra (double) 1.e6

#define alpha_psi (double) 0.1
#define alpha_w (double) 0.01
#define alpha_T (double) 0.5

#define which_h 1

#endif