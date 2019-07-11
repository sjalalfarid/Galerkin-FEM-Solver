#ifndef FEM_BASISFUNCTIONS

#include <stdlib.h>
#include <iostream>
#include "FEM_Misc.h"

using namespace std;

// Create a general function pointer variable type 
typedef void (*function)(int* ,double*, double*, double**);				
typedef double (*jacobian)(int*, double*, double*);

// -------------------------------- TASK 1: 1D Linear element ------------------------------------------------
void My_Linear1DBasisFunction_uw_11(int *nodesInElem, double *x, double *y, double **term);
void My_Linear1DBasisFunction_wdudx_11(int *nodesInElem, double *x, double *y, double **term);
void My_Linear1DBasisFunction_dudxdwdx_11(int *nodesInElem, double *x, double *y, double **term);
double My_Linear1DJacobianDet_11(int *nodesInElem, double *x, double *y);

// --------------------------------- TASK 7: 1D Quadratic element --------------------------------------------
void My_Quadratic1DBasisFunction_uw_11(int *nodes, double *x, double *y, double **term);
void My_Quadratic1DBasisFunction_wdudx_11(int *nodes, double *x, double *y, double **term);
void My_Quadratic1DBasisFunction_dudxdwdx_11(int *nodes, double *x, double *y, double **term);
double My_Quadratic1DJacobianDet_11(int *nodes, double *x, double *y);

// --------------------------------- TASK 9: 2D Bi-linear element --------------------------------------------
double My_Linear2DJacobianDet_11(int *nodes, double *x, double *y);
void My_Linear2DBasisFunction_uw_11(int *nodes, double *x, double *y, double **term);
void My_Linear2DBasisFunction_wdudx_11(int *nodes, double *x, double *y, double **term);
void My_Linear2DBasisFunction_dudxdwdx_11(int *nodes, double *x, double *y, double **term);


#define FEM_BASISFUNCTIONS
#endif
