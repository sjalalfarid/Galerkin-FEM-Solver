#ifndef FEM_MATRIXASSEMBLY

#include <stdlib.h>
#include <iostream>

using namespace std;

#include "FEM_BasisFunctions.h"
#include "FEM_Misc.h"

//--------------------------------------------------------------------------------------------------------

void My_AssembleLocalElementMatrix(int *nodes, double *x, double *y, double *c, int elementType, int dimension, double **E);
void My_AssembleGlobalElementMatrix(int numP, int numE, int nodesPerElem, int **elem, double ***E, double **K);
void My_ApplyEssentialBoundaryConditions(int numP, int numBCe, int *BCenode, double *BCevalue, double *f, double **K);
void My_ApplyNaturalBoundaryConditions(int numP, int numBCn, int *BCnnode, double *BCnvalue, double *c, double *f);
void AssembleRightHandSide(int numE, int **elem, double *x, double *f); //Added to Evaluate the source Term F(x).
#define FEM_MATRIXASSEMBLY
#endif