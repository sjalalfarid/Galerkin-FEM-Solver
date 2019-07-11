#ifndef FEM_INTEGRALFUNCTIONS

//typedef void (*integral)(double ***terms, double *c, int nodesPerElem, // inputs
//													double **E);	   // outputs

//void SteadyStateDiffusion(double ***terms, double *c, int nodesPerElem, // inputs
//													double **E);	    // outputs

void GeneralSecondOrderEquation(double ***terms, double *c, int nodesPerElem, // inputs
																 double **E);  // outputs

#define FEM_INTEGRALFUNCTIONS
#endif