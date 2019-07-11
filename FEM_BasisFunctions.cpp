#include "FEM_BasisFunctions.h"

//
// ------------------- Code for solving using linear Lagrange 1D elements -----------------------------
//	Code for for -1 < xi < 1 
//

double My_Linear1DJacobianDet_11(int *nodes, double *x, double *y) 
{
	return ((x[nodes[1]] - x[nodes[0]])/2);
}
void My_Linear1DBasisFunction_uw_11(int *nodes, double *x, double *y, double **term)
{
	double J = My_Linear1DJacobianDet_11(nodes, x, y);
	term[0][0] = term[1][1] = ((2.0*J) / 3);
	term[0][1] = term[1][0] = (J / 3);
}
void My_Linear1DBasisFunction_wdudx_11(int *nodes, double *x, double *y, double **term)
{
	term[0][0] = -0.5 + 0;	// ****** The second term after the " + / - " sign is the bubble function value (Beta/2).
							// ****** If it is zero it means no bubble function is added and standard GFEM Method is used.
							// ****** If it is not zero it means bubble function is added and Streamline Upwind Petrov GFEM Method is used.
	term[0][1] = 0.5 - 0;
	term[1][0] = -0.5 - 0;
	term[1][1] = 0.5 + 0;
}
void My_Linear1DBasisFunction_dudxdwdx_11(int *nodes, double *x, double *y, double **term)
{
	double J = My_Linear1DJacobianDet_11(nodes, x, y);
	term[0][0] = term[1][1] = (0.5) / J;
	term[0][1] = term[1][0] = (-0.5) / J;
}

// ------------------- Code for solving using  1D Quadratic element -----------------------------
double My_Quadratic1DJacobianDet_11(int *nodes, double *x, double *y)
{
	return ((x[nodes[1]] - x[nodes[0]]) / 1);
}

void My_Quadratic1DBasisFunction_uw_11(int *nodes, double *x, double *y, double **term)
{
	double J = My_Quadratic1DJacobianDet_11(nodes, x, y);
	term[0][0] = (4.0*J) / 15;
	term[0][1] = term[1][0] = (2.0*J) / 15;
	term[0][2] = term[2][0] = (-J) / 15;
	term[1][1] = (16.0*J) / 15;
	term[1][2] = term[2][1] = (2.0*J) / 15;
	term[2][2] = (4.0*J) / 15;
}

void My_Quadratic1DBasisFunction_wdudx_11(int *nodes, double *x, double *y, double **term)
{
	term[0][0] = -0.5;
	term[0][1] = (2.0) / (3);
	term[0][2] = (-1.0) / (6);
	term[1][0]= (-2.0) / (3);
	term[1][1] = 0;
	term[1][2] =  (2.0) / (3);
	term[2][0] = (1.0) / (6);
	term[2][1] = (-2.0) / (3);
	term[2][2] = 0.5;
}

void My_Quadratic1DBasisFunction_dudxdwdx_11(int *nodes, double *x, double *y, double **term)
{
	double J = My_Quadratic1DJacobianDet_11(nodes, x, y);
	term[0][0] =  (7.0 )/(6*J) ;
	term[0][1] = term [1][0] = (-4.0)/(3 * J);
	term[0][2] = term [2][0] = (1.0)/(6 * J);
	term[1][1] = (8.0)/(3 * J);
	term[1][2] = term[2][1] = (-4.0)/(3 * J);
	term[2][2] = (7.0) / (6 * J);
}
// ------------------- Code for solving using  2D Linear element -----------------------------
double My_Linear2DJacobianDet_11(int *nodes, double *x, double *y)
{
	return ((x[nodes[1]] - x[nodes[0]])*(y[nodes[0]] - y[nodes[2]]) / 4);
}
void My_Linear2DBasisFunction_uw_11(int *nodes, double *x, double *y, double **term)
{
		
	double delx = x[nodes[1]] - x[nodes[0]];	// Evaluating the terms in terms of delx, dely.
	double dely = y[nodes[0]] - y[nodes[2]];

	term[0][0] = (delx*dely) / 9;
	term[0][1] = term[1][0] = (delx*dely) / 18;
	term[0][2] = term[2][0] = (delx*dely) / 18;
	term[0][3] = term[3][0] = (delx*dely) / 36;

	term[1][1] = (delx*dely) / 9;
	term[1][2] = term[2][1] = (delx*dely) / 36;
	term[1][3] = term[3][1] = (delx*dely) / 18;

	term[2][2] = (delx * dely) / 9;
	term[2][3] = term[3][2] = (delx*dely) / 18;

	term[3][3] = (delx*dely) / 9;					
}
						//Imp!! no need to calculate c1 terms for 2-D as heat conduction equation has c1 = 0.

void My_Linear2DBasisFunction_dudxdwdx_11(int *nodes, double *x, double *y, double **term)
{
	double delx = x[nodes[1]] - x[nodes[0]];		// Evaluating the terms in terms of delx, dely.
	double dely = y[nodes[0]] - y[nodes[2]];

	term[0][0] = ((delx/dely) + (dely/ delx)) / 3;
	term[0][1] = term[1][0] = ((delx/(6*dely))-(dely/(3* delx)));
	term[0][2] = term[2][0] = ((dely / (6 * delx)) - (delx / (3 * dely)));
	term[0][3] = term[3][0] = (((-dely) / (6 * delx)) - (delx / (6 * dely)));

	term[1][1] = ((delx / (3 * dely)) + (dely / (3 * delx)));
	term[1][2] = term[2][1] = (((-dely) / (6 * delx)) - (delx / (6 * dely)));
	term[1][3] = term[3][1] = (((-delx) / (3 * dely)) + (dely / (6 * delx)));

	term[2][2] = ((delx / (3 * dely)) + (dely / (3 * delx)));
	term[2][3] = term[3][2] = (((-dely) / (3 * delx)) + (delx / (6 * dely)));

	term[3][3] = ((dely / (3 * delx)) + (delx / (3 * dely)));
}