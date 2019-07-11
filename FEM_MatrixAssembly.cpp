#include "FEM_MatrixAssembly.h"
#include "FEM_BasisFunctions.h"
#include "FEM_SourceTerm.h"

typedef void(*function)(int*, double*, double*, double**);

void My_AssembleLocalElementMatrix(int *nodes, double *x, double *y, double *c, int elementType, int dimension, double **E)
{
		double** terms = NULL;
		terms = new double*[((elementType+1)*dimension)];		//dynamic allocation for 2-D array of Terms
		

		for (int i = 0; i < ((elementType + 1)*dimension); i++)
		{
			terms[i] = new double[((elementType + 1)*dimension)];	//dynamic allocation 2-D array of Terms
		}

	if (elementType==1 && dimension ==1)		// Depending on Element Type and Dimension, we create the required 
											    //Local Element Stifness Matrix
	{ 
		for (int i = 0; i < (elementType + 1); i++)		//	fill 2-D array to create local stiffness matrix
		{
			for (int j = 0; j < (elementType + 1); j++)
			{

				My_Linear1DBasisFunction_uw_11(nodes, x, y, terms);
				double a = terms[i][j];	
				My_Linear1DBasisFunction_wdudx_11(nodes, x, y, terms);
				double b = terms[i][j];						 //Putting every term seperately and then Combining it.
				My_Linear1DBasisFunction_dudxdwdx_11(nodes, x, y, terms);
				double d = terms[i][j];

				E[i][j] = (c[2] * d) - (c[0] * a) - (c[1] * b);
			}
			
		}

		for (int i = 0; i < (elementType + 1); ++i)
			{	
				delete[] terms[i];			// Free up Terms Matrix; Preventing Memory Leaks
			}
		delete[] terms;
	}
	else if (elementType == 2)			// If element is Quadratic
	{
		for (int i = 0; i < (elementType + 1); i++)		//	fill 2-D array to create local stiffness matrix
		{
			for (int j = 0; j < (elementType + 1); j++)
			{

				My_Quadratic1DBasisFunction_uw_11(nodes, x, y, terms);
				double a = terms[i][j];
				My_Quadratic1DBasisFunction_wdudx_11(nodes, x, y, terms);
				double b = terms[i][j];
				My_Quadratic1DBasisFunction_dudxdwdx_11(nodes, x, y, terms);
				double d = terms[i][j];

				E[i][j] = (c[2] * d) - (c[0] * a) - (c[1] * b);
			}

		}

		for (int i = 0; i < (elementType + 1); ++i)
		{
			delete[] terms[i];			// Free up Terms Matrix; Preventing Memory Leaks
		}
		delete[] terms;
	}

		if (elementType == 1 && dimension == 2)			// If Element is 2-D
	{
		for (int i = 0; i < ((elementType + 1)*dimension); i++)		//	fill 2-D array to create local stiffness matrix
		{
			for (int j = 0; j < ((elementType + 1)*dimension); j++)
			{

				My_Linear2DBasisFunction_uw_11(nodes, x, y, terms);
				double a = terms[i][j];

				My_Linear2DBasisFunction_dudxdwdx_11(nodes, x, y, terms);
				double d = terms[i][j];

				E[i][j] = ((c[2] * d) - (c[0] * a));
			}

		}

		for (int i = 0; i < ((elementType + 1)*dimension); ++i)
		{
			delete[] terms[i];				// Free up Terms Matrix; Preventing Memory Leaks
		}
		delete[] terms;

	}
}

void My_AssembleGlobalElementMatrix(int numP, int numE, int nodesPerElem, int **elem, double ***E, double **K)
{
	for (int l = 0; l < numP; l++)
	{
		for (int m = 0; m < numP; m++)		//Making sure all the values in K Matrix is Zero.
			K[l][m] = 0;
	}

	for (int e = 0; e < numE; e++)				//Assemblying the K matrix with the help of element Matrix "elem" to put every 
												//term in its appropriate position.
	{
		for (int i = 0; i < nodesPerElem; i++)
		{
			for (int j = 0; j < nodesPerElem; j++)
			{
				K[elem[e][i]][elem[e][j]] = K[elem[e][i]][elem[e][j]] + E[e][i][j];	
			}
		}
	}
}

void AssembleRightHandSide(int numE, int **elem, double *x, double *f)
{
	for (int i = 0; i < numE + 1; i++)		// Initializing the f array with 0's to disposed of any
											//garbage values.
	{
		f[i] = 0;
	}

	double checkfunc1;					// A check function is added to determine if the terms of f(x)=0 or not.
	double checkfunc11;					// It is important to check that both the fuction terms are zero at all the
										// nodes to be sure that there is no source term.

	double *check = NULL;	 // The check array store the values of both the terms of f(x) to determine if f(x)=0.
	check = new double[numE + 1];

	for (int i = 0; i < numE + 1; i++)	 // Initializing the check array with 0's to disposed of any
											//garbage values. 										
	{
		check[i] = 0;
	}

	double *assemble = NULL;			// If f(x) is not equal to zero, then the assemble array will be used to get
										// values of f(x) at different global node positions.
	assemble = new double[numE + 1];

	int l = 0;
	double length = x[l + 1] - x[l];	// The length is simply the total length of each element used for computation 
										// of f(x) at different global node positions.

	for (int e = 0; e < numE; e++)		//check if f(x)=0 for all elemental nodes.
	{
		for (int i = 0; i < 1; i++)
		{
			checkfunc1 = Eval_SourceTerm(x[elem[e][i]]);	//Storing the values of terms of f(x) at the respecitve two nodes using the Source Term
																// Function (Depending on what is defined in cpp file.)
			checkfunc11 = Eval_SourceTerm(x[elem[e][i + 1]]);
		}

		check[e] = checkfunc1 + checkfunc11;	//Putting the terms in check array to make sure both the terms are zero.

		if (check[e] != 0)	//If any value of the check array is non-zero this mean f(x) is not zero.
			goto assembly;	//Therefore, Source Term exists and the routine is directed to assemble the f array.
	}

	goto cout;	 //Otherwise the program displays the message that f(x)=0.

assembly:				//  *************** When F(x) is not zero *******************

	for (int i = 0; i < numE + 1; i++)			// Initializing the assemble array with 0's to disposed of any
											//garbage values.
	{

		assemble[i] = 0;
	}

	for (int e = 0; e < numE; e++)	//For All elements we need to find the local load vector and then assemble it to  
									//find the global vector.
	{
		assemble[e] = f[l];		// The first entry of local load vector is stored so that it can be added with the
						// respective first entry of the connecting node (load vector) to calculate the global load vector.

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 1; j++)
			{
				double func1 = Eval_SourceTerm(x[elem[e][j]]);		// The two function terms are calculated at the two nodes of the element with the help of Eval function.
				double func11 = Eval_SourceTerm(x[elem[e][j + 1]]);

				if (i == 0)		// The value is added to the respective position in the f vector
				{				//****** The second term with "0" in it is for the bubble function addition if Consistent SUPG Method is desired.
								//****** Replace the "0" in both eq. for any value of Beta to incorporate for the bubble function.
					f[i] = -((((func1*length)/3)-((func1)*((0*length)/4))) + (((func11*length)/6)+ ((func11)*((0*length)/4))));
				}						
				else if (i == 1)
				{
					f[i] = -((((func1*length)/6) - ((func1)*((0*length)/4))) + (((func11*length) / 3) + ((func11)*((0*length)/4))));
				}
			}
		}
		assemble[e] = assemble[e] + f[l];	//The values are then assembled to the "assemble array" to create global load vector.
		assemble[e + 1] = f[l + 1];
	}
	for (int e = 0; e < numE + 1; e++)
	{
		f[e] = assemble[e];			//All the global values of assemble array are then transfered 
									//to the f vector which is then returned by the function for the 
									// application of essential and natural boundary conditions.
	}

	delete[] check ;				// Free up Terms Matrix; Preventing Memory Leaks
	delete[] assemble;				// Free up Terms Matrix; Preventing Memory Leaks

	return;

cout:			// If f(x)=0 the routine will directed to this statement.
	
	std::cout << "No Source Term Detected (f(x)=0)" << "\n";

	delete[] check;				// Free up Terms Matrix; Preventing Memory Leaks
	delete[] assemble;			// Free up Terms Matrix; Preventing Memory Leaks

}

void My_ApplyNaturalBoundaryConditions(int numP, int numBCn, int *BCnnode, double *BCnvalue, double *c, double *f)
{
	for (int i = 0; i < numBCn; i++)			 
		//Adding the respective Natural Boundary Conditions to the first and last term of the RHS (Global f vector). 						
	{	
		if (BCnnode[i] == 0)    
			f[BCnnode[i]] = f[BCnnode[i]] + (-c[2] * BCnvalue[i]);	
		else if (BCnnode[i] == (numP - 1))
			f[BCnnode[i]] = f[BCnnode[i]] + (c[2] * BCnvalue[i]);
	}
}

void My_ApplyEssentialBoundaryConditions(int numP, int numBCe, int *BCenode, double *BCevalue, double *f, double **K)
{
		for (int i = 0; i < numBCe; i++)		 
			//Adding the respective Essential Boundary Conditions at the respective nodes in the RHS (Global f vector). 										
		{
			f[BCenode[i]] = BCevalue[i];
		}
		

		for (int i = 0; i < numP; i++)	
		{
			for (int c = 0; c < numBCe; c++)
			{
				if (BCevalue[c] == f[i])	// This serves as a check on which nodes the BC's is applied.
				{
					for (int j = 0; j < numP-1; j++)
					{
						K[BCenode[c]][BCenode[c]] = 1;	//Decoupling the respective row from global stifness matrix. 
						K[BCenode[c]][j] = 0;
					}
				}
				else
					continue;
			}
		}
}


