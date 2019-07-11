/*

	2-D Steady State Heat Conduction FEM Solver using Galerkin FEM Scheme with SUPG Scheme and Bubble Functions Capability

*************************** Student Name: Syed Jalal Farid *************************** 
***************************			UoA ID: 814903243      *************************** 

	Building the program requires the following files:
		- FEM_CommandLineFilter.cpp/.h
		- FEM_BasisFunctions.cpp/.h
		- FEM_MatrixAssembly.cpp/.h
		- FEM_MatrixOperations.cpp/.h
		- FEM_SourceTerm.cpp/.h

	Running the program requires that you supply:
	  - FEM_config.txt, a plain text file containing the program running options.  The flags in this file are:
				-n nodefilename -> this is where you specify the number of nodes and their global coordinates
				-e elementfilename -> this is where you specify the number of elements, the number of nodes per elem, and the global connectivity
				-b bcfilename -> number of Dirichlet nodes, the nodes and the values associated with each
				-s elementType -> an example is given for a linear 2node element for elementType = 1.  You will write others.
				-d number of dimensions
				-x if a 2D problem with a rectangular domain, the number of points in the x direction
				-y if a 2D problem with a rectangular domain, the number of points in the y direction


// **************************************************  IMP NOTES !!!!  **************************************************

//***** The SourceTerm.cpp must return a "0" value for homogenous differential equations or non-zero value for non-homogenous differential equations.**************
//***** Make Sure to add the Math.h Header in FEM_SourceTerm.cpp in order to define generic function using C++ Builtin Functions.**************					

//***** The Bubble function feature have been disabled by default for current solver settings but If you want to use the SUPG / Inconsistent Streamline upwind Scheme
		just add the respective bubble functions terms (Beta Values instead of "0") in Basis Function.cpp (In the Linear1DBasisfunction) & Matrix Assembly.cpp (In the Assemble func. for f(x)).  **************

*/

#include <iostream>
#include "FEM_BasisFunctions.h"
#include "FEM_MatrixAssembly.h"
#include "FEM_SourceTerm.h"
#include "FEM_CommandLineFilter.h"
#include "FEM_MatrixOperations.h"
#include "FEM_Misc.h"

using namespace std; 

int main(int argc, char **argv)
{
	// Declarations:
	int numP = 0,							// total number of nodes in the mesh
		numE = 0,							// total number of elements in the mesh
		numBCe = 0,							// number of essential boundary nodes in the mesh
		numBCn = 0,							// **** number of natural boundary condition ****
		numNodesPerElem = 0,				// number of nodes per element
		dimensions = 0,						// number of dimensions in the problem
		elementType = 0,					// the elements to be used (linear, quadratic, cubic hermite etc)
		numX,								// the number of nodes in the X direction (2D only)
		numY;								// the number of nodes in the Y direction (2D only)

	int *BCenode = NULL;					// array (numP) long with 1 meaning there is an essential boundary condition
											// specified at this node, 0 means internal (solved for) node.
	int *BCnnode = NULL;				//**** array (numP) long with 1 meaning there is an natural boundary condition****
	
	int **nodesInElem = NULL;				// ordered list of nodes in each element

	double *x = NULL,						// x coordinates in global space
		   *y = NULL,						// y coordinates in global space
		   *k = NULL,						// conductivity for each element (assumes isotropic)
		   *BCevalue = NULL,				// value of essential boundary conditions at BCnodes
		   *BCnvalue = NULL,				//**** value of natural boundary conditions at BCnodes ****
		   *f = NULL,						// the RHS load vector
		   *u = NULL;						// the result vector

			
	double **K = NULL,						// Global element stiffness matrix
			**origK = NULL,					// copy of K for output
			***E = NULL;					// Array of local element stiffness matrices

	ofstream ofid;
	ifstream ifid;

	char aconfig_path[] = "FEM_Config.txt";
	CommandLineFilter *cmd = new CommandLineFilter;
		cmd->ReadFromFile(aconfig_path);

	// Reading in the problem parameters:
	dimensions = cmd->GetInt('d');
	elementType = cmd->GetInt('s');
	
	// -----------------------------------------------------------------------------------------------------------
	// 
	//	Part 1: Reading in the problem data from the files specified as follows:
	//			-n nodefilename
	//			-e elementfilename
	//			-b bcfilename
	//
	// -----------------------------------------------------------------------------------------------------------

	// Reading in the node positions from the file specified by the -n flag in FEM_Config.txt:
	ifid.open(cmd->GetString('n'));
	ifid >> numP;

	// Allocating space to store these points:
	x = new double [numP];
	if(dimensions >= 2) y = new double [numP];
	for(int i = 0; i < numP; i++) {
		ifid >> x[i];
		if(dimensions >= 2) 
			ifid >> y[i];
	}
	ifid.close();

	// Reading in the BC information from the file specified by the -b flag in FEM_Config.txt:
	ifid.open(cmd->GetString('b'));
	ifid >> numBCe;
	// Allocating space for the essential boundary conditions:
	BCenode = new int[numBCe];
	BCevalue = new double[numBCe];
	for (int i = 0; i < numBCe; i++) {
		ifid >> BCenode[i] >> BCevalue[i];
	}
	ifid.close();

	// Reading in the natural BC information:
	ifid.open(cmd->GetString('v'));
	ifid >> numBCn;
	// Allocating space for the essential boundary conditions:
	BCnnode = new int[numBCn];
	BCnvalue = new double[numBCn];
	for (int i = 0; i < numBCn; i++) {
		ifid >> BCnnode[i] >> BCnvalue[i];
	}
	ifid.close();

	// Allocating space for the other arrays:
	K = new double* [numP];
	origK = new double* [numP];
	for(int i = 0; i < numP; i++) {
		K[i] = new double [numP];
		origK[i] = new double [numP];
	}

	u = new double [numP];
	f = new double [numP];

	for(int i = 0; i < numP; i++) f[i] = u[i] = 0.0;

	// We need to read in the connectivity between nodes within an element: 
	ifid.open(cmd->GetString('e'));		// Opening the file specified by the -e flag in the FEM_Config.txt file
	ifid >> numE >> numNodesPerElem;

	// Allocating space for the element connectivity array:
	nodesInElem = new int* [numE];
	E = new double** [numE];
	for(int e = 0; e < numE; e++) {
		nodesInElem[e] = new int [numNodesPerElem];
		E[e] = new double* [numNodesPerElem];
		for(int i = 0; i < numNodesPerElem; i++) E[e][i] = new double [numNodesPerElem];
	}
	for(int e = 0; e < numE; e++) {
		for(int i = 0; i < numNodesPerElem; i++)
			ifid >> nodesInElem[e][i];
	}
	ifid.close();

	
	// Allocating space for the conductivity array:
	k = new double [MAXORDER+1];
	for(int i = 0; i <= MAXORDER; i++)
		k[i] = cmd->GetDouble('k',i+1);

	// -----------------------------------------------------------------------------------------------------------
	//
	//	Part 2: Setting up the local stiffness (E), global stiffness (K) and load (f) arrays according to the
	//					specifications in the config file:
	//							-f functionNumber
	//							-s elementType
	//							-d dimensions
	//
	// -----------------------------------------------------------------------------------------------------------

	// Element loop to fill the local stiffness E[e][m][n] matrix for each element:
	for(int e = 0; e < numE; e++) 
	{
		My_AssembleLocalElementMatrix(nodesInElem[e],x,y,k,elementType,dimensions,E[e]);
	}
	
	// Combining the local element stiffness matrices into the global stiffness matrix:
	My_AssembleGlobalElementMatrix(numP,numE,numNodesPerElem,nodesInElem,E,K);


	//Applying Essential and Natural Boundary Conditions. 

	AssembleRightHandSide(numE, nodesInElem, x, f); //This function is added to incorporate for the source term f(x)
													// depending on the source term given in .cpp file.
												//If no source term is given (Eval_SourceTerm returns 0 value)
												// then the program will output a statement stating f(x)=0.

	My_ApplyNaturalBoundaryConditions(numP, numBCn, BCnnode, BCnvalue, k, f);

	My_ApplyEssentialBoundaryConditions(numP, numBCe, BCenode, BCevalue, f, K);

	// Storing a copy for output:
	for(int i = 0; i < numP; i++) {
		for(int j = 0; j < numP; j++) {
			origK[i][j] = K[i][j];
		}
	}

	// Copying the f array to u (the array passed to the SolveAxB routine gets overwritten with the solved 'u' values
	for(int i = 0; i < numP; i++) {
		u[i] = f[i];
	}

	// -----------------------------------------------------------------------------------------------------------
	//
	//	Part 3: Solving the system
	//
	// -----------------------------------------------------------------------------------------------------------

	// RHS values passed in through u, and overwritten with new values 
	SolveAxB(K,u,numP);

	// Printing the data to a file for checking and analysis:
	ofid.open("FEM_LinearSystem.xls");
	ofid << cmd->GetLine() << endl;

	for(int i = 0; i < numP; i++) {
		for(int j = 0; j < numP; j++) {
			ofid << origK[i][j] << TAB;
		}
		ofid << TAB << u[i] << TAB << TAB << f[i] << endl;
	}
	ofid.close();

	// Printing the results to a file for visualisation using Excel:
	if(dimensions == 2) {
		if(cmd->IsPresent('x')) numX = cmd->GetInt('x');
		if(cmd->IsPresent('y')) numY = cmd->GetInt('y');
		
		// This writes an tab-delimited text file for visualisation using Excel.  Import into the meshing wizard spreadsheet using Ctrl + L
		ofid.open("FEM_Results.xls");
		for(int j = 0; j < numY; j++) {
			for(int i = 0; i < numX; i++) {
				ofid << u[numX*j+i]<<TAB;
			}
			ofid << endl;
		}
		ofid.close();
	}

	
	// -----------------------------------------------------------------------------------------------------------
	//
	//	Part 4: Housekeeping 
	//
	// -----------------------------------------------------------------------------------------------------------

	for(int i = 0; i < numP; i++) {
		delete [] K[i];
		delete [] origK[i];
	}
	delete [] K;
	delete [] origK;
	delete [] u;
	delete [] x;
	if (dimensions >= 2) delete [] y;
	delete [] BCenode;
	delete [] BCevalue;
	for(int e = 0; e < numE; e++) {
		for(int i = 0; i < numNodesPerElem; i++) delete [] E[e][i];
		delete [] nodesInElem[e];
		delete [] E[e];
	}
	delete [] nodesInElem;
	delete [] E;

	cout << endl;
	cout << "Finished! You can view the results by moving the file called \"FEM_Results.xls\" "<<endl;
	cout << "next to the \"FEM_Wizard.xls\" file, and hitting Ctrl+L to load them into the"<<endl;
	cout << "wizard for visualisation."<<endl<<endl;
	PAUSE
}