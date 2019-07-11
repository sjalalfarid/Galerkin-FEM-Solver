/*

	The source function f(x)

*/

#include "FEM_SourceTerm.h"
#include <math.h>   // ************ The Math.h Header Must be added in order to define generic function *******

//
//

double Eval_SourceTerm(double point)
{
	return(0);			// For f(x)=0, That is for all homogeneous differential equations.


	//****** Simple Test Functions *******

	//return (1.0);     // for f(x)=1, When Function on RHS is a constant, here it is 1. If zero then its the same as homogenous diff. eq.
	//return (point);   // for f(x)=x. Simple Test Function
	
	//****** For Any other generic function just type the function directly (Make sure math.h header is added) *******
	//****** Below is The f(x)=10e^5x-4e^-x. Assignment Task 1 Question, Section 2.3. 

	//return ((10 * (pow(2.71828, (-5 * point)))) - ((4 * (pow(2.71828, (-point))))));  //For Advection Diffusion Problem Task 2.3 of assignment f(x)=10e^5x-4e^-x.
																		
}
