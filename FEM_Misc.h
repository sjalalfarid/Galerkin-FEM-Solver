#ifndef FEM_MISC

#include <stdlib.h>
#include <iostream>

using namespace std;

// Defining a macro for use in pausing the run-time of the code:
#ifndef PAUSE
#define PAUSE {char pause; cout << "PAUSED: Hit any character and enter to continue."<<endl; cin >> pause;}
#endif 

// Defining a macro for highest order differential equation that can be solved using this code:
#define MAXORDER 2	

// Defining the TAB macro for easy output to Excel:
#define TAB "\t"

#define FEM_MISC
#endif