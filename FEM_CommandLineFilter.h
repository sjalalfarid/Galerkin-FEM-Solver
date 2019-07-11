#ifndef COMMANDLINEFILTER

#define COMMANDLINEFILTER

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <stdio.h>
#include "FEM_Misc.h"

using namespace std;

class CommandLineFilter {
	
public:
	CommandLineFilter(){Initialise(); SetDefaults();}
	~CommandLineFilter(){Destroy();}

	void Initialise();
	void SetDefaults();
	void Destroy();

	void SetInput(int argc, char **argv);
	void SetUsage(string infile);
	void SetUsage(char *infile);
	void Usage();
	void ReadFromFile(char *filename);

	int IsPresent(char get);
	int GetInt(char get);
	int GetInt(char get, int n);
	double GetDouble(char get);
	double GetDouble(char get, int n);
	char* GetString(char get);
	char* GetString(char get, int n);
	char* GetLine();

protected:
	int num,
			readFromFile;
	char **line;
	string usagefile;

};

#endif