#ifndef ENTESS_H
#define ENTESS_H
#include <iostream>
#include <string.h>
#include <time.h>
#include  "apstring.h"
#include "Dataset.h"
#include "dTess.h"
#include "Descriptor.h"


using namespace std;

typedef apstring A_STRING;


void usage();
void getParameters(int , char **);
float ** trans2OldTess(Molecule &, Molecule &, int&);
void freePTS(float **, int);
void checkPath(A_STRING& );

A_STRING p_list = "";
A_STRING p_pdb = "";
A_STRING p_mol2 = "";
A_STRING p_mode = "entess";
A_STRING p_path = "";
float p_dcut = 8;

#endif

