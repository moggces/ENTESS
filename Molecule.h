#ifndef Molecule_H
#define Molecule_H
#include <vector>
#include "apstring.h"
#include "Atom.h"
#include "Bond.h"
/*
when to delete the pointers??
why I can't use * atoms?? in EnTess it will become 0
*/

#define DUMMY_VAL	-999	

using namespace std;
typedef apstring A_STRING;

class Molecule {

public:

	//constructors , deconstructor
	Molecule(int); 
	Molecule();
	~Molecule();

	vector<Atom_n>  atoms; 
	vector<Bond>  bonds;
	

private:
	int isLigand;
	

};

#endif

