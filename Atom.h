#ifndef Atom_H
#define Atom_H
#include <vector>
#include "apstring.h"

#define DUMMY_VAL	-999	

using namespace std;
typedef apstring A_STRING;

struct Protein_S {
	float temper;
};

struct Ligand_S {
	float p_charge;
	
};

class Atom_n {

public:

	//constructors , deconstructor
	Atom_n(int); //0 is ligand ; the other numbers are protein ... but can redefine to other types 
	Atom_n();
	~Atom_n();

	//overloaded operators
	Atom_n& operator = ( const Atom_n&);


	int atom_number;
	int isLigand;
	A_STRING atomType; // CA, O, or C.3 in mol2
	A_STRING atom_name;
	A_STRING res_name;
	A_STRING chain;
	int res_num;
	float x, y, z;

	Protein_S p_s;
	Ligand_S l_s;

private:
	


};
#endif

