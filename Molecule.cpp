#include "Molecule.h"
using namespace std;

Molecule::Molecule(int isL)
{
	isLigand = isL;
	// * atoms = new vector<Atom_n>;
	 //* bonds = new vector<Bond>;
	
}

Molecule::Molecule()
{
}


Molecule::~Molecule()
{
	atoms.clear();
	bonds.clear();
}

