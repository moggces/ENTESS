#include "Bond.h"
using namespace std;

/*
is there a better way for initialization?
*/

Bond::Bond(int isL)
{
	isLigand = isL;
	Atom_n atom1(isLigand);
	conn [0] = atom1; // is there a better way to do that?
	conn [1] = atom1;
	bond_type = "";
}

Bond::~Bond()
{
}

