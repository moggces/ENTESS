#ifndef Bond_H
#define Bond_H
#include <vector>
#include "apstring.h"
#include "Atom.h"


#define DUMMY_VAL	-999	

using namespace std;
typedef apstring A_STRING;

class Bond {
public:
	//constructors
	Bond(int);
	~Bond();

	//members
	Atom_n conn [2]; //when you do this, the conn is initialized ; therefore you need the default constructor
	A_STRING bond_type;

	private:
	int isLigand;
};

#endif

