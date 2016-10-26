#ifndef DATASET_H
#define DATASET_H
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "apstring.h"
#include "Atom.h"
#include "Bond.h"
#include "dTess.h"
#include "Descriptor.h"
#include "Molecule.h"

using namespace std;
typedef apstring A_STRING;

/*
haven't remove the H bond
haven't fixed the problem of occupency
*/

class Dataset {

public:
	//constructors & destructors
	Dataset();
	Dataset(const Dataset&);
	~Dataset();
	//functions
	void importPDB(A_STRING&, A_STRING ); //path, pdb_name
	void importMOL2(A_STRING&, A_STRING ); //path, mol2_name
	void importMOL2_tessellate(A_STRING&, A_STRING, A_STRING, A_STRING, double, Molecule& );  //path, mol2_name, mode, protein
	static A_STRING basename;
	vector<Molecule> moles; //maybe Molecule mole is better ...?

protected:
	void checkExt(A_STRING, A_STRING);
	void checkPath(A_STRING& );
	void removeH(vector<Atom_n>&);
	

private:
	void tessellation_wrap(Molecule &, Molecule &, double, Dtess&);
	float ** trans2OldTess(Molecule &, Molecule &, int&);
	void freePTS(float **, int);
	
};
#endif
