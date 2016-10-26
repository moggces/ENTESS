#ifndef Tetran_H
#define Tetran_H
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Atom.h"

/*
volume funtion hasn't fully checked
*/

#define DUMMY_VAL	-999	

using namespace std;
typedef apstring A_STRING;

struct Edges {
	float dist[6]; // 01, 02, 03, 12, 13, 23
};

class Tetran {

public:
	Tetran();
	~Tetran();

	// essential of each tetrahedron
	Atom_n nodes[4]; // the order is based on v[4]
	Edges edges;
	int v[4]; // for internal use; sorted based on atom_number (and ligand is behind protein) ; start from 1...
	//i would like to set it as private but it won't be accessed in the Dtess class

	//Descriptor purpose
	double EN[4]; // could be either MCT or oldEN; order is same as Atom_n[4]'s order

	int tetran_type[4]; // sorted based on atom_type (int) ; independent fingerprint of this tetran;  didn't match the order in nodes[4], which is based on v[4]

	

	//functions
	
	int get_Tetran_edges_distsq(Tetran&, float);
	float get_Tetran_volume();
	
	bool isPLtetran(Tetran&, int);
	void assignOLDtetranType(); //try this->
	void assignPaulingEN();
	void assignDFT_MCT();

	void sort_quad(Tetran &); // sort v[4]
	void sort_TetranType(Tetran &); // sort tetran_type[4]

protected:
	
	
	int getOLDAtomType(Atom_n&);
	float getAtomPaulingEN(Atom_n&);
	float getAtomDFT_MCT(Atom_n&);
	float get_edge_distsq(Atom_n&, Atom_n&);
	float deter(float [][5],int );
	int chckdgnl(float [][5],int );

private:
	



};

#endif
