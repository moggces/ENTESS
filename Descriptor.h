#ifndef Descriptor_H
#define Descriptor_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include "apstring.h"
#include "Tetran.h"
#include "Molecule.h"

#define DUMMY_VAL	-999

using namespace std;
typedef apstring A_STRING;

struct Desc {
	int tetran_type[4]; 
	float value; 
};

class Descriptor {

public:
	Descriptor(int, A_STRING); 
	~Descriptor();

	//static members
	static int desc_type; // 0: 554 ...
	static int charge; //0: occurence (no charge); 1:oldEN; 2: DFT_EN 3; DFT_MCT ...
	static int proc; //0: tetra-summation; 1: pl_edge-summation 2: volume summation
	static int counter;

	// functions
	void genDesc(vector<Tetran>&);
	void combDesc(A_STRING, A_STRING, A_STRING, A_STRING, A_STRING, bool);
	void printDesc(A_STRING, A_STRING, A_STRING, A_STRING);
	void displayOLDTetranPymol(Molecule& , Molecule& , A_STRING , A_STRING, vector<Tetran>& );



protected:
	void initializeOldDesc();
	void setDescType(A_STRING);
	A_STRING printTitle();
	A_STRING getDescName(int []);
	vector<Desc> descs;
	static A_STRING desc_names;

private:
	
	
};
#endif
