#ifndef DTESS_H
#define DTESS_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <set>
#include "Tetran.h"
#include "Molecule.h"

using namespace std;

#define SQ(x)           (x) * (x)
#define BIGNUM          1E37
#define EPSILON         0.00001
#define TSIZE           100                  /* temporary storage size factor*/
#define RANGE           15.0                /* factor (>=1) for radius of control points */
#define AND             &&
#define OR              ||
#define EQ              ==
#define NE              !=
#define PI              3.14159
#define DUMMY_VAL	-999	


class Dtess {
public:
		//constructors & destructors
		Dtess();
		~Dtess();

		//functions
		void tessellate(float ** &, int);
		void tessellate(float ** &, int, int); // add an cutoff point for protein ligand to quick filter-out
		void filter_tess(Molecule& , Molecule&);  // filter out a) edge longer than cutoff 
		void displayDTPymol(Molecule& , Molecule&, A_STRING&, A_STRING&);
		vector<Tetran> tetrans;
		static float cut_off; 
		
protected:
	// functions
		double **DoubleMatrix(int , int );
		int ** IntMatrix(int , int );
		int * IntVect(int );
		void FreeVecti(int *);
		void FreeMatrixi(int **);
		void FreeMatrixd(double **);
			
private:
			
	
};

#endif

