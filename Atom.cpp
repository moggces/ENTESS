#include "Atom.h"
using namespace std;




Atom_n::Atom_n(int isL)
{
	isLigand= isL;
	
	atom_number = DUMMY_VAL;
	atomType = "";
	atom_name = "";
	res_name = "";
	chain = "";
	res_num =DUMMY_VAL;
	//EN =DUMMY_VAL;
	 x =y = z=DUMMY_VAL;
	
	 if (isLigand == 0)
	 {
		 l_s.p_charge = DUMMY_VAL;
	 } else
	 {
		 p_s.temper = DUMMY_VAL;
	 }

}
Atom_n::Atom_n()
{
	atom_number = DUMMY_VAL;
	atomType = "";
	atom_name = "";
	res_name = "";
	chain = "";
	res_num =DUMMY_VAL;
	//EN =DUMMY_VAL;
	 x =y = z=DUMMY_VAL;
}

Atom_n::~Atom_n()
{
	
}

Atom_n& Atom_n::operator = ( const Atom_n& atom1)
{
	if(this == &atom1 ) { return *this; }// check for self assignment

		this->atom_number = atom1.atom_number;
		this->atomType = atom1.atomType;
		this->atom_name = atom1.atom_name;
		this->res_name = atom1.res_name;
		this->chain = atom1.chain;
		this->res_num = atom1.res_num;
		//this->EN = atom1.EN; 
		this->x = atom1.x; this->y = atom1.y; this->z = atom1.z;
		this->isLigand = atom1.isLigand;
		if (atom1.isLigand == 0)
		{
			this->l_s = atom1.l_s;
		} else { this->p_s = atom1.p_s;}

		return *this;
}
