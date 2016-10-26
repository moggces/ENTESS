#include "Tetran.h"
using namespace std;

Tetran::Tetran()
{
	Atom_n atom1(0); //assume all of them are ligand in order to initialize
	for (int i=0; i<4; ++i)
	{
		v[i] = DUMMY_VAL;
		nodes[i] = atom1;
		EN[i] = DUMMY_VAL;
		tetran_type[i] = DUMMY_VAL;
	}
	for (int i = 0; i< 6; ++i)
	{
		edges.dist[i] = DUMMY_VAL;
	}
}

Tetran::~Tetran()
{

}

void Tetran::sort_quad(Tetran &qt)
{
  int vtemp,i,j;
 
  for(i=0;i<3;i++)
     for(j=0;j<3-i;j++)
	if(qt.v[j+1] < qt.v[j])
	   {
	     vtemp=qt.v[j];
	     qt.v[j]=qt.v[j+1];
	     qt.v[j+1]=vtemp;
	   }
}


int Tetran::get_Tetran_edges_distsq(Tetran& qt, float cut_off)
{
	int c=0;
	int j=0;
	//int result = 0;
	for (int i = 0; i< 3; ++i)
	{
		++j;
		for (int k=j; k<4; ++k)
		{
			float dist = qt.get_edge_distsq(qt.nodes[i], qt.nodes[k]);
			if (  dist > cut_off) { 
				return 1;
			} else { qt.edges.dist[c++] = dist; }
		}
	}
	return 0; // all edges are generated
}

float Tetran::get_edge_distsq(Atom_n& atom1, Atom_n& atom2)
{
	return ( pow(atom1.x-atom2.x, 2) + pow(atom1.y-atom2.y, 2) + pow(atom1.z-atom2.z, 2) );
}

bool Tetran::isPLtetran(Tetran& qt, int n_p_atoms)
{
	if ( (qt.v[0]-1 < n_p_atoms && qt.v[1]-1 < n_p_atoms && qt.v[2]-1 < n_p_atoms && qt.v[3]-1 < n_p_atoms ) ||
			(qt.v[0]-1 >= n_p_atoms && qt.v[1]-1 >= n_p_atoms && qt.v[2]-1 >= n_p_atoms && qt.v[3]-1 >= n_p_atoms ) )
	{ return false;} else {return true;}
}

void Tetran::assignOLDtetranType()
{
	for (int i = 0; i < 4; ++i)
	{
		this->tetran_type[i] = getOLDAtomType(this->nodes[i]);
	}
	sort_TetranType(*this);
	
}

void Tetran::sort_TetranType(Tetran& qt)
{
	//cout << "here to sort vertex type"<< endl;    //vertex_type from small to large
	int tmp, i, j;
	for (i=0; i < 4; i++)
	{
		for (j=i+1; j<4; j++)
		{
			if (qt.tetran_type[i] >= qt.tetran_type[j])
			{
				tmp = qt.tetran_type[i];
				qt.tetran_type[i] = qt.tetran_type[j];
				qt.tetran_type[j] = tmp;
			}
		}
	}
}


int Tetran::getOLDAtomType(Atom_n& atom)
{
	if (atom.atomType == "") { cout << "Ligand?" << atom.isLigand << " atomType hasn't defined" << atom.atom_number << endl; exit(1);} 
	if (atom.isLigand == 0)
	{
		if (atom.atomType.find("Na") == 0 || atom.atomType.find("K") == 0 || atom.atomType.find("Zn") == 0 || atom.atomType.find("Cu") == 0 ||
			atom.atomType.find("Ca") == 0 || atom.atomType.find("Mn") == 0 || atom.atomType.find("Fe") == 0 || atom.atomType.find("Co") == 0 ||
			atom.atomType.find("Mg") == 0 || atom.atomType.find("Ni") == 0 || atom.atomType.find("Au") == 0 || atom.atomType.find("Al") == 0 ||
			atom.atomType.find("Sn") == 0 || atom.atomType.find("Li") == 0 || atom.atomType.find("Hg") == 0 || atom.atomType.find("B") == 0 || atom.atomType.find("Cd") == 0)
		{ return 1; }
		else if (atom.atomType.find("H") == 0 || atom.atomType.find("Du") == 0 || atom.atomType.find("LP") == 0 ) {return 0;}
		else if (atom.atomType.find("Cl") == 0 || atom.atomType.find("Br") == 0 || atom.atomType.find("I") == 0 || atom.atomType.find("F") == 0 ||
			atom.atomType.find("P.") == 0 || atom.atomType.find("P") == 0 ) {return 2;}
		else if (atom.atomType.find("S.") == 0 ) {return 3;}
		else if (atom.atomType.find("C.") == 0 ) {return 4;}
		else if (atom.atomType.find("N.") == 0 ) {return 5;}
		else if (atom.atomType.find("O.") == 0 ) {return 6;}
		else 
		{
			cout << "unknown atomType: "<< atom.atomType << " Ligand? " << atom.isLigand << " atom-number: " << atom.atom_number << endl; 
			return 0; //H
		}
	} else 
	{
		if ( *atom.atomType.c_str()  == 'H'  || isdigit(*atom.atomType.c_str() )) { return 99; }
		else if (*atom.atomType.c_str() == 'S'  || *atom.atomType.c_str() == 'P' ) {return 11;}
		else if (*atom.atomType.c_str()  == 'C') {return 12;}
		else if (*atom.atomType.c_str()  == 'N') {return 13;}
		else if (*atom.atomType.c_str()  == 'O') {return 14;}
		else 
		{
			cout << "unknown atomType: "<< atom.atomType << " Ligand? " << atom.isLigand << " atom-number: " << atom.atom_number << endl; 
			return 99; //H
		}
	}

	
}

void Tetran::assignPaulingEN()
{
	for (int i = 0; i < 4; ++i)
	{
		this->EN[i] = getAtomPaulingEN(this->nodes[i]);
	}
}

float Tetran::getAtomPaulingEN(Atom_n& atom)
{
	if (atom.atomType == "") { cout << "Ligand? " << atom.isLigand << " atomType hasn't defined " << atom.atom_number << endl; exit(1);}
	if (atom.atomType.find("Na") == 0 || atom.atomType.find("Ca") == 0 || atom.atomType.find("Li") == 0 || atom.atomType.find("Ca") == 0 ) 
	{return 1.0;}
	else if (atom.atomType.find("K") == 0 || atom.atomType.find("Cs") == 0 ) {return 0.9;}
	else if (atom.atomType.find("Mg") == 0 ) {return 1.2;}
	else if (atom.atomType.find("Au") == 0 || atom.atomType.find("Ag") == 0 || atom.atomType.find("Pt") == 0 || atom.atomType.find("Pd") == 0 )
	{ return 1.4;}
	else if (atom.atomType.find("Hg") == 0 || atom.atomType.find("Al") == 0 ) {return 1.5;}
	else if (atom.atomType.find("Mn") == 0 || atom.atomType.find("Fe") == 0 || atom.atomType.find("Cr") == 0 || atom.atomType.find("Pb") == 0 || atom.atomType.find("Cd") == 0 ) 
	{return 1.6;}
	else if (atom.atomType.find("Zn") == 0 || atom.atomType.find("Sn") == 0 || atom.atomType.find("Si") == 0 || atom.atomType.find("Co") == 0 || atom.atomType.find("Bi") == 0 ) 
	{return 1.7;}
	else if (atom.atomType.find("Ni") == 0 || atom.atomType.find("Cu") == 0 ) {return 1.8;}
	else if (atom.atomType.find("B") == 0 ) {return 2.0;}
	else if (atom.atomType.find("P") == 0 || atom.atomType.find("P.") == 0 ) {return 2.1;}
	else if (atom.atomType.find("H") == 0 || atom.atomType.find("Du") == 0 || atom.atomType.find("LP") == 0 || isdigit(*atom.atomType.c_str() ) ||  atom.atomType.find("I") == 0 || atom.atomType.find("As") == 0 )
	{ return 2.2;}
	else if (atom.atomType.find("Cl") == 0 || atom.atomType.find("Br") == 0 ) {return 2.8;}
	else if (atom.atomType.find("C") == 0 || atom.atomType.find("Se") == 0 || atom.atomType.find("C.") == 0) {return 2.5;}
	else if (atom.atomType.find("S") == 0 || atom.atomType.find("S.") == 0 ) {return 2.4;}
	else if (atom.atomType.find("N") == 0 || atom.atomType.find("N.") == 0 || atom.atomType.find("AD") == 0 || atom.atomType.find("AE") == 0 )
	{ return 3.1;}
	else if (atom.atomType.find("O") == 0 || atom.atomType.find("O.") == 0 ) {return 3.5;}
	else if (atom.atomType.find("F") == 0 ) {return 4.1;}
	else { 
		cout << "Ligand? " << atom.isLigand << " don't match any of known atomtypes " << atom.atomType << " number: " << atom.atom_number << " assigned as H" << endl;
		return 2.2; //H
	}
}

void Tetran::assignDFT_MCT()
{
	for (int i = 0; i < 4; ++i)
	{
		this->EN[i] = getAtomDFT_MCT(this->nodes[i]);
	}
}

float Tetran::getAtomDFT_MCT(Atom_n& atom)
{
	if (atom.atomType == "") { cout << "Ligand? " << atom.isLigand << " atomType hasn't defined " << atom.atom_number << endl; exit(1);}
	if (atom.atomType.find("Ag") == 0 ) {return 0.71;}
	else if (atom.atomType.find("Al") == 0) { return 0.58;}
	else if (atom.atomType.find("As") == 0) {return 0.59;}
	else if (atom.atomType.find("Au") == 0) {return 0.83;}
	else if (atom.atomType.find("Ba") == 0) {return 0.53;}
	else if (atom.atomType.find("Br") == 0) {return 0.9;}
	else if (atom.atomType.find("B") == 0) {return 0.54;}
	else if (atom.atomType.find("Ca") == 0) {return 0.5;}
	else if (atom.atomType.find("Mg") == 0) 
	{cout << "No Mg; use Ca instead" << endl; return 0.5;} // add Mg ... which paper doesn't provide the value
	else if (atom.atomType.find("Co") == 0) {return 0.59;}
	else if (atom.atomType.find("Cr") == 0) {return 0.61;}
	else if (atom.atomType.find("Cs") == 0) {return 0.64;}
	else if (atom.atomType.find("Cu") == 0) {return 0.69;}
	else if (atom.atomType.find("Zn") == 0 ) 
	{cout << "No Zn; use Cu instead" << endl; return 0.69;} // add Zn ... which paper doesn't provide the value
	else if (atom.atomType.find("Fe") == 0) {return 0.52;}
	else if (atom.atomType.find("K") == 0 || atom.atomType.find("Li") == 0 ) {return 0.63;}
	else if (atom.atomType.find("Mn") == 0 ) {return 0.5;}
	else if (atom.atomType.find("Na") == 0 ) {return 0.62;}
	else if (atom.atomType.find("Ni") == 0 ) {return 0.68;}
	else if (atom.atomType.find("Pb") == 0 ) {return 0.55;}
	else if (atom.atomType.find("Pd") == 0 ) {return 0.57;}
	else if (atom.atomType.find("Pt") == 0 ) {return 0.81;}
	else if (atom.atomType.find("Si") == 0 ) {return 0.71;}
	else if (atom.atomType.find("Sn") == 0 ) {return 0.68;}

	else if (atom.atomType.find("Cl") == 0) {return 0.89;} // bug in old program: XlXlClCr ; if X == Cl ; then the related values are not correct
	else if (atom.atomType.find("F") == 0) {return 0.74;}
	else if (atom.atomType.find("I") == 0 ) {return 0.91;}
	else if (atom.atomType.find("P") == 0 ||  atom.atomType.find("P.") == 0 ) {return 0.58;}
	else if (atom.atomType.find("N") == 0 || atom.atomType.find("N.") == 0 || atom.atomType.find("AD") == 0 || atom.atomType.find("AE") == 0 ) 
	{return 0.5;}
	else if (atom.atomType.find("O") == 0 || atom.atomType.find("O.") == 0 ) {return 0.62;}
	else if (atom.atomType.find("S") == 0 || atom.atomType.find("S.") == 0 ) {return 0.75;}
	else if (atom.atomType.find("C") == 0 || atom.atomType.find("C.") == 0 ) {return 0.63;}
	else if (atom.atomType.find("H") == 0 || atom.atomType.find("Du") == 0 || atom.atomType.find("LP") == 0 || isdigit (*atom.atomType.c_str()) ) 
	{return 0.56;}
	else {
		cout << "Ligand? " << atom.isLigand << " don't match any of known atomtypes " << atom.atomType << " number: " << atom.atom_number << " assigned as H" << endl;
		return 0.56;
	}
}

float Tetran::get_Tetran_volume()
{
	float a[5][5], value;
	int i, z = 0;
	
	for (i=0; i<5; ++i) 
	{
		for (z=0; z<5; ++z)
		{
			if (i==z)
				a[i][i]=0;
			else if ( (i==0 || z==0) && i!=z ) 
				a[i][z]=1;
			else
				a[i][z] =  get_edge_distsq(this->nodes[i-1], this->nodes[z-1]);
		}
	}
	
	if(chckdgnl(a,5)==0)
		value=0;
	else
		value=deter(a,5);

	return sqrt((float) value/288);
}

float Tetran::deter(float a[][5],int forder)
{
  int i,j,k;
  float mult;
  float deter=1;
  for(i=0;i<forder;i++)
  {
        for(j=0;j<forder;j++)
        {
          if (a[i][i] == 0)    //Jui consider 3X3 matrix 1 2 3 4 5 6 7 8 9
                mult=0;
          else
                mult=a[j][i]/a[i][i];
          for(k=0;k<forder;k++)
          {
                if(i==j) break;
                a[j][k]=a[j][k]-a[i][k]*mult;
          }
        }
  }
  for(i=0;i<forder;i++)
  {
        //printf("a[i][i]: %f \n", a[i][i]);
        deter=deter*a[i][i];
  }
  return(deter);
}


int Tetran::chckdgnl(float array[][5],int forder)
{
  int i,j,k=0;
  for(i=0;i<forder;i++)
  {
         if(array[i][i]==0)
         {
                for(j=0;j<forder;j++)
                {
                  if(array[i][j]!=0)
                  {
                         k=j;
                         break;
                  }
                  if(j==(forder)) //forder-1
                         return(0);
                }
                for(j=0;j<forder;j++)
                {
                  array[j][i]=array[j][i]-array[j][k];
                }
         }
  }
  return(1);
}

