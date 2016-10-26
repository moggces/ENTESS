#include "Dataset.h"
#include "Atom.h"
#include "Bond.h"
using namespace std;


Dataset::Dataset()
{
	//need a vector to store atoms
	

}

Dataset::~Dataset()
{
	moles.clear();
	
}

void Dataset::checkExt (A_STRING file, A_STRING suffix)
{
	int pos_suffix = file.find(suffix); //starting from 0
	if (pos_suffix != 0)
	{
		basename = file.substr(0, pos_suffix-1);
	} else {
		cout << "The file extension indicates it's not a " << suffix << endl;
		exit (1);
	}

}

void Dataset::checkPath(A_STRING& path)
{
	int pos;
	int ll =path.length();
	if ( (pos = path.find('/') ) >= 0 )
	{
		// it is a linux path)
		if (path[ll-1] != '/')
		{
			path += '/';
		}
	} else // it is a windows path 
	{
		if (path[ll-1] != '\\')
		{
			path += '\\';
		}
	}
}

A_STRING Dataset::basename = "";

void Dataset::removeH(vector<Atom_n>& atoms) 
{
	  //cout << atoms.size() << endl;
	  vector<Atom_n>::iterator it;
	  for (it=atoms.begin(); it<atoms.end(); )
	  {
		  if ((*it).atomType.find("H") == 0 || isdigit(*((*it).atomType.c_str())))
		  {
			  it = atoms.erase(it);
		  } else { ++it;}
	  }
	  //cout << atoms.size() << endl;
}

void Dataset::importPDB (A_STRING& path, A_STRING p_pdb)
{
	
	checkPath(path);
	checkExt(p_pdb, "pdb");
	
	A_STRING line;
	vector<Atom_n> atoms;

	A_STRING file_name = path+p_pdb;
	ifstream mypdbfile (file_name.c_str());
	if (mypdbfile.is_open())
	{
		while (! mypdbfile.eof() )
		{
			line.getline(mypdbfile);
			if ( line.find("ATOM") == 0 ) // only the ATOM at the begining Metal should be in ligand or protein?
			{
				if ( line.find("HOH") < 0) // no water
				{
					Atom_n pp(1); 
					pp.atomType = line.substr(12, 4); pp.atomType.parse_string();
					pp.atom_name = pp.atomType;
					if (pp.atomType.find("H") == 0 || isdigit(*(pp.atomType.c_str())) )
					{ continue;}
					else 
					{
						pp.atom_number = atoi(line.substr(6,5).c_str());
						pp.res_name = line.substr(17, 3); pp.res_name.parse_string();
						pp.chain = (line.substr(21, 1)); pp.chain.parse_string(); 
						pp.res_num = atoi(line.substr(22,4).c_str());
						pp.x = atof(line.substr(30,8).c_str()); pp.y = atof(line.substr(38,8).c_str());  pp.z = atof(line.substr(46,8).c_str());
						pp.p_s.temper = atof(line.substr(60, 6).c_str());
						//pp.assignOLDAtomType();
						atoms.push_back(pp);
					}
				}
			}
			
		}
		//removeH(atoms);
		Molecule mole(1);
		mole.atoms = atoms;
		moles.push_back(mole);
		atoms.clear();
	} else { cout << "cannot open pdb file" <<endl; exit (1);}
	
}

void Dataset::importMOL2 (A_STRING& path, A_STRING p_mol2)
{
	checkPath(path);
	checkExt(p_mol2, "mol2");
	basename = p_mol2.substr(0, p_mol2.find('.'));
	vector<Atom_n> atoms;
	vector<Bond> bonds;

	A_STRING line;
	A_STRING file_name = path+p_mol2;
	ifstream mymol2file (file_name.c_str());
	if (mymol2file.is_open())
	{
		while (! mymol2file.eof() )
		{
			line.getline(mymol2file);
			if (line.find("@<TRIPOS>ATOM") >= 0)
			{
				line.getline(mymol2file);
				do {
					Atom_n ll(0);
					line += " ";
					A_STRING column;
					stringstream ss(line.c_str());
					for ( unsigned int i=1;  ss >> column; ++i)
					{
						switch (i)
						{
							case 1:
								ll.atom_number = atoi(column.c_str()); break;
							case 2:
								ll.atom_name = column; break;
							case 3:
								ll.x = atof(column.c_str()); break;
							case 4:
								ll.y = atof(column.c_str()); break;
							case 5:
								ll.z = atof(column.c_str()); break;
							case 6:
								ll.atomType = column; break;
							case 9:
								ll.l_s.p_charge = atof(column.c_str()); break;
						}
					} 
					//if (i == 9) {ll.l_s.p_charge = atof(column.c_str());}
					ll.res_name = "LIG";
					ll.res_num = 0;
					ll.chain = '0';
					//ll.assignOLDAtomType();
					atoms.push_back(ll); // why ll.atomType != NULL doesn't work??
					line.getline(mymol2file);
				} while ( line.find("@<TRIPOS>BOND") < 0);
			} 
			if (line.find("@<TRIPOS>BOND") >= 0)
			{
				line.getline(mymol2file);				
				do {
					Bond lb(0);
					line += " ";
					A_STRING column;
					stringstream ss(line.c_str());
					unsigned int i = 1;
				
					while ( ss >> column )
					{
						switch (i)
						{
							case 2:	
								lb.conn[0] = atoms[atoi(column.c_str())-1];  break;
							case 3: 
								lb.conn[1]  = atoms[atoi(column.c_str())-1]; break;
							case 4: //column4
								lb.bond_type = column; break;  // can't jump to this if there's no space after the last column... why?
						}
						++i;
					}
					//if (i == 4) {lb.bond_type = column;}
					bonds.push_back(lb); 
					line.getline(mymol2file);
				} while ((! mymol2file.eof() ) && line.find("@<TRIPOS>MOLECULE") < 0 && line.find("@<TRIPOS>SUBSTRUCTURE") < 0);
			}
			
			if (line.find("@<TRIPOS>MOLECULE") >= 0 && atoms.size() != 0)
			{
				
				Molecule mole(0);
				removeH(atoms);
				if (atoms.size() == 0) { cout << "There is no ligand atoms. Please check the mol2 format" << endl; exit(1);}
				mole.atoms = atoms;
				mole.bonds = bonds;
				moles.push_back(mole);
				atoms.clear();
				bonds.clear();
			}

		}
		
		Molecule mole(0);
		removeH(atoms);
		if (atoms.size() == 0) { cout << "There is no ligand atoms. Please check the mol2 format" << endl; exit(1);}
		mole.atoms = atoms;
		mole.bonds = bonds;
		moles.push_back(mole);
		atoms.clear();
		bonds.clear();
	} else {  cout << "can't open the mol2 file ... " << endl; exit (1);}

}

void Dataset::importMOL2_tessellate (A_STRING& path, A_STRING p_mol2, A_STRING p_mode, A_STRING p_list, double p_dcut, Molecule& protein)
{
	checkPath(path);
	checkExt(p_mol2, "mol2");
	basename = p_mol2.substr(0, p_mol2.find('.'));
	vector<Atom_n> atoms;
	vector<Bond> bonds;

	A_STRING line;
	A_STRING file_name = path+p_mol2;
	A_STRING pose_name;
	bool isSingleMol;
	int counter=0;

	ifstream mymol2file (file_name.c_str());
	if (mymol2file.is_open())
	{
		while (! mymol2file.eof() )
		{
			line.getline(mymol2file);
			if (line.find("@<TRIPOS>ATOM") >= 0)
			{
				line.getline(mymol2file);
				do {
					Atom_n ll(0);
					line += " ";
					A_STRING column;
					stringstream ss(line.c_str());
					for ( unsigned int i=1;  ss >> column; ++i)
					{
						switch (i)
						{
							case 1:
								ll.atom_number = atoi(column.c_str()); break;
							case 2:
								ll.atom_name = column; break;
							case 3:
								ll.x = atof(column.c_str()); break;
							case 4:
								ll.y = atof(column.c_str()); break;
							case 5:
								ll.z = atof(column.c_str()); break;
							case 6:
								ll.atomType = column; break;
							case 9:
								ll.l_s.p_charge = atof(column.c_str()); break;
						}
					} 
					//if (i == 9) {ll.l_s.p_charge = atof(column.c_str());}
					ll.res_name = "LIG";
					ll.res_num = 0;
					ll.chain = '0';
					//ll.assignOLDAtomType();
					atoms.push_back(ll); // why ll.atomType != NULL doesn't work??
					line.getline(mymol2file);
				} while ( line.find("@<TRIPOS>BOND") < 0);
			} //if (line.find("@<TRIPOS>ATOM") == 0)

			if (line.find("@<TRIPOS>BOND") >= 0)
			{
				line.getline(mymol2file);				
				do {
					Bond lb(0);
					line += " ";
					A_STRING column;
					stringstream ss(line.c_str());
					unsigned int i = 1;
				
					while ( ss >> column )
					{
						switch (i)
						{
							case 2:	
								lb.conn[0] = atoms[atoi(column.c_str())-1];  break;
							case 3: 
								lb.conn[1]  = atoms[atoi(column.c_str())-1]; break;
							case 4: //column4
								lb.bond_type = column; break;  // can't jump to this if there's no space after the last column... why?
						}
						++i;
					}
					//if (i == 4) {lb.bond_type = column;}
					bonds.push_back(lb); 
					line.getline(mymol2file);
				} while ( (! mymol2file.eof() ) && line.find("@<TRIPOS>MOLECULE") < 0 && line.find("@<TRIPOS>SUBSTRUCTURE") < 0);
			}//if (line.find("@<TRIPOS>BOND") == 0)

			if (line.find("@<TRIPOS>MOLECULE") >= 0 )
			{
				if (atoms.size() != 0)
				{
					Molecule mole(0);
					counter++;
					removeH(atoms);
					mole.atoms = atoms;
					mole.bonds = bonds;
					//moles.push_back(mole);
					atoms.clear();
					bonds.clear();

					//tessellation
					Dtess dt;
					tessellation_wrap(protein, mole, p_dcut, dt);

					// descriptor generation
					Descriptor dst(0, p_mode); 
					cout << "Pose ID: " << dst.counter << endl;
					//dst.setDescType(p_mode);
					dst.genDesc(dt.tetrans);
					isSingleMol = false;
					dst.combDesc(path, basename, pose_name, p_mode, p_list, isSingleMol);

					pose_name.getline(mymol2file);
				}  //if (atoms.size() != 0)
				else {
					pose_name.getline(mymol2file);
				}
			}//if (line.find("@<TRIPOS>MOLECULE") == 0 && atoms.size() != 0)

		}//main loop, eof()...

		
		Molecule mole(0);
		counter++;
		removeH(atoms);
		if (atoms.size() == 0) { cout << "There is no ligand atoms. Please check the mol2 format" << endl; exit(1);}
		mole.atoms = atoms;
		mole.bonds = bonds;
		//moles.push_back(mole);
		atoms.clear();
		bonds.clear();

		//tessellation
				Dtess dt;
				tessellation_wrap(protein, mole, p_dcut, dt);
				
				// descriptor generation
				Descriptor dst(0, p_mode); 
				cout << "Pose ID: " << dst.counter << endl;
				//dst.setDescType(p_mode);
				dst.genDesc(dt.tetrans);

				if ( counter > 1 )
				{
					isSingleMol = false;
				} else { 
					isSingleMol = true;
					dt.displayDTPymol(protein,mole, path, basename);
					dst.displayOLDTetranPymol (protein,mole, path, basename, dt.tetrans);
				}
				dst.combDesc(path, basename, pose_name, p_mode, p_list, isSingleMol);

				dst.printDesc(path, basename, p_mode, p_list);
	} else {  cout << "can't open the mol2 file ... " << endl; exit (1);}

}

float ** Dataset::trans2OldTess(Molecule& protein, Molecule& ligand, int& n_atoms)
{
	//pointer initialization
	float ** pts;
	try { pts = new (float *[n_atoms+4]); }
	catch (bad_alloc&) { cout << "Error allocating memory." << endl;}

	for (int i = 0; i<n_atoms+4; ++i) //+4 is necessary , could be the bugs in original tessellation program
	{
		try {pts[i] = new float[3];}
		catch (bad_alloc&) { cout << "Error allocating memory." << endl;}
		pts[i][0] = 0.0;
		pts[i][1] = 0.0; 
		pts[i][2] = 0.0;
	}

	vector<Atom_n>::iterator it;
	int j=0;
	for (it=protein.atoms.begin(); it<protein.atoms.end(); ++it)
	{
		pts[j][0] = (*it).x; pts[j][1] = (*it).y; pts[j][2] = (*it).z;
		++j;
	}
	for (it=ligand.atoms.begin(); it<ligand.atoms.end(); ++it)
	{
		pts[j][0] = (*it).x; pts[j][1] = (*it).y; pts[j][2] = (*it).z;
		++j;
	}

	return pts;
	
}

void Dataset::freePTS(float ** pts, int n_atoms)
{
	for (int k=0; k<n_atoms+4; k++)
	{
			if (pts[k])
			{
				delete[] pts[k];
				pts[k]= NULL;
			}
	}
	delete[] pts;
	pts = NULL;
}

void Dataset::tessellation_wrap(Molecule& protein, Molecule& mole, double p_dcut, Dtess& dt)
{
	float ** pts;
	int n_atoms = (protein.atoms.size()+ mole.atoms.size());
	pts = trans2OldTess(protein, mole, n_atoms); 
	dt.tessellate(pts, n_atoms, protein.atoms.size());
	freePTS(pts, n_atoms); 

	dt.cut_off=p_dcut;
	dt.filter_tess(protein, mole);
}

