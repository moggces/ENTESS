/*
EnTess
Authors: Jui-Hua Hsieh (moggces@gmail.com), Shuxing Zhang
Date: 04/11
Organization of the C-style EnTess codes to C++ version ; Current implementation: 554 tetrahedral types / (entess|plmct|count)
Debug: -protein=2VT4_glide-pocket.pdb -ligand=first.mol2 -path=C:\Users\Jui-Hua\Documents\Project\4Tang -mode=plmct
-list=tmp.list -mode=plmct -path=C:\Users\Jui-Hua\Documents\Project\CSAR\NRC\set1
-protein=241.pdb -ligand=241_multi.mol2 -mode=plmct -path=C:\Users\Jui-Hua\Documents\Project\CSAR\NRC\set1
*/
#include "EnTess.h"
#include <exception>
using namespace std;


void usage() {
   cout << "Usage: ./EnTess -protein=aaa.pdb -ligand=bbb.mol2 -path=? [-dcut=8(d) -mode=entess(d)|plmct|count]" << endl;
   cout <<"\t[bbb.mol2 could be a multi-mol2 file]" << endl;
   cout << "\tOR" << endl;
   cout << "\t-list=ccc.list -path=? [-dcut=8(d) -mode=entess(d)plmct|count]]" <<endl;
   cout << "\t[colum1: aaa.pdb; column2: bbb.mol2; bbb.mol2 is a single-mol2 file]" << endl;
   exit(1);
}

// will have problems if user forgot to put argument at -protein ...
void getParameters(int argc, char * argv[])
{
	if (argc < 2)
	{
		usage();
	} else 
	{ 
		for ( int i=0; i< argc; ++i)
		{
			cout <<  argv[i] << " ";
		}
		cout << endl;
	}

	for (int i=1; i < argc; ++i)
	{
		if (argv[i][0] == '-')
		{
			if (memcmp(argv[i], "-protein=", 9) == 0)
			{
				A_STRING temp (argv[i]);
				p_pdb = temp.substr(9, temp.length()-9);
				continue;
			} else if ( memcmp(argv[i], "-ligand=", 8)== 0 ) 
			{
				A_STRING temp(argv[i]);
				p_mol2 = temp.substr(8, temp.length()-8);
				continue;
			} else if (memcmp(argv[i], "-list=", 6) == 0) 
			{
				A_STRING temp(argv[i]);
				p_list = temp.substr(6, temp.length()-6);
				continue;
			} else if (memcmp(argv[i], "-path=", 6) == 0 )
			{
				A_STRING temp(argv[i]);
				p_path = temp.substr(6, temp.length()-6);
				continue;
			} else if (memcmp(argv[i], "-mode=", 6) == 0)
			{
				A_STRING temp(argv[i]);
				p_mode = temp.substr(6, temp.length()-6);
				
				if ( p_mode != "entess" && p_mode != "plmct" && p_mode != "count" && p_mode != "volume") 
				{
					cout << p_mode << " is not implemented" << endl;
					exit(1);
				}
				cout << "mode type is changed to: "<< p_mode << endl;
				
				continue;
			} else if ( memcmp(argv[i], "-dcut=", 6) == 0 )
			{
				A_STRING  temp(argv[i]);
				p_dcut = atof(temp.substr(6, temp.length()-6).c_str());
				continue;
			} else { usage(); }
		}
	}

	if (p_path.length() > 0 )
	{
		if (p_pdb.length() > 0)
		{
			if (p_mol2.length() == 0) {usage();}
			cout << "mode: " << p_mode << endl;
			cout << "Tessellation cutoff: " << p_dcut << endl;
		} else if (p_list.length() == 0) 
		{
			usage();
		} else
		{
			cout << "mode: " << p_mode<< endl;
			cout << "Tessellation cutoff: " << p_dcut << endl;
		}
	} else { usage();}
}

float ** trans2OldTess(Molecule& protein, Molecule& ligand, int& n_atoms)
{
	//pointer initialization
	float ** pts;
	try { pts = new (float *[n_atoms+4]); }
	catch (bad_alloc&) { cout << "Error allocating memory." << endl;}

	for (int i = 0; i<n_atoms+4; ++i) //+4 is necessary , could be the bugs in original tessellation program
	{
		try {pts[i] = new float[3];}
		catch (bad_alloc&) { cout << "Error allocating memory." << endl;}
		pts[i][0], pts[i][1], pts[i][2] = 0.0;
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

void freePTS(float ** pts, int n_atoms)
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

void checkPath(A_STRING& path)
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


void importList(vector<A_STRING>& proteins, vector<A_STRING>& ligands)
{
	checkPath(p_path);
	ifstream list_file ( (p_path+p_list).c_str());
	if (list_file.is_open())
	{
		A_STRING line;
		while (! list_file.eof() )
		{
			A_STRING column;
			line.getline(list_file);
			line += " ";
			stringstream ss(line.c_str());
			for ( unsigned int i=1;  ss >> column; i++)
			{
				switch (i)
				{
					case 1:
						proteins.push_back(column); break;
					case 2:
						ligands.push_back(column); break;
				}
			}
		}
	} else { cout << "can't open file: " << p_path << p_list << endl; exit(1);}
}
	

int main( int argc, char ** argv ) {

	time_t start, end;
	double diff;
	time (&start);
	//get the parameters on the commandline
	getParameters(argc, argv); // either p_list or p_mol2 && p_mol2 is set
	
	if (p_list.length() == 0)
	{
		Dataset  protein;
		protein.importPDB(p_path, p_pdb);
		Dataset  ligand;
		ligand.importMOL2_tessellate(p_path, p_mol2, p_mode, p_list, p_dcut, protein.moles[0]); 
		
	} else 
	{
		vector<A_STRING> proteins;
		vector<A_STRING> ligands;
		importList(proteins, ligands);

		if (proteins.size() != ligands.size() ) { cout << "# of proteins and # of ligands doesn't match." << endl;  exit(1);}
		else 
		{
			
			for (unsigned int i=0; i<proteins.size(); ++i)
			{
				Dataset  protein;
				Dataset  ligand;
				
				protein.importPDB(p_path, proteins[i]);
				ligand.importMOL2(p_path, ligands[i]);

				cout << "P-L complex ID: " << i+1  << " Tesselation" << endl;

				//tessellation
				float ** pts;
				int n_atoms = (protein.moles[0].atoms.size()+ ligand.moles[0].atoms.size());
				pts = trans2OldTess(protein.moles[0], ligand.moles[0], n_atoms); 
				Dtess dt;
				dt.tessellate(pts, n_atoms, protein.moles[0].atoms.size());
				freePTS(pts, n_atoms); 

				dt.cut_off=p_dcut;
				dt.filter_tess(protein.moles[0],ligand.moles[0]);
				//dt.displayDTPymol(protein.moles[0],ligand.moles[0], p_path, ligand.basename);

				//descriptor generation
				Descriptor dst(0, p_mode); // only King's type is implemented now
				//dst.setDescType(p_mode);
				dst.genDesc(dt.tetrans);
				dst.combDesc(p_path, ligand.basename, "", p_mode, p_list, true);
				
				if ( i == proteins.size() -1)
				{
					dst.printDesc(p_path, ligand.basename, p_mode, p_list);
				}
			}
		}
	}

	time (&end);
	diff = difftime (end,start);
	int hours = (int) diff/3600;
	int minutes = (int) diff%3600;
	int seconds = (int) minutes%60;
	minutes /=60;


	cout << "It took " << hours << " hours " << minutes << " minutes " << seconds << " seconds."<< endl;

	return 0;

}
