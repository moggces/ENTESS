#include "Descriptor.h"
using namespace std;

void Descriptor::initializeOldDesc()
{
	int i, i1, i2, i3, i4;
	i = 0;
    desc_names = "";   
	
	// Generate LLLR types
	for (i1 = 1; i1 <= 6; i1++)
	{
		for (i2 = 1; i2 <= 6; i2++)
		{
			for (i3 = 1; i3 <= 6; i3++)
			{
				for (i4 = 11; i4 <= 14; i4++)
				{
					if (i1 <= i2 && i2 <= i3)
					{                   
						Desc dst;
						dst.value = 0;
						dst.tetran_type[0] = i1;
						dst.tetran_type[1] = i2;
						dst.tetran_type[2] = i3;
						dst.tetran_type[3] = i4;
						descs.push_back(dst);
						desc_names = desc_names + getDescName(dst.tetran_type) + " ";
						i++;
					} // end if
				} //end the fourth lay loop
			} // end the third loop
		} // end the second loop
	} // end the first lay loop
	
	// Generate LLRR types
	for (i1 = 1; i1 <= 6; i1++)
	{
		for (i2 = 1; i2 <= 6; i2++)
		{
			for (i3 = 11; i3 <= 14; i3++)
			{
				for (i4 = 11; i4 <= 14; i4++)
				{
					if (i1 <= i2 && i3 <= i4)
					{
						Desc dst;
						dst.value = 0;
						dst.tetran_type[0] = i1;
						dst.tetran_type[1] = i2;
						dst.tetran_type[2] = i3;
						dst.tetran_type[3] = i4;
						descs.push_back(dst);
						desc_names = desc_names + getDescName(dst.tetran_type) + " ";
						i++;
					} // end if
				} //end the fourth lay loop
			} // end the third loop
		} // end the second loop
	} // end the first lay loop
	
	// Generate LRRR types
	for (i1 = 1; i1 <= 6; i1++)
	{
		for (i2 = 11; i2 <= 14; i2++)
		{
			for (i3 = 11; i3 <= 14; i3++)
			{
				for (i4 = 11; i4 <= 14; i4++)
				{
					if (i2 <= i3 && i3 <= i4)
					{
						Desc dst;
						dst.value = 0;
						dst.tetran_type[0] = i1;
						dst.tetran_type[1] = i2;
						dst.tetran_type[2] = i3;
						dst.tetran_type[3] = i4;
						descs.push_back(dst);
						desc_names = desc_names + getDescName(dst.tetran_type) + " ";
						i++;
					} // end if
				} //end the fourth lay loop
			} // end the third loop
		} // end the second loop
	} // end the first lay loop
}

//static member
int Descriptor::charge =0;
int Descriptor::proc =0;
int Descriptor::desc_type = 0;
int Descriptor::counter = 0;
A_STRING Descriptor::desc_names ="";

Descriptor::Descriptor(int dst_type, A_STRING p_mode)
{
	desc_type = dst_type;
	if (desc_type == 0)
	{
			initializeOldDesc();
			setDescType(p_mode);
	}
	counter++;
}

Descriptor::~Descriptor()
{
}

void Descriptor::setDescType(A_STRING p_mode)
{
	//char  p_mode_f = p_mode[0];
	//char  p_mode_b = p_mode[1];
	//char * pp_mode_f = &p_mode_f;
	//char * pp_mode_b = &p_mode_b;
	//int f = atoi(pp_mode_f);
	//if ( f == NULL)
	//{
		if (p_mode == "entess")
		{
			charge = 1;
			proc = 0;
		} else if (p_mode == "plmct")
		{
			charge = 3; 
			proc = 1;
		} else if (p_mode == "count")
		{
			charge = 0;
			proc = 0;
		} else if (p_mode == "volume" )
		{
			charge = 0;
			proc = 2;
		}
		//else { cout << "mode type is changed to: "<< p_mode << endl; exit(1);}
	//} else {
	//	charge = f;
	//	proc = atoi(pp_mode_b);
	//}

	//delete pp_mode_f; 
	//delete pp_mode_b; 
}

A_STRING Descriptor::getDescName(int tetran_type [4] )
{
	A_STRING name;
	for (int i = 0; i< 4; ++i)
	{
		if( tetran_type[i] > 10) {
        if ( tetran_type[i] == 11 ) name += "Sr";
        if ( tetran_type[i] == 12 ) name += "Cr";
        if ( tetran_type[i] == 13 ) name += "Nr";
        if ( tetran_type[i] == 14 ) name += "Or";
		} else {
			if ( tetran_type[i] == 1 ) name += "Ml";
			if ( tetran_type[i] == 2 ) name += "Xl";
			if ( tetran_type[i] == 3 ) name += "Sl";
			if ( tetran_type[i] == 4 ) name += "Cl";
			if ( tetran_type[i] == 5 ) name += "Nl";
			if ( tetran_type[i] == 6 ) name += "Ol";
		}
	}

	return name;
}

void Descriptor::genDesc(vector<Tetran>& tetrans)
{

	for (unsigned int k=0; k < tetrans.size(); ++k)
	{
		if (desc_type == 0) 
		{ 
			tetrans[k].assignOLDtetranType(); 
		} // for other new descriptor types

		switch (charge)
		{
			case 0:
				break;
			case 1:
				tetrans[k].assignPaulingEN();
				break;
			case 2:
				//tetrans[k].assignDFT_EN();
				break;
			case 3:
				tetrans[k].assignDFT_MCT();
				break;
		}
	}

	for (unsigned int i = 0; i < descs.size(); ++i)
	{
		for (unsigned int k=0; k < tetrans.size(); ++k)
		{
			if (proc == 0)
			{
				if ( memcmp(descs[i].tetran_type, tetrans[k].tetran_type, sizeof(descs[i].tetran_type) ) == 0 ) // memory compare
				//if ((descs[i].tetran_type[0] == tetrans[k].tetran_type[0]) && (descs[i].tetran_type[1] == tetrans[k].tetran_type[1]) && (descs[i].tetran_type[2] == tetrans[k].tetran_type[2]) && (descs[i].tetran_type[3] == tetrans[k].tetran_type[3]))
				{
					if (charge == 0)
					{descs[i].value = descs[i].value + 1;}
					else {
						descs[i].value = descs[i].value + tetrans[k].EN[0] + tetrans[k].EN[1] + tetrans[k].EN[2] + tetrans[k].EN[3];
					} 
				}
			} else if (proc == 1)
			{
				if ( memcmp(descs[i].tetran_type, tetrans[k].tetran_type, sizeof(descs[i].tetran_type) ) == 0 ) // memory compare
				{
					int b=0;
					int d= 0;
					for (int a = 0; a< 3; ++a)
					{
						++b;
						for (int c=b; c<4; ++c)
						{
							if  ( (tetrans[k].nodes[a].isLigand + tetrans[k].nodes[c].isLigand) == 1)
							{
								descs[i].value = descs[i].value + (tetrans[k].EN[a]*tetrans[k].EN[c])/sqrt(tetrans[k].edges.dist[d++]);
							} else { d++;}
						}
					}
				}
			} else if (proc ==2 )
			{
				if ( memcmp(descs[i].tetran_type, tetrans[k].tetran_type, sizeof(descs[i].tetran_type) ) == 0 ) // memory compare
				{
					descs[i].value += tetrans[k].get_Tetran_volume();
				}
			}
		}
			if (proc == 0 && charge !=0) { descs[i].value = descs[i].value/100;} else if (proc == 1) {descs[i].value = descs[i].value*2;} 
		// not really necessary for descriptor; just be consistent with old version
	}
}




void Descriptor::combDesc(A_STRING path, A_STRING base, A_STRING pose_name, A_STRING mode, A_STRING p_list, bool isSingle)
{
	//ifstream temp_file;
	//if (p_list.length() != 0 ) { temp_file.open((path+".temp."+p_list+"."+mode).c_str());} else { temp_file.open((path+".temp."+base+"."+mode).c_str()); }
	if ( counter == 1 ) 
	{
		//temp_file.close();
		//if ( remove((path+".temp."+base+"."+mode).c_str()) != 0) { perror( "Error deleting file" ); }
		if (p_list.length() != 0 )
		{ remove((path+".temp."+p_list+"."+mode).c_str()); }  else { remove((path+".temp."+base+"."+mode).c_str());}
	}
//	temp_file.close();
	if (pose_name.length() == 0 || pose_name.find("*") == 0) {pose_name = base;}
	ofstream temp2_file;
	if (p_list.length() != 0) { temp2_file.open((path+".temp."+p_list+"."+mode).c_str(), ios::app); } else { temp2_file.open((path+".temp."+base+"."+mode).c_str(), ios::app); }
	if (temp2_file.is_open() && temp2_file.good())
	{
		temp2_file.precision(3); // the meaning of precision is different from what I understood before...
		temp2_file.setf(ios::fixed,ios::floatfield);

		if (isSingle)
		{
			temp2_file << counter << " " << base << ".pdb ";
		} else 
		{
			temp2_file << counter << " " << pose_name << "." << counter << ".pdb ";
		}
		for (unsigned int i = 0; i< descs.size(); ++i)
		{
			temp2_file << descs[i].value << " ";
		}
		temp2_file << endl;
	} else { cout << "can't open temporary file"<< endl; exit(1); } 
	temp2_file.close();

}

void Descriptor::printDesc(A_STRING path, A_STRING base, A_STRING mode, A_STRING p_list)
{
	ofstream desc_file;
	if (p_list.length() == 0)
	{
		remove((path+base + "-" + mode+".x").c_str());
		desc_file.open((path+base + "-" + mode+".x").c_str(),  ios::app);
	} else {
		remove((path+p_list + "-" + mode+".x").c_str());
		desc_file.open((path+p_list + "-" + mode+".x").c_str(),  ios::app);
	}

	if (desc_file.is_open())
	{
		desc_file << printTitle() << endl;

		ifstream temp_desc;
		if (p_list.length() == 0)
		{	temp_desc.open((path+".temp."+base+"."+mode).c_str()); } else {temp_desc.open((path+".temp."+p_list+"."+mode).c_str()); }

		if (temp_desc.is_open())
		{
			A_STRING line;
			line.getline(temp_desc);
			while ( line != "")
			{
				desc_file << line << endl;
				line.getline(temp_desc);
			}

		} else {cout << "can't open temporary file"<< endl; exit(1); } 
	} else {cout << "can't open descriptor file" << endl; exit(1); }
}

A_STRING Descriptor::printTitle()
{
	A_STRING title;
	stringstream convert;
	convert << counter;

	if (desc_type == 0)
	{
		title = A_STRING (convert.str().c_str());
		title += " 554\n";
		title += desc_names;

	}
	//delete p_counter;
	return title;

}

void Descriptor::displayOLDTetranPymol(Molecule& protein, Molecule& ligand, A_STRING path, A_STRING base, vector<Tetran>& tetrans)
{
	ofstream pymolfile;
	pymolfile.open((path+base+"-tetran.py").c_str());
	
	if (pymolfile.is_open())
	{
		int n_atoms, break_point = 0;
		n_atoms = (protein.atoms.size() + ligand.atoms.size());
		break_point = protein.atoms.size();

		pymolfile << "from pymol import cmd" << endl;
		pymolfile << "cmd.set(\"auto_show_spheres\", \"on\")" << endl;
		pymolfile << "cmd.set(\"sphere_scale\", .25)" << endl;
		pymolfile << "cmd.set(\"auto_show_lines\", \"on\")" << endl;

		//get unique atom ID ...
		int *ptemp = new int [n_atoms]; 
		for (int c = 0 ; c<n_atoms; ++c) { ptemp[c] = 0;}
		for (unsigned int i =0; i<tetrans.size(); ++i)
		{
			for (int k = 0; k<4; ++k)
			{
				ptemp[tetrans[i].v[k]-1] = ptemp[tetrans[i].v[k]-1]+1;
			}
		}
		//print nodes
		for (int c =0 ; c<n_atoms; ++c)
		{
			if (ptemp[c] >= 1)
			{
				if ( c< break_point)
				{
					pymolfile << "cmd.pseudoatom(\"" << base << "-tess\",pos=[" << protein.atoms[c].x << "," << protein.atoms[c].y << "," <<
						protein.atoms[c].z << "],segi=\"pro\"" << ",resn=\"" << protein.atoms[c].res_name << "\",resi=" << protein.atoms[c].res_num << ",chain=\"" <<
						protein.atoms[c].chain << "\",name=\"" << protein.atoms[c].atom_name << "\")" << endl;
					//pymolfile 
				} else { 
					pymolfile << "cmd.pseudoatom(\"" << base << "-tess\",pos=[" << ligand.atoms[c-break_point].x << "," << ligand.atoms[c-break_point].y << "," <<
						ligand.atoms[c-break_point].z << "], segi=\"lig\"" << ",resn=\"" << ligand.atoms[c-break_point].res_name << "\",resi=" << ligand.atoms[c-break_point].res_num << ",chain=" <<
						ligand.atoms[c-break_point].chain << ",name=\"" << ligand.atoms[c-break_point].atom_name << "\")" << endl;
				}
			}
		}
		delete [] ptemp;

		// select tetrans
		for (unsigned int i = 0; i < descs.size(); ++i)
		{
			for (unsigned int k=0; k < tetrans.size(); ++k)
			{
				if ( memcmp(descs[i].tetran_type, tetrans[k].tetran_type, sizeof(descs[i].tetran_type) ) == 0 ) // memory compare
				{
					A_STRING tetran_name = getDescName(descs[i].tetran_type);
					
					//print nodes # residue id ....
					for (unsigned int j=0; j<4; ++j)
					{
						int c = tetrans[k].v[j]-1;
						if (  c < break_point)
						{
							pymolfile << "cmd.select(\"" << tetran_name << "\",\""  << "/" << base << "-tess/pro/" << 
								protein.atoms[c].chain << "/" << protein.atoms[c].res_name << "`"<< protein.atoms[c].res_num <<
							 "/" <<  protein.atoms[c].atom_name << "\", -1, 1, 1)" << endl;
						} else {
							pymolfile << "cmd.select(\"" << tetran_name << "\",\""<< "/"<< base << "-tess/lig/" << ligand.atoms[c-break_point].chain << 
								"/" <<  ligand.atoms[c-break_point].res_name  << "`" << ligand.atoms[c-break_point].res_num << 
								"/" << ligand.atoms[c-break_point].atom_name << "\", -1, 1, 1)" << endl ;
						}
					} //for (unsigned int j=0; j<4; ++j )

					// print edges
					int b=0;
					for (int a = 0; a< 3; ++a)
					{
						++b;
						for (int e=b; e<4; ++e)
						{
							int ft = tetrans[k].v[a]-1;
							int sd = tetrans[k].v[e]-1;
							if ( ft < break_point && sd >= break_point )
							{
								pymolfile << "cmd.bond(\"/" << tetran_name << "/pro/" << protein.atoms[ft].chain << "/" << protein.atoms[ft].res_num << "/" << 
								protein.atoms[ft].atom_name << "\",\"/" << base << "-tess/lig/" << ligand.atoms[sd-break_point].chain << "/" << 
								ligand.atoms[sd-break_point].res_num << "/" << ligand.atoms[sd-break_point].atom_name << "\")" << endl; 
							} else if ( ft >=break_point && sd >= break_point) 
							{
								pymolfile << "cmd.bond(\"/" << tetran_name << "/lig/" << ligand.atoms[ft-break_point].chain << "/" << ligand.atoms[ft-break_point].res_num << "/" << 
								ligand.atoms[ft-break_point].atom_name << "\",\"/" << base << "-tess/lig/" << ligand.atoms[sd-break_point].chain << "/" << 
							ligand.atoms[sd-break_point].res_num << "/" << ligand.atoms[sd-break_point].atom_name << "\")" << endl; 
							
							} else if ( ft < break_point && sd < break_point)
							{
								pymolfile << "cmd.bond(\"/" << tetran_name << "/pro/" << protein.atoms[ft].chain << "/" << protein.atoms[ft].res_num << "/" << 
								protein.atoms[ft].atom_name << "\",\"/" << base << "-tess/pro/" << protein.atoms[sd].chain << "/" << 
								protein.atoms[sd].res_num << "/" << protein.atoms[sd].atom_name << "\")" << endl; 
							}  
						}
					} // edges done

				} //memcmp
			} //tetrans.size(); 
		} //descs.size(); 
		pymolfile << "cmd.set_bond(\"stick_radius\", 0.15, \"all\")" << endl;
		pymolfile << "cmd.orient(\"" << base << "-tess\")" << endl; 
		pymolfile.close();
	} else { cout << "can't open file"  << path << base << "-tetran.py" << endl; exit(1);}
	
}

