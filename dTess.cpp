#include "dTess.h"
using namespace std;

/**********************************************************************************
Modified from Bala's program
Programmer:  David Bostick
FUNCTION: tessellate.c
Date:  5/1/00
Description: finds the 3D delaunay triangulation of a set of points.  It is tweaked
	for use with points about 3 units apart (angstroms)   
Inputs: pts- an array full of points in 3d that you wish to tessellate
        nrs- an integer containing the number of points to be tessellated i. e. 
 	for residues, nrs is the number of residues.
	qtess- a tetran pointer sent to the function to be updated with the tetrahedra
	found by the tessellation..
	 
Outputs: the function returns an integer value equalling the number of quadruplets
	in the array quad that are returned by the function. 

 * additional information about this algorithm can be found in -
 *    CONTOURING: A guide to the analysis and display of spatial data,
 *    by David F. Watson, Pergamon Press, 1992, ISBN 0 08 040286 0
 *
************************************************************************************/

Dtess::Dtess(){
	cut_off = 8;
}
Dtess::~Dtess(){
}

float Dtess::cut_off = 8;

// original tessellate functions

int * Dtess::IntVect(int ncols)
{  int *vectptr;
   if ((vectptr = (int *) malloc(ncols * sizeof(int))) EQ NULL)
   {  fprintf(stderr,"\nnsort: Unable to allocate storage for ivector");
      exit(1);
   }
   return vectptr;
}

void Dtess::FreeVecti(int *vectptr)
{  free(vectptr);
}

int ** Dtess::IntMatrix(int nrows, int ncols)
{  int i0;
   int **matptr;
   if (nrows<2) nrows = 2;
   if (ncols<2) ncols = 2;
   if ((matptr = (int **) malloc(nrows * sizeof(int *))) EQ NULL)
   {  fprintf(stderr,"\nnsort: Unable to allocate storage for **imatrix");
      exit(1);
   }
   if ((matptr[0] = (int *) malloc(nrows * ncols * sizeof(int))) EQ NULL)
   {  fprintf(stderr,"\nnsort: Unable to allocate storage for imatrix[]");
      exit(1);
   }
   for (i0=1; i0<nrows; i0++) matptr[i0] = matptr[0] + i0 * ncols;
   return matptr;
}

void Dtess::FreeMatrixi(int **matptr)
{  free(matptr[0]);
   free(matptr);
}

double ** Dtess::DoubleMatrix(int nrows, int ncols)
{  int i0;
   double **matptr;
   if (nrows<2) nrows = 2;
   if (ncols<2) ncols = 2;
   if ((matptr = (double **) malloc(nrows * sizeof(double *))) EQ NULL)
   {  fprintf(stderr,"\nnsort: Unable to allocate storage for **dmatrix");
      exit(1);
   }
   if ((matptr[0] = (double *) malloc(nrows * ncols * sizeof(double))) EQ NULL)
   {  fprintf(stderr,"\nnsort: Unable to allocate storage for dmatrix[]");
      exit(1);
   }
   for (i0=1; i0<nrows; i0++) matptr[i0] = matptr[0] + i0 * ncols;
   return matptr;
}

void Dtess::FreeMatrixd(double **matptr)
{  free(matptr[0]);
 free(matptr);
}

void Dtess::tessellate(float **& pts, int NumOfAtoms)
{ 
   double xx, bgs, **mxy, **wrk, **ccr;
   int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i11, ii[3], 
       dim, dm, dim1, nts, tsz, chl, *id, **tmp, **a3s;
   
   
   int nsimplices = 0;
   dim = 3;
   chl = 0;

   mxy = DoubleMatrix(2, dim);

   for (i0=0; i0<dim; i0++) mxy[0][i0] = - (mxy[1][i0] = BIGNUM);
   dim1 = dim + 1;

   wrk = DoubleMatrix(dim, dim1);
   for (i0=0; i0<dim; i0++) for (i1=0; i1<dim1; i1++) wrk[i0][i1] = -RANGE;
   for (i0=0; i0<dim; i0++) wrk[i0][i0] = RANGE * (3 * dim - 1);



/******************** assigning coordinates**************************************/
   for (i0=0; i0<NumOfAtoms; i0++)
   {  
      for (i1=0; i1<dim; i1++)
      {
	 if (mxy[0][i1] < *(pts[i0]+i1)) mxy[0][i1] = *(pts[i0]+i1);
         if (mxy[1][i1] > *(pts[i0]+i1)) mxy[1][i1] = *(pts[i0]+i1);
      }
   }

   for (bgs=0, i0=0; i0<dim; i0++) 
   {  mxy[0][i0] -= mxy[1][i0];
      if (bgs < mxy[0][i0]) bgs = mxy[0][i0];
   }

   
   bgs *= EPSILON;
   srand(367);   
   for (i0=0; i0<NumOfAtoms; i0++) for (i1=0; i1<dim; i1++)
      *(pts[i0]+i1) += bgs * (0.5 - (double)rand() / 0x7fffffff); /* random numbers [0, 1] */

   for (i0=0; i0<dim1; i0++) for (i1=0; i1<dim; i1++)               /* line number 100 */
      *(pts[NumOfAtoms+i0]+i1) = mxy[1][i1] + wrk[i1][i0] * mxy[0][i1];

   for (i1=1, i0=2; i0<dim1; i0++) i1 *= i0;
   tsz = TSIZE * i1;
   tmp = IntMatrix(tsz + 1, dim);
   i1 *= (NumOfAtoms + 50 * i1)*10;  /* storage allocation - increase value of `i1' for 3D if necessary ********/
   id = IntVect(i1);
   for (i0=0; i0<i1; i0++) id[i0] = i0;
   a3s = IntMatrix(i1,dim1);
   ccr = DoubleMatrix(i1,dim1);
   for (a3s[0][0]=NumOfAtoms, i0=1; i0<dim1; i0++) a3s[0][i0] = a3s[0][i0-1] + 1;
   for (ccr[0][dim]=BIGNUM, i0=0; i0<dim; i0++) ccr[0][i0] = 0;
   nts = i4 = 1;
   dm = dim - 1;

   for (i0=0; i0<NumOfAtoms; i0++)
   {  i1 = i7 = -1;
      i9 = 0;
      for (i11=0; i11<nts; i11++)
      {  i1++;
         while (a3s[i1][0] < 0) i1++;
         xx = ccr[i1][dim];
         for (i2=0; i2<dim; i2++)
         {  xx -= SQ( *(pts[i0]+i2) - ccr[i1][i2]);
            if (xx<0) goto Corner3;
         }
         i9--;
         i4--;
         id[i4] = i1;
         for (i2=0; i2<dim1; i2++)
         {  ii[0] = 0;
            if (ii[0] EQ i2) ii[0]++;
            for (i3=1; i3<dim; i3++)
            {  ii[i3] = ii[i3-1] + 1;
               if (ii[i3] EQ i2) ii[i3]++;
            }
            if (i7>dm)
            {  i8 = i7;
               for (i3=0; i3<=i8; i3++)
               {  for (i5=0; i5<dim; i5++) if (a3s[i1][ii[i5]] NE tmp[i3][i5]) goto Corner1;
                  for (i6=0; i6<dim; i6++) tmp[i3][i6] = tmp[i8][i6];
                  i7--;
                  goto Corner2;
Corner1:;
               }
            }
            if (++i7 > tsz)
            {  fprintf(stderr,"\nnnsort: Temporary storage exceeded - increase TSIZE");
               exit(1);
            }
            for (i3=0; i3<dim; i3++) tmp[i7][i3] = a3s[i1][ii[i3]];
Corner2:;
         }
         a3s[i1][0] = -1;
Corner3:;
      }
      for (i1=0; i1<=i7; i1++)
      {  if (!(chl AND tmp[i1][0] < NumOfAtoms))
         {  for (i2=0; i2<dim; i2++)
            {  for (wrk[i2][dim]=0, i3=0; i3<dim; i3++)
               {  wrk[i2][i3] = *( pts[tmp[i1][i2]]+i3 ) - *(pts[i0]+i3);
                  wrk[i2][dim] += wrk[i2][i3] * ( *(pts[tmp[i1][i2]]+i3) + *(pts[i0]+i3) ) / 2;
               }
            }
            if (dim < 3)
            {  xx = wrk[0][0] * wrk[1][1] - wrk[1][0] * wrk[0][1];
               ccr[id[i4]][0] = (wrk[0][2] * wrk[1][1] - wrk[1][2] * wrk[0][1]) / xx;
               ccr[id[i4]][1] = (wrk[0][0] * wrk[1][2] - wrk[1][0] * wrk[0][2]) / xx;
            }
            else
            {  xx = (wrk[0][0] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                  - (wrk[0][1] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                  + (wrk[0][2] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1]));
               ccr[id[i4]][0] = ((wrk[0][3] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                               - (wrk[0][1] * (wrk[1][3] * wrk[2][2] - wrk[2][3] * wrk[1][2])) 
                               + (wrk[0][2] * (wrk[1][3] * wrk[2][1] - wrk[2][3] * wrk[1][1]))) / xx;
               ccr[id[i4]][1] = ((wrk[0][0] * (wrk[1][3] * wrk[2][2] - wrk[2][3] * wrk[1][2])) 
                               - (wrk[0][3] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                               + (wrk[0][2] * (wrk[1][0] * wrk[2][3] - wrk[2][0] * wrk[1][3]))) / xx;
               ccr[id[i4]][2] = ((wrk[0][0] * (wrk[1][1] * wrk[2][3] - wrk[2][1] * wrk[1][3])) 
                               - (wrk[0][1] * (wrk[1][0] * wrk[2][3] - wrk[2][0] * wrk[1][3])) 
                               + (wrk[0][3] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1]))) / xx;
            }
            for (ccr[id[i4]][dim]=0, i2=0; i2<dim; i2++) 
            {  ccr[id[i4]][dim] += SQ(pts[i0][i2] - ccr[id[i4]][i2]);
               a3s[id[i4]][i2] = tmp[i1][i2];
            }
            a3s[id[i4]][dim] = i0;
            i4++;
            i9++;
         }
      }
      nts += i9;
   }

   
   //cout << "Middle of the tesselate" << endl;
   
   
   if (chl)
   {  i0 = -1;
   }
/********* printing the vertices *******************/

   else
   {  i0 = -1;
      for (i11=0; i11<nts; i11++)
      {  i0++;
         while (a3s[i0][0] < 0) i0++;
         if (a3s[i0][0] < NumOfAtoms)
         {  for (i1=0; i1<dim; i1++) for (i2=0; i2<dim; i2++) 
               wrk[i1][i2] = *(pts[a3s[i0][i1]]+i2) - *(pts[a3s[i0][dim]]+i2);

              xx = ((wrk[0][0] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                  -  (wrk[0][1] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                  +  (wrk[0][2] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1])));
	      
	      //start = 0;
               if (fabs(xx) > EPSILON)
               {  
					if (xx < 0)
					{
							Tetran qtess;
							qtess.v[0] = a3s[i0][0]+1;
							qtess.v[1] = a3s[i0][1]+1;
							qtess.v[2] = a3s[i0][3]+1;
							qtess.v[3] = a3s[i0][2]+1;
							qtess.sort_quad(qtess);
				
							//get_vertex_type(quad_new, atom1);
							nsimplices++;
							tetrans.push_back(qtess);
						
					} else 
					{
							Tetran qtess;
							qtess.v[0] = a3s[i0][0]+1; 
							qtess.v[1] = a3s[i0][1]+1;
							qtess.v[2] = a3s[i0][2]+1;
							qtess.v[3] = a3s[i0][3]+1;
							qtess.sort_quad(qtess);
							
							//get_vertex_type(qtess, atom1);
			
							nsimplices++;
							tetrans.push_back(qtess);
					} 
				}
			}
		 }
	 }
	

	FreeMatrixi(a3s);
	FreeMatrixd(ccr);
	FreeMatrixi(tmp);
	FreeVecti(id);

	cout << "Number of Tetrahedra (Original): " << nsimplices << endl;

}

void Dtess::tessellate(float **& pts, int NumOfAtoms, int check_point)
{ 
   double xx, bgs, **mxy, **wrk, **ccr;
   int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i11, ii[3], 
       dim, dm, dim1, nts, tsz, chl, *id, **tmp, **a3s;
   
   int nsimplices = 0;
   dim = 3;
   chl = 0;

   mxy = DoubleMatrix(2, dim);

   for (i0=0; i0<dim; i0++) mxy[0][i0] = - (mxy[1][i0] = BIGNUM);
   dim1 = dim + 1;

   wrk = DoubleMatrix(dim, dim1);
   for (i0=0; i0<dim; i0++) for (i1=0; i1<dim1; i1++) wrk[i0][i1] = -RANGE;
   for (i0=0; i0<dim; i0++) wrk[i0][i0] = RANGE * (3 * dim - 1);



/******************** assigning coordinates**************************************/
   for (i0=0; i0<NumOfAtoms; i0++)
   {  
      for (i1=0; i1<dim; i1++)
      {
	 if (mxy[0][i1] < *(pts[i0]+i1)) mxy[0][i1] = *(pts[i0]+i1);
         if (mxy[1][i1] > *(pts[i0]+i1)) mxy[1][i1] = *(pts[i0]+i1);
      }
   }

   for (bgs=0, i0=0; i0<dim; i0++) 
   {  mxy[0][i0] -= mxy[1][i0];
      if (bgs < mxy[0][i0]) bgs = mxy[0][i0];
   }

   
   bgs *= EPSILON;
   srand(367);   
   for (i0=0; i0<NumOfAtoms; i0++) for (i1=0; i1<dim; i1++)
      *(pts[i0]+i1) += bgs * (0.5 - (double)rand() / 0x7fffffff); /* random numbers [0, 1] */

   for (i0=0; i0<dim1; i0++) for (i1=0; i1<dim; i1++)               /* line number 100 */
      *(pts[NumOfAtoms+i0]+i1) = mxy[1][i1] + wrk[i1][i0] * mxy[0][i1];

   for (i1=1, i0=2; i0<dim1; i0++) i1 *= i0;
   tsz = TSIZE * i1;
   tmp = IntMatrix(tsz + 1, dim);
   i1 *= (NumOfAtoms + 50 * i1)*10;  /* storage allocation - increase value of `i1' for 3D if necessary ********/
   id = IntVect(i1);
   for (i0=0; i0<i1; i0++) id[i0] = i0;
   a3s = IntMatrix(i1,dim1);
   ccr = DoubleMatrix(i1,dim1);
   for (a3s[0][0]=NumOfAtoms, i0=1; i0<dim1; i0++) a3s[0][i0] = a3s[0][i0-1] + 1;
   for (ccr[0][dim]=BIGNUM, i0=0; i0<dim; i0++) ccr[0][i0] = 0;
   nts = i4 = 1;
   dm = dim - 1;

   for (i0=0; i0<NumOfAtoms; i0++)
   {  i1 = i7 = -1;
      i9 = 0;
      for (i11=0; i11<nts; i11++)
      {  i1++;
         while (a3s[i1][0] < 0) i1++;
         xx = ccr[i1][dim];
         for (i2=0; i2<dim; i2++)
         {  xx -= SQ( *(pts[i0]+i2) - ccr[i1][i2]);
            if (xx<0) goto Corner3;
         }
         i9--;
         i4--;
         id[i4] = i1;
         for (i2=0; i2<dim1; i2++)
         {  ii[0] = 0;
            if (ii[0] EQ i2) ii[0]++;
            for (i3=1; i3<dim; i3++)
            {  ii[i3] = ii[i3-1] + 1;
               if (ii[i3] EQ i2) ii[i3]++;
            }
            if (i7>dm)
            {  i8 = i7;
               for (i3=0; i3<=i8; i3++)
               {  for (i5=0; i5<dim; i5++) if (a3s[i1][ii[i5]] NE tmp[i3][i5]) goto Corner1;
                  for (i6=0; i6<dim; i6++) tmp[i3][i6] = tmp[i8][i6];
                  i7--;
                  goto Corner2;
Corner1:;
               }
            }
            if (++i7 > tsz)
            {  fprintf(stderr,"\nnnsort: Temporary storage exceeded - increase TSIZE");
               exit(1);
            }
            for (i3=0; i3<dim; i3++) tmp[i7][i3] = a3s[i1][ii[i3]];
Corner2:;
         }
         a3s[i1][0] = -1;
Corner3:;
      }
      for (i1=0; i1<=i7; i1++)
      {  if (!(chl AND tmp[i1][0] < NumOfAtoms))
         {  for (i2=0; i2<dim; i2++)
            {  for (wrk[i2][dim]=0, i3=0; i3<dim; i3++)
               {  wrk[i2][i3] = *( pts[tmp[i1][i2]]+i3 ) - *(pts[i0]+i3);
                  wrk[i2][dim] += wrk[i2][i3] * ( *(pts[tmp[i1][i2]]+i3) + *(pts[i0]+i3) ) / 2;
               }
            }
            if (dim < 3)
            {  xx = wrk[0][0] * wrk[1][1] - wrk[1][0] * wrk[0][1];
               ccr[id[i4]][0] = (wrk[0][2] * wrk[1][1] - wrk[1][2] * wrk[0][1]) / xx;
               ccr[id[i4]][1] = (wrk[0][0] * wrk[1][2] - wrk[1][0] * wrk[0][2]) / xx;
            }
            else
            {  xx = (wrk[0][0] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                  - (wrk[0][1] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                  + (wrk[0][2] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1]));
               ccr[id[i4]][0] = ((wrk[0][3] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                               - (wrk[0][1] * (wrk[1][3] * wrk[2][2] - wrk[2][3] * wrk[1][2])) 
                               + (wrk[0][2] * (wrk[1][3] * wrk[2][1] - wrk[2][3] * wrk[1][1]))) / xx;
               ccr[id[i4]][1] = ((wrk[0][0] * (wrk[1][3] * wrk[2][2] - wrk[2][3] * wrk[1][2])) 
                               - (wrk[0][3] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                               + (wrk[0][2] * (wrk[1][0] * wrk[2][3] - wrk[2][0] * wrk[1][3]))) / xx;
               ccr[id[i4]][2] = ((wrk[0][0] * (wrk[1][1] * wrk[2][3] - wrk[2][1] * wrk[1][3])) 
                               - (wrk[0][1] * (wrk[1][0] * wrk[2][3] - wrk[2][0] * wrk[1][3])) 
                               + (wrk[0][3] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1]))) / xx;
            }
            for (ccr[id[i4]][dim]=0, i2=0; i2<dim; i2++) 
            {  ccr[id[i4]][dim] += SQ(pts[i0][i2] - ccr[id[i4]][i2]);
               a3s[id[i4]][i2] = tmp[i1][i2];
            }
            a3s[id[i4]][dim] = i0;
            i4++;
            i9++;
         }
      }
      nts += i9;
   }

   
   //cout << "Middle of the tesselate" << endl;
   
   
   if (chl)
   {  i0 = -1;
   }
/********* printing the vertices *******************/

   else
   {  i0 = -1;
      for (i11=0; i11<nts; i11++)
      {  i0++;
         while (a3s[i0][0] < 0) i0++;
         if (a3s[i0][0] < NumOfAtoms)
         {  for (i1=0; i1<dim; i1++) for (i2=0; i2<dim; i2++) 
               wrk[i1][i2] = *(pts[a3s[i0][i1]]+i2) - *(pts[a3s[i0][dim]]+i2);

              xx = ((wrk[0][0] * (wrk[1][1] * wrk[2][2] - wrk[2][1] * wrk[1][2])) 
                  -  (wrk[0][1] * (wrk[1][0] * wrk[2][2] - wrk[2][0] * wrk[1][2])) 
                  +  (wrk[0][2] * (wrk[1][0] * wrk[2][1] - wrk[2][0] * wrk[1][1])));
	      
	      //start = 0;
               if (fabs(xx) > EPSILON)
               {  
					if (xx < 0)
					{
							Tetran qtess;
							qtess.v[0] = a3s[i0][0]+1;
							qtess.v[1] = a3s[i0][1]+1;
							qtess.v[2] = a3s[i0][3]+1;
							qtess.v[3] = a3s[i0][2]+1;
							qtess.sort_quad(qtess);
				
							if ( qtess.isPLtetran(qtess, check_point) )
							{
								nsimplices++;
								tetrans.push_back(qtess);
							}
						
					} else 
					{
							Tetran qtess;
							qtess.v[0] = a3s[i0][0]+1; 
							qtess.v[1] = a3s[i0][1]+1;
							qtess.v[2] = a3s[i0][2]+1;
							qtess.v[3] = a3s[i0][3]+1;
							qtess.sort_quad(qtess);
							
							if ( qtess.isPLtetran(qtess, check_point) )
							{
								nsimplices++;
								tetrans.push_back(qtess);
							}
					} 
				}
			}
		 }
	 }
	

	FreeMatrixi(a3s);
	FreeMatrixd(ccr);
	FreeMatrixi(tmp);
	FreeVecti(id);

	cout << "Number of P-L Tetrahedra (Original): " << nsimplices << endl;

}
/*
Based on Bala's program
*/
void Dtess::filter_tess(Molecule& protein, Molecule& ligand)
{
   cut_off=cut_off*cut_off;
   int n_p_atoms = protein.atoms.size();
/*
   vector<Tetran>::iterator it;
   for (it=tetrans.begin(); it<tetrans.end(); )
   {
		   for (int k=0; k< 4; ++k)
		   {
			   if ((*it).v[k]-1 < n_p_atoms)
			   {
				   (*it).nodes[k] = protein.atoms[(*it).v[k]-1];
			   } else { (*it).nodes[k] = ligand.atoms[(*it).v[k]-1-n_p_atoms];}
		   }

		   if ( it->get_Tetran_edges_distsq(*it, cut_off) != 0 ) 
		   { 
			   it = tetrans.erase(it);
		   } else {
			   ++it;
		   }
		
   }
*/

   
  int j = 0;
  for (unsigned int i = 0; i < tetrans.size(); ++i) 
  {
	   for (int k=0; k< 4; ++k)
	   {
		   if (tetrans[i].v[k]-1 < n_p_atoms)
		   {
				   tetrans[i].nodes[k] = protein.atoms[tetrans[i].v[k]-1];
			} else { tetrans[i].nodes[k] = ligand.atoms[tetrans[i].v[k]-1-n_p_atoms];}
		}

		if ( tetrans[i].get_Tetran_edges_distsq(tetrans[i], cut_off) == 0) 
		{
			tetrans[j++] = tetrans[i];
		}
  }
  tetrans.resize(j);

  cout << "Number of tetrahedra remain after cutoff filtering:  " << tetrans.size() << endl;
   
}

void Dtess::displayDTPymol(Molecule& protein, Molecule& ligand, A_STRING& path, A_STRING& base)
{
	ofstream pymolfile;
	pymolfile.open((path+base+".py").c_str());	
	
	if (pymolfile.is_open())
	{
		pymolfile << "from pymol import cmd" << endl;
		pymolfile << "cmd.set(\"auto_show_spheres\", \"on\")" << endl;
		pymolfile << "cmd.set(\"sphere_scale\", .25)" << endl;
		pymolfile << "cmd.set(\"auto_show_lines\", \"on\")" << endl;

		 
		int n_atoms, break_point = 0;
		n_atoms = (protein.atoms.size() + ligand.atoms.size());
		break_point = protein.atoms.size();

		int *ptemp = new int [n_atoms]; 
		for (int c = 0 ; c<n_atoms; ++c) { ptemp[c] = 0;}
		
		set<A_STRING> inter_edges; // however, the ligand-ligand interaction will also be included
		set<A_STRING> intra_edges;

		for (unsigned int i =0; i<tetrans.size(); ++i)
		{
			for (int k = 0; k<4; ++k)
			{
				//get unique atom ID ...
				ptemp[tetrans[i].v[k]-1] = ptemp[tetrans[i].v[k]-1]+1;
			}

			int b=0;
			for (int a = 0; a< 3; ++a)
			{
				++b;
				for (int e=b; e<4; ++e)
				{
					//get unique bond
					//char ft [6];
					//char sd [6];
					//sprintf_s(ft, "%d", tetrans[i].v[a]-1);
					//sprintf_s(sd, "%d", tetrans[i].v[e]-1);
					stringstream ft;
					stringstream sd;
					ft << (tetrans[i].v[a]-1);
					sd << tetrans[i].v[e]-1;
					A_STRING key = A_STRING(ft.str().c_str()) + "-" + A_STRING(sd.str().c_str());
					
					if ( tetrans[i].v[e]-1 >= break_point)
					{
						inter_edges.insert(key);
					} else { intra_edges.insert(key);}
				}
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
		pymolfile << "cmd.select(\"c_pro\",\"/" << base << "-tess/pro///c*\")" << endl; 
		pymolfile << "cmd.color(\"tv_green\",\"c_pro\")" << endl; 
		pymolfile << "cmd.select(\"n_pro\",\"/" << base << "-tess/pro///n*\")" << endl; 
		pymolfile << "cmd.color(\"tv_blue\",\"n_pro\")" << endl; 
		pymolfile << "cmd.select(\"o_pro\",\"/" << base << "-tess/pro///o*\")" << endl; 
		pymolfile << "cmd.color(\"tv_red\",\"o_pro\")" << endl; 
		pymolfile << "cmd.select(\"s_pro\",\"/" << base << "-tess/pro///s*\")" << endl; 
		pymolfile << "cmd.color(\"tv_yellow\",\"s_pro\")" << endl; 

		delete [] ptemp;


		
		// print inter_edges
		set<A_STRING>::iterator it;
		for (it=inter_edges.begin(); it!=inter_edges.end(); it++)
		{
			int ft =  atoi((*it).substr(0, (*it).find("-")).c_str());
			int sd =  atoi((*it).substr((*it).find("-")+1, (*it).length()-(*it).find("-")-1).c_str());
			if ( ft < break_point && sd >= break_point )
			{
				pymolfile << "cmd.bond(\"/" << base << "-tess/pro/" << protein.atoms[ft].chain << "/" << protein.atoms[ft].res_num << "/" << 
					protein.atoms[ft].atom_name << "\",\"/" << base << "-tess/lig/" << ligand.atoms[sd-break_point].chain << "/" << 
					ligand.atoms[sd-break_point].res_num << "/" << ligand.atoms[sd-break_point].atom_name << "\")" << endl; 
			} else if ( ft >=break_point && sd >= break_point) 
			{
				pymolfile << "cmd.bond(\"/" << base << "-tess/lig/" << ligand.atoms[ft-break_point].chain << "/" << ligand.atoms[ft-break_point].res_num << "/" << 
					ligand.atoms[ft-break_point].atom_name << "\",\"/" << base << "-tess/lig/" << ligand.atoms[sd-break_point].chain << "/" << 
					ligand.atoms[sd-break_point].res_num << "/" << ligand.atoms[sd-break_point].atom_name << "\")" << endl; 
			}
		}
		pymolfile << "cmd.set_bond(\"line_color\",\"forest\",\"" << base << "-tess\")" << endl; 
		pymolfile << "cmd.set_bond(\"line_width\",0.25,\"" << base << "-tess\")" << endl; 

		// print intra_edges
		for (it=intra_edges.begin(); it!=intra_edges.end(); it++)
		{
			int ft =  atoi((*it).substr(0, (*it).find("-")).c_str());
			int sd =  atoi((*it).substr((*it).find("-")+1, (*it).length()-(*it).find("-")-1).c_str());

			if ( ft < break_point && sd < break_point)
			{
				pymolfile << "cmd.bond(\"/" << base << "-tess/pro/" << protein.atoms[ft].chain << "/" << protein.atoms[ft].res_num << "/" << 
					protein.atoms[ft].atom_name << "\",\"/" << base << "-tess/pro/" << protein.atoms[sd].chain << "/" << 
					protein.atoms[sd].res_num << "/" << protein.atoms[sd].atom_name << "\")" << endl; 
			}  
		}

		pymolfile << "cmd.orient(\"" << base << "-tess\")" << endl; 
		
		pymolfile.close();
		
	} else { cout << "can't open file"  << path << base << ".py" << endl; exit(1);}

}


