/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_cis.cc                                                             *
 *                                                                              *
 * contains functions for especially needed for cis calculations                * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 *                                                    Stefan Klinkusch 2015     *
 ********************************************************************************/ 

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;


//Functions
double calc_1p_op_cis_od(int a, int b, double* MOs, int nroao, double* mat, 
			 double* tmpvec);
double calc_1p_op_cis_d(int i, int f, double* MOs, int nroao, double* opmat, 
			double* Pmat, int nroe);
void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
                     double* cismat_x, double* cismat_y, double* cismat_z, double* mumat_x, 
		     double* mumat_y, double* mumat_z, double* mucore, 
                     double* MOs, double* Pmat, double* tmpvec_ao, ofstream* outf);
void calc_mu_mat_cis_crs(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
		     double *cismat, double* mumat, double mucore, 
		     double* MOs, double* Pmat, double* tmpvec_ao, ofstream* outf);

//Extern Functions
extern double   calc_op_1el(int nroao, double* opmat, double* Pmat);

/*******************************************************************************
 *                                                                             *
 * calc_1p_op_cis_od                                                           *
 *                                                                             *
 *                                                                             *
 * one particle operator between cis basis functions   (off doagonal elements) *
 *******************************************************************************/

double calc_1p_op_cis_od(int a, int b, double* MOs, int nroao, double* mat, 
			 double* tmpvec){
  double op = 0.;
  
  for(int x = 0; x < nroao; x++){
    tmpvec[x] = 0.;
    for(int y = 0; y < nroao; y++)
      tmpvec[x] += mat[x*nroao+y] * MOs[b*nroao+y];
    op += tmpvec[x]*MOs[a*nroao+x];
  }

  return(op);
}


/*******************************************************************************
 *                                                                             *
 * calc_1p_op_cis_d                                                            *
 *                                                                             *
 *                                                                             *
 * one particle operator between cis basis functions  (diagonal elements)      *
 *******************************************************************************/

double calc_1p_op_cis_d(int i, int f, double* MOs, int nroao, double* opmat, 
			 double* Pmat, int nroe){
  for(int x = 0; x <nroao*nroao; x++) Pmat[x] = 0.;
  
  for(int e = 0; e < nroe/2; e++){
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++)
	Pmat[x*nroao+y] += 2.*MOs[e*nroao+x]*MOs[e*nroao+y];
    }
  }

  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++)
      Pmat[x*nroao+y] +=  -MOs[i*nroao+x]*MOs[i*nroao+y]+MOs[f*nroao+x]*MOs[f*nroao+y];
  }
  
  double opval = calc_op_1el( nroao,  opmat,  Pmat);
  return(opval);
}


/*******************************************************************************
 *                                                                             *
 * calc_mu_mat_cis                                                             *
 *                                                                             *
 *                                                                             *
 * calculates dipole matrix in cis_basis                                       *
 *******************************************************************************/

void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
                     double* cismat_x, double* cismat_y, double* cismat_z, double* mumat_x, 
		     double* mumat_y, double* mumat_z, double* mucore, 
                     double* MOs, double* Pmat, double* tmpvec_ao, ofstream* outf){

  
  long long int Cis_size = cis_size;
  for(long long int x = 0; x < Cis_size*Cis_size; x++){
      cismat_x[x] = 0.;
      cismat_y[x] = 0.;
      cismat_z[x] = 0.;
  }
  
  cismat_x[0] = -calc_1p_op_cis_d( 0, 0,  MOs,  nroao,  mumat_x, Pmat, nroe)+mucore[0];
  cismat_y[0] = -calc_1p_op_cis_d( 0, 0,  MOs,  nroao,  mumat_y, Pmat, nroe)+mucore[1];
  cismat_z[0] = -calc_1p_op_cis_d( 0, 0,  MOs,  nroao,  mumat_z, Pmat, nroe)+mucore[2];
  
  for(long long int x = 1; x < Cis_size; x++){
    long long int i = (x-1)/umo+llim;
    long long int f = (x-1)%umo+omo+llim;
    cismat_x[x] =  sqrt(2.)*calc_1p_op_cis_od( i, f, MOs,  nroao, mumat_x, tmpvec_ao);
    cismat_y[x] =  sqrt(2.)*calc_1p_op_cis_od( i, f, MOs,  nroao, mumat_y, tmpvec_ao);
    cismat_z[x] =  sqrt(2.)*calc_1p_op_cis_od( i, f, MOs,  nroao, mumat_z, tmpvec_ao);
  }
  long long int i1 = 0;
  long long int i2 = 0;
  long long int f1 = 0;
  long long int f2 = 0;
  for(long long int x = 1 ; x < Cis_size; x++){
    i1 = (x-1)/umo+llim;
    f1 = (x-1)%umo+omo+llim;
    *outf << x << "\t";
    if(x%10 == 9) *outf << "\n";
    outf->flush();
#pragma omp parallel for default(shared) private(i2,f2) 
    for(long long int y = x ; y < Cis_size; y++){
      i2 = (y-1)/umo+llim;
      f2 = (y-1)%umo+omo+llim;
      if(i1==i2 && f1!=f2){
        cismat_x[x*Cis_size+y] =  -calc_1p_op_cis_od( f2, f1, MOs,  nroao, mumat_x, tmpvec_ao);
        cismat_y[x*Cis_size+y] =  -calc_1p_op_cis_od( f2, f1, MOs,  nroao, mumat_y, tmpvec_ao);
        cismat_z[x*Cis_size+y] =  -calc_1p_op_cis_od( f2, f1, MOs,  nroao, mumat_z, tmpvec_ao);
      }
      if(i1!=i2 && f1==f2){
        cismat_x[x*Cis_size+y] =   calc_1p_op_cis_od( i1, i2, MOs,  nroao, mumat_x, tmpvec_ao);
        cismat_y[x*Cis_size+y] =   calc_1p_op_cis_od( i1, i2, MOs,  nroao, mumat_y, tmpvec_ao);
        cismat_z[x*Cis_size+y] =   calc_1p_op_cis_od( i1, i2, MOs,  nroao, mumat_z, tmpvec_ao);
      }
      if(y == x){
        cismat_x[x*Cis_size+y] =  -calc_1p_op_cis_d( i1, f1,  MOs,  nroao,  mumat_x, Pmat, nroe)+mucore[0];   
        cismat_y[x*Cis_size+y] =  -calc_1p_op_cis_d( i1, f1,  MOs,  nroao,  mumat_y, Pmat, nroe)+mucore[1];   
        cismat_z[x*Cis_size+y] =  -calc_1p_op_cis_d( i1, f1,  MOs,  nroao,  mumat_z, Pmat, nroe)+mucore[2];   
      }
    }
  }

  for(long long int x = 0; x < Cis_size; x++){
    for(long long int y = x; y < Cis_size; y++){
      cismat_x[y*Cis_size+x] = cismat_x[x*Cis_size+y];
      cismat_y[y*Cis_size+x] = cismat_y[x*Cis_size+y];
      cismat_z[y*Cis_size+x] = cismat_z[x*Cis_size+y];
    }
  }
}

void calc_mu_mat_cis_crs(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
		     double *cismat, double* mumat, double mucore, 
		     double* MOs, double* Pmat, double* tmpvec_ao, double* dip_val, ofstream* outf){
    int z = 0;
    int f1, f2, i1, i2;
    for(int x = 0; x < cis_size; x++){
	i1 = (x-1)/umo+llim;
	f1 = (x-1)%umo+omo+llim;
	for(int y = 0; y < cis_size; y++){
	    i2 = (y-1)/umo+llim;
	    f2 = (y-1)%umo+omo+llim;
	    if(x == 0 && y == 0){
		cismat[x*cis_size+y] = -calc_1p_op_cis_d(0, 0, MOs, nroao, mumat, Pmat, nroe) + mucore;
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }else if(x == 0){
		cismat[x*cis_size+y] = sqrt(2.)*calc_1p_op_cis_od(i2, f2, MOs, nroao, mumat, tmpvec_ao);
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }else if(y == 0){
		cismat[x*cis_size+y] = sqrt(2.)*calc_1p_op_cis_od(i1, f1, MOs, nroao, mumat, tmpvec_ao);
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }else if(x == y){
		cismat[x*cis_size+y] = -calc_1p_op_cis_d(i1, f1, MOs, nroao, mumat, Pmat, nroe) + mucore;
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }else if(i1 == i2 && f1 != f2){
		cismat[x*cis_size+y] = -calc_1p_op_cis_od(f2, f1, MOs, nroao, mumat, tmpvec_ao);
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }else if(i1 != i2 && f1 == f2){
		cismat[x*cis_size+y] = calc_1p_op_cis_od(i1, i2, MOs, nroao, mumat, tmpvec_ao);
		dip_val[z] = cismat[x*cis_size+y];
		z++;
	    }
	}
    }
}

