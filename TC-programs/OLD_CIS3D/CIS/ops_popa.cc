/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_popa.cc                                                            *
 *                                                                              *
 * contains ops need for td population analysis                                 * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 
#include <math.h>
#include <complex>


using namespace std;
#define Complex complex<double>
//Functions
void transformationCIS_SE(Complex* wav, Complex* SEwav, double* cisvecs, int cis_size);
void calc_MOpops(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* pops);
void calc_dens(int nroe, int nroao, int llim, int ulim, Complex* SEwav, Complex* Pmat, double* MOs);
void calc_dens_MO(int nroe, int nroao, int llim, int ulim, Complex* SEwav, Complex* Pmat);
void mat_mat(int np, Complex* mat_i1, double* mat_i2, Complex* mat_f);

//Extern Functions




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* TRANSFOMRS FROM CIS -> CSF                                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void transformationCIS_SE(Complex* wav, Complex* SEwav, double* cisvecs, int cis_size){
  for(int x = 0; x < cis_size; x++){
    SEwav[x]  = 0.;
    for(int y = 0; y < cis_size; y++)
      SEwav[x] += wav[y]*cisvecs[y*cis_size+x];
  }
}

//for large matrices
void transformationCIS_SE(Complex* wav, Complex* SEwav, double* cisvecs, long long int cis_size){
  for(long long int x = 0; x < cis_size; x++){
    SEwav[x]  = 0.;
    for(long long int y = 0; y < cis_size; y++)
      SEwav[x] += wav[y]*cisvecs[y*cis_size+x];
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* CALC POPULATIONS OF MOS                                                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_MOpops(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* pops){
  //Derived sizes
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  int cis_size = omo*umo+1;
  
  //HF GROUND STATE PART
  for(int x = 0; x < nroe/2; x++)
    pops[x] = norm(SEwav[0])*2.;
  for(int x = nroe/2; x < nroao; x++)
    pops[x] = 0.;
  
  for(int x = 1; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    double fac =norm(SEwav[x]);
    
    for(int y = 0; y < omo+llim; y++)
      pops[y] += fac*2.;
    pops[i] -= fac;
    pops[f] += fac; 
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Reconstruct density matrix  AO-SPACE                                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_dens(int nroe, int nroao, int llim, int ulim, Complex* SEwav, Complex* Pmat, double* MOs){
  for(int x = 0; x < nroao*nroao; x++)
    Pmat[x] = 0.;

  //Derived sizes
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  int cis_size = omo*umo+1;
  
  
  //Diagonal elements 
  double fac = real(conj(SEwav[0])*SEwav[0]);
  for(int x = 0; x < omo+llim; x++){
    for(int z = 0; z < nroao; z++){
      for(int Z = 0; Z < nroao; Z++)
	Pmat[z*nroao+Z] += 2.*fac*MOs[x*nroao+z]*MOs[x*nroao+Z];
    }
  }
  
  
  for(int x =1 ; x < cis_size; x++){
    double fac =  real(conj(SEwav[x])*SEwav[x]);
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    for(int x = 0; x < omo+llim; x++){
      for(int z = 0; z < nroao; z++){
	for(int Z = 0; Z < nroao; Z++)
	  Pmat[z*nroao+Z] += 2.*fac*MOs[x*nroao+z]*MOs[x*nroao+Z];
      }
    }
    for(int z = 0; z < nroao; z++){
      for(int Z = 0; Z < nroao; Z++){
	Pmat[z*nroao+Z] += fac*MOs[f*nroao+z]*MOs[f*nroao+Z];
	Pmat[z*nroao+Z] -= fac*MOs[i*nroao+z]*MOs[i*nroao+Z];
      }
    }
  }
  

  //off diagonal elements (0-other)
  for(int x = 1 ; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    Complex FAC = -sqrt(2.)*conj(SEwav[0])*SEwav[x];
    for(int z = 0; z < nroao; z++){
      for(int Z = 0; Z < nroao; Z++){
	Pmat[z*nroao+Z] +=      FAC *MOs[i*nroao+z]*MOs[f*nroao+Z];
	Pmat[z*nroao+Z] += conj(FAC)*MOs[f*nroao+z]*MOs[i*nroao+Z];
      }
    }
}
  
  //off diagonal elements (other-other)
for(int x = 1 ; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    for(int X = x+1; X < cis_size; X++){
      int I = (X-1)/umo+llim;
      int F = (X-1)%umo+omo+llim;
      Complex FAC = conj(SEwav[x])*SEwav[X];
      if(I==i){
	for(int z = 0; z < nroao; z++){
	  for(int Z = 0; Z < nroao; Z++){
	    Pmat[z*nroao+Z] +=      FAC *MOs[f*nroao+z]*MOs[F*nroao+Z];
	    Pmat[z*nroao+Z] += conj(FAC)*MOs[F*nroao+z]*MOs[f*nroao+Z];
	  }
	}
      }
      if(F==f){
	for(int z = 0; z < nroao; z++){
	  for(int Z = 0; Z < nroao; Z++){
	    Pmat[z*nroao+Z] += -conj(FAC)*MOs[i*nroao+z]*MOs[I*nroao+Z];
	    Pmat[z*nroao+Z] += -     FAC *MOs[I*nroao+z]*MOs[i*nroao+Z];
	  }
	}
      }
    }
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Reconstruct density matrix  MO-SPACE                                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_dens_MO(int nroe, int nroao, int llim, int ulim, Complex* SEwav, Complex* Pmat){
  for(int x = 0; x < nroao*nroao; x++)
    Pmat[x] = 0.;

  //Derived sizes
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  int cis_size = omo*umo+1;
  
  
  //Diagonal elements 
  double fac = real(conj(SEwav[0])*SEwav[0]);
  for(int x = 0; x < omo+llim; x++){
    Pmat[x*nroao+x] += 2.*fac;
  }
  
  
  for(int x =1 ; x < cis_size; x++){
    double fac = real(conj(SEwav[x])*SEwav[x]);
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    for(int k = 0; k < omo+llim; k++){
      Pmat[k*nroao+k] += 2.*fac;
    }
    Pmat[f*nroao+f] += fac;
    Pmat[i*nroao+i] -= fac;
  }
  
  
  // < psi_o | rho | 1pis_i^f>
  //off diagonal elements (0-other)   //LOWER DIAGONAL
  for(int x = 1 ; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    Complex FAC = -sqrt(2.)*conj(SEwav[0])*SEwav[x];
    Pmat[i*nroao+f] +=       FAC;
    //Pmat[f*nroao+i] += conj(FAC);
  }
  
 

  // <  1pis_i^f | rho | 1pis_I^F>
  //off diagonal elements (other-other)  //LOWER DIAGONAL
  for(int x = 1 ; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    for(int X = x+1; X < cis_size; X++){
      int I = (X-1)/umo+llim;
      int F = (X-1)%umo+omo+llim;
      Complex FAC = conj(SEwav[x])*SEwav[X];
      if(I==i){
	//Pmat[F*nroao+f] += conj(FAC);
	Pmat[f*nroao+F] +=  FAC;
      }
      if(F==f){
	//Pmat[I*nroao+i] += -conj(FAC);
	Pmat[i*nroao+I] +=  -conj(FAC);
      }
    }
  }
  //MAKE HERM
  for(int x = 0; x < nroao; x++){
    for(int y = x+1; y < nroao; y++){
      Pmat[y*nroao+x] = conj(Pmat[x*nroao+y]);
    }
  }

}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     SIMPLE MATRIX ROUTINES   (mixed COMPLEX double)           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void mat_mat(int np, Complex* mat_i1, double* mat_i2, Complex* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}
