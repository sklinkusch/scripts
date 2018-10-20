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
void transformationCIS_SE(double* wav, double* SEwav, double* cisvecs, int cis_size);
void calc_MOpops(int nroe, int nroao, int llim, int ulim, double* SEwav, double* pops);
void calc_dens(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* Pmat);

//Extern Functions




/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* TRANSFOMRS FROM CIS -> CSF                                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void transformationCIS_SE(double* wav, double* SEwav, double* cisvecs, int cis_size){
  for(int x = 0; x < cis_size; x++){
    SEwav[x]  = 0.;
    for(int y = 0; y < cis_size; y++)
      SEwav[x] += wav[y]*cisvecs[y*cis_size+x];
  }
}

//for large matrices
void transformationCIS_SE(double* wav, double* SEwav, double* cisvecs, long long int cis_size){
  for(long long int x = 0; x < cis_size; x++){
    SEwav[x]  = 0.;
    for(long long int y = 0; y < cis_size; y++)
      SEwav[x] += wav[y]*cisvecs[y*cis_size+x];
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* CALC POPULATIONS OF MOS                                                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_MOpops(int nroe, int nroao, int llim, int ulim, double* SEwav, double* pops){
  //Derived sizes
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  int cis_size = omo*umo+1;
  
  //HF GROUND STATE PART
  for(int x = 0; x < nroe/2; x++)
    pops[x] = pow(SEwav[0],2)*2.;
  for(int x = nroe/2; x < nroao; x++)
    pops[x] = 0.;
  
  for(int x = 1; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    double fac =pow(SEwav[x],2);
    
    for(int y = 0; y < omo+llim; y++)
      pops[y] += fac*2.;
    pops[i] -= fac;
    pops[f] += fac; 
  }
}

/*********************************************************************************
 * calc_CSFs                                                                     *
 * calculates Wavefunction in CSF basis from the MO populations                  *
 *********************************************************************************/



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Reconstruct density matrix                                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_dens(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* Pmat, double* MOs){
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
    double fac = real(conj(SEwav[x])*SEwav[x]);
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
    double FAC = -sqrt(2.)*(2.*real(SEwav[x])*real(SEwav[0])+ 2*imag(SEwav[x])*imag(SEwav[0]));
    for(int z = 0; z < nroao; z++){
      for(int Z = 0; Z < nroao; Z++)    
	Pmat[z*nroao+Z] += FAC*MOs[i*nroao+z]*MOs[f*nroao+Z];
    }
  }
  
  //off diagonal elements (other-other)
  for(int x = 1 ; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    for(int X = x+1; X < cis_size; X++){
      int I = (X-1)/umo+llim;
      int F = (X-1)%umo+omo+llim;
      double FAC = 2.*real(SEwav[x])*real(SEwav[X])+ 2*imag(SEwav[x])*imag(SEwav[X]);
      if(I==i){
	for(int z = 0; z < nroao; z++){
	  for(int Z = 0; Z < nroao; Z++)    
	    Pmat[z*nroao+Z] += FAC*MOs[F*nroao+z]*MOs[f*nroao+Z];
	}
      }
      if(F==f){
	for(int z = 0; z < nroao; z++){
	  for(int Z = 0; Z < nroao; Z++)    
	    Pmat[z*nroao+Z] += -FAC*MOs[I*nroao+z]*MOs[i*nroao+Z];
	}
      }
    }
  }
}
