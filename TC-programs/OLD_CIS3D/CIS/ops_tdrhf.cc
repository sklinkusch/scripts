/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_tdrhf.cc                                                           *
 *                                                                              *
 * contains functions for  time depentend   operations  needed for RHF          *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/ 
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

using namespace std;
#define Complex complex<double>

//Functions
// double calc_Sp12(int nroao, double* Smat, double* Som2, double* tmpmat, 
// 		 double* tmpvecs, double* tmpvals);
void build_Pmat(int nroao, int nroe, Complex* Pmat, Complex* MOs);
void build_Fmat(int nroao, Complex* Fmat, Complex* Pmat, double* Hmat,
		double* intval, unsigned short* intnums,
		long long int* sortcount, long long int nrofint); 
Complex  calc_e_el(int nroao, Complex* Fmat, double* Hmat, Complex* Pmat);
void calc_mu(int nroao,  Complex* Pmat, 
             Complex* mu, double *mu_x, double* mu_y, double* mu_z);
void calc_MOenergies(int nroao, int nroe, Complex* Fmat,
                     Complex* MOs, Complex* MOenergien, Complex *cdumvec);
void calc_TDMOs(int nroao, int nroe,  Complex* Fmat, 
                Complex* MOs, Complex* MOs_dt, Complex* dumvec, Complex* dummat,
                double* Som12, double* So12, int nrofc);
void ruku4(int nroao, int nroe, double dt, double* Hmat,  Complex* Pmat,  Complex* Fmat, 
           Complex* MOs, Complex* MOs_dt, Complex *alt_MOs, Complex *sum_MOs, 
           Complex* dumvec, Complex* dummat,
           double* Som12, double* So12, double* E_field, double *Dx, double* Dy, double* Dz,
           double* intval, unsigned short* intnums,
	   long long int* sortcount, long long int nrofint, int nrofc);
Complex calc_one_p_op(int nroao,  Complex* Pmat, double *op);
void check_id_td(double* S, Complex* MOs, Complex*  dummat, int nroao, int nroe,  double* norms);
void calc_pops(double* S, double* initMOs, Complex* MOs, Complex*  dumvec, int nroao, int nroe, double* pops);

//Extern Functions
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void symortho_mat(int nroao, Complex* mat,  double *tmat, Complex* dummat);
extern void symortho_MOs(int nroao,int nroe, Complex* MOs, double* tmat, Complex* dumvec);
extern void pmv(double* mat, Complex* vi, Complex* vo, int nroao);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_Sp12                                                        */
/*                                                                               */
/* calculate S^+1/2 for symmetric orthogonalisation                              */
/* returns mimum eigenvalue of S                                                 */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// double calc_Sp12(int nroao, double* Smat, double* So12, double* tmpmat, 
// 		double* tmpvecs, double* tmpvals){
  
//   diag_mat(nroao, Smat, tmpvals, tmpvecs);
  
//   double min_val =10.;
//   for(int x = 0; x < nroao; x++){
//     if(min_val > tmpvals[x]) min_val = tmpvals[x];
//   }
  
//   for(int x = 0; x <  nroao*nroao; x++){
//     So12[x] = 0.; tmpmat[x] = 0.;
//   }
//   for(int x = 0; x <  nroao; x++)
//     So12[x*nroao+x] = sqrt(tmpvals[x]);    //1/2 !!!!!!!!!
//   for(int x = 0; x  < nroao; x++){
//     for(int y = 0; y < nroao;y++){
//       for(int z = 0; z < nroao; z++)  //Adjungiert !!!!!!
//         tmpmat[x*nroao+y] +=  So12[z*nroao+y] * tmpvecs[z*nroao+x];
//     }
//   }
//   for(int x = 0; x  < nroao; x++){
//     for(int y = 0; y < nroao;y++){
//       So12[x*nroao+y] = 0.;
//       for(int z = 0; z < nroao; z++)
//         So12[x*nroao+y] += tmpvecs[z*nroao+y] * tmpmat[x*nroao+z];
//     }
//   }  
  
//   return(min_val);
// }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*              bulid_Pmat   (for td)                                             */
/*                                                                               */
/* builds complex pmat for time dependent calculations                           */
/*                                                                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void build_Pmat(int nroao, int nroe, Complex* Pmat, Complex* MOs){
  for(int x = 0; x < nroao*nroao; x++)
    Pmat[x] = 0.;
  
  for(int e = 0; e < nroe/2;e++){
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao;y++){
        Pmat[x*nroao+y] += 2.*conj(MOs[e*nroao+y])*MOs[e*nroao+x];
      }
    }
  }
}

/*******************************************************************************
 * Permutations for fock matrix build    (td)                                  *
 *                                                                             *
 ******************************************************************************/


inline void perm_all(unsigned short a, unsigned short b, 
		     unsigned short c, unsigned short d, double integral, 
                     Complex* Fmat, Complex* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
  Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
  Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6
  Fmat[c*nroao+d] += Pmat[b*nroao+a]*integral; //7
  Fmat[d*nroao+c] += Pmat[b*nroao+a]*integral; //8
     
  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
  Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
  Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
  Fmat[c*nroao+b] -= 0.5*Pmat[d*nroao+a]*integral; //7
  Fmat[d*nroao+b] -= 0.5*Pmat[c*nroao+a]*integral; //8
}

inline void perm_1234(unsigned short a, unsigned short b, 
		      unsigned short c, unsigned short d, double integral, 
                      Complex* Fmat, Complex* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[b*nroao+a] += Pmat[c*nroao+d]*integral; //3
  Fmat[b*nroao+a] += Pmat[d*nroao+c]*integral; //4
     
  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[b*nroao+c] -= 0.5*Pmat[a*nroao+d]*integral; //3
  Fmat[b*nroao+d] -= 0.5*Pmat[a*nroao+c]*integral; //4
}


inline void perm_1256(unsigned short a, unsigned short b, 
		      unsigned short c, unsigned short d, double integral, 
                      Complex* Fmat, Complex* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[a*nroao+b] += Pmat[d*nroao+c]*integral; //2
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
  Fmat[d*nroao+c] += Pmat[a*nroao+b]*integral; //6

  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[a*nroao+d] -= 0.5*Pmat[b*nroao+c]*integral; //2
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
  Fmat[d*nroao+a] -= 0.5*Pmat[c*nroao+b]*integral; //6
}

inline void perm_15(unsigned short a, unsigned short b, 
		    unsigned short c, unsigned short d, double integral, 
                    Complex* Fmat, Complex* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
  Fmat[c*nroao+d] += Pmat[a*nroao+b]*integral; //5
     
  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
  Fmat[c*nroao+a] -= 0.5*Pmat[d*nroao+b]*integral; //5
}

inline void perm_1(unsigned short a, unsigned short b, 
		   unsigned short c, unsigned short d, double integral, 
                   Complex* Fmat, Complex* Pmat, int nroao){
  //J-Terme
  Fmat[a*nroao+b] += Pmat[c*nroao+d]*integral; //1
     
  //K-Terme
  Fmat[a*nroao+c] -= 0.5*Pmat[b*nroao+d]*integral; //1
}

/*******************************************************************************
 * Fock matrix build        (td)                                               *
 *                                                                             *
 ******************************************************************************/

void build_Fmat(int nroao, Complex* Fmat, Complex* Pmat, double* Hmat,
		double* intval, unsigned short* intnums,
		long long int* sortcount, long long int nrofint){
  for(int mu = 0; mu < nroao; mu++){
    for(int nu = 0; nu < nroao;nu++){
      Fmat[mu*nroao+nu] = Hmat[mu*nroao+nu];
    }
  }


  //PERM_1
  for(long long int x = 0; x < sortcount[0]; x++)
    perm_1(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3], 
	   intval[x], Fmat, Pmat, nroao);
    

  //PERM_15
  for(long long int x = sortcount[0]; x < sortcount[1]; x++)
    perm_15(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3], 
	    intval[x], Fmat, Pmat, nroao);

  //PERM_1234
  for(long long int x = sortcount[1]; x < sortcount[2]; x++)
    perm_1234(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3], 
	      intval[x], Fmat, Pmat, nroao);
  

  //PERM_1256
  for(long long int x = sortcount[2]; x < sortcount[3]; x++)  
    perm_1256(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3], 
	      intval[x], Fmat, Pmat, nroao);
  
  
  //PERM_ALL
  for(long long int x = sortcount[3]; x < nrofint; x++)
    perm_all(intnums[x*4+0],  intnums[x*4+1],  intnums[x*4+2],  intnums[x*4+3], 
	      intval[x], Fmat, Pmat, nroao);
}

/*******************************************************************************
 * Calc elecrtonic energy   (td)                                               *
 *                                                                             *
 ******************************************************************************/

Complex  calc_e_el(int nroao, Complex* Fmat, double* Hmat, Complex* Pmat){
  Complex E_tot_el = 0.;
  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao;y++){
      E_tot_el += Pmat[x*nroao+y]*(Fmat[y*nroao+x] + Hmat[y*nroao+x]);
    }
  }
  E_tot_el /= 2.;

  return(E_tot_el);
}

/*******************************************************************************
 * Calc mu   (td)                                                              *
 *                                                                             *
 ******************************************************************************/

void calc_mu(int nroao,  Complex* Pmat, 
             Complex* mu, double *mu_x, double* mu_y, double* mu_z){
  for(int x = 0; x < 3; x++)
    mu[x] = 0.;
  for(int x = 0; x < nroao*nroao; x++){
    mu[0] += Pmat[x]*mu_x[x];
    mu[1] += Pmat[x]*mu_y[x];
    mu[2] += Pmat[x]*mu_z[x];
  }    
}


/*******************************************************************************
 * Calc mu   (td)                                                              *
 *                                                                             *
 ******************************************************************************/
Complex calc_one_p_op(int nroao,  Complex* Pmat, double *op){

  Complex val = 0.;
  for(int x = 0; x < nroao*nroao; x++)
    val += Pmat[x]*op[x];

  return(val);
}

/*******************************************************************************
 * Calc MO energies (td)                                                       *
 *                                                                             *
 ******************************************************************************/

void calc_MOenergies(int nroao, int nroe, Complex* Fmat,
                     Complex* MOs, Complex* MOenergien, Complex *cdumvec){
  for(int e = 0; e < nroe/2; e++){
    MOenergien[e] = 0.;
    for(int x = 0; x < nroao; x++){
      cdumvec[x] = 0.;
      for(int y = 0; y < nroao; y++){
        cdumvec[x] += Fmat[x*nroao+y]*MOs[e*nroao+y];
      }
    }
    for(int x = 0; x < nroao; x++)
      MOenergien[e] += conj(MOs[e*nroao+x])*cdumvec[x];
  }
}


/*******************************************************************************
 * Calc MO time derivative of MOs                                              *
 *                                                                             *
 ******************************************************************************/

void calc_TDMOs(int nroao, int nroe,  Complex* Fmat, 
                Complex* MOs, Complex* MOs_dt, Complex* dumvec, Complex* dummat,
                double* Som12, double* So12, int nrofc){
  symortho_mat(nroao,  Fmat, Som12, dummat);
  symortho_MOs(nroao, nroe,  MOs, So12, dumvec );
  for(int e = nrofc; e < nroe/2; e++){
    for(int x = 0; x < nroao; x++)
      MOs_dt[e*nroao+x] = 0;
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++)
        MOs_dt[e*nroao+x] += Fmat[x*nroao+y]*MOs[e*nroao+y];
      MOs_dt[e*nroao+x] *= Complex(0.,1.);
    }
  }
  symortho_MOs(nroao, nroe, MOs,    Som12, dumvec );
  symortho_MOs(nroao, nroe, MOs_dt, Som12, dumvec );
}

/*******************************************************************************
 * ruku4 propagator for TDSCF                                                  *
 *                                                                             *
 ******************************************************************************/

void ruku4(int nroao, int nroe, double dt, double* Hmat,  Complex* Pmat,  Complex* Fmat, 
           Complex* MOs, Complex* MOs_dt, Complex *alt_MOs, Complex *sum_MOs, 
           Complex* dumvec, Complex* dummat,
           double* Som12, double* So12, double* E_field, double *Dx, double* Dy, double* Dz,
           double* intval, unsigned short* intnums,
	   long long int* sortcount, long long int nrofint, int nrofc){
  
  static double dt2 = dt/2.;
  static double dt3 = dt/3.;
  static double dt6 = dt/6.;
  static int msc     = nroao*nroe*8;
  
  //0 Initialisrerung
  memcpy(  alt_MOs,  MOs, (size_t) msc );
  memcpy(  sum_MOs,  MOs, (size_t) msc );
  
 
  /*----------------------------------------------
   * 1) Benutze Steigung am linken Rand 
   *    Erste Vorhersage f"ur Mittelpunkt 
   *    Aufaddieren in QSUM,PSUM (mit Faktor h/6)
   *----------------------------------------------*/

  //Build actual fock matrix
  build_Pmat(nroao, nroe, Pmat, MOs); 
  build_Fmat( nroao,  Fmat,  Pmat,  Hmat, intval,  intnums, sortcount,  nrofint);

  //ADD Field
  for(int x = 0; x < nroao*nroao; x++) 
    Fmat[x] += Dx[x]*E_field[0*3+0] + Dy[x]*E_field[0*3+1] + Dz[x]*E_field[0*3+2];

  //Calc derivative 
  calc_TDMOs(nroao,  nroe,  Fmat,  MOs, MOs_dt, dumvec, dummat, Som12, So12, nrofc);

  for(int x = nrofc*nroao; x < nroe/2*nroao; x++){
    MOs[x]      +=  MOs_dt[x]*dt2;
    sum_MOs[x]  +=  MOs_dt[x]*dt6;
  }


  /*----------------------------------------------
   * 2) Benutze (erste) Steigung am Mittelpunkt
   *    Zweite Vorhersage f"ur Mittelpunkt 
   *    Aufaddieren in QSUM,PSUM (mit Faktor h/3)
   *----------------------------------------------*/

  //Build actual fock matrix
  build_Pmat(nroao, nroe, Pmat, MOs);
  build_Fmat( nroao,  Fmat,  Pmat,  Hmat, intval,  intnums, sortcount,  nrofint);

  //ADD Field
  for(int x = 0; x < nroao*nroao; x++) 
    Fmat[x] += Dx[x]*E_field[1*3+0] + Dy[x]*E_field[1*3+1] + Dz[x]*E_field[1*3+2];

  //Calc derivative 
  calc_TDMOs(nroao,  nroe,  Fmat,  MOs, MOs_dt, dumvec, dummat, Som12, So12, nrofc);
  
  for(int x = nrofc*nroao; x < nroe/2*nroao; x++){
    MOs[x]       =  alt_MOs[x] + MOs_dt[x]*dt2;
    sum_MOs[x]  +=               MOs_dt[x]*dt3;
  }

  /*----------------------------------------------
   * 3) Benutze (zweite) Steigung am Mittelpunkt
   *    Vorhersage f"ur rechten Rand
   *    Aufaddieren in QSUM,PSUM (mit Faktor h/3)
   *----------------------------------------------*/
  
  //Build actual fock matrix
  build_Pmat(nroao, nroe, Pmat, MOs);
  build_Fmat( nroao,  Fmat,  Pmat,  Hmat, intval,  intnums, sortcount,  nrofint);

  //ADD Field
  for(int x = 0; x < nroao*nroao; x++) 
    Fmat[x] += Dx[x]*E_field[1*3+0] + Dy[x]*E_field[1*3+1] + Dz[x]*E_field[1*3+2];

  //Calc derivative 
  calc_TDMOs(nroao,  nroe,  Fmat,  MOs, MOs_dt, dumvec, dummat, Som12, So12, nrofc);

  for(int x = nrofc*nroao; x < nroe/2*nroao; x++){
    MOs[x]       =  alt_MOs[x] + MOs_dt[x]*dt;
    sum_MOs[x]  +=               MOs_dt[x]*dt3;
  }
  
  /*-----------------------------------------------
   * 4) Steigung am rechten Rand
   *    Aufaddieren in QSUM,PSUM (mit Faktor h/6)
   *-----------------------------------------------*/

  //Build actual fock matrix
  build_Pmat(nroao, nroe, Pmat, MOs);
  build_Fmat( nroao,  Fmat,  Pmat,  Hmat, intval,  intnums, sortcount,  nrofint);

  //ADD Field
  for(int x = 0; x < nroao*nroao; x++) 
    Fmat[x] += Dx[x]*E_field[2*3+0] + Dy[x]*E_field[2*3+1] + Dz[x]*E_field[2*3+2];

  //Calc derivative 
  calc_TDMOs(nroao,  nroe,  Fmat,  MOs, MOs_dt, dumvec, dummat, Som12, So12, nrofc);

  for(int x = nrofc*nroao; x < nroe/2*nroao; x++){
    sum_MOs[x]  +=               MOs_dt[x]*dt6;
  }
  
  /*-----------------------------------------------
   * Ergebnisse zurueckspeichern
   *-----------------------------------------------*/
  memcpy( MOs,  sum_MOs,   (size_t) msc );
  
}


/*******************************************************************************
 * Check identity: checks td MO-Basis                                          *
 *                                                                             *
 ******************************************************************************/

void check_id_td(double* S, Complex* MOs, Complex*  dummat, int nroao, int nroe,  double* norms){
  for(int x = 0; x < 6; x++) norms[x] = 0.;

  Complex tot_norm = 1.;
  
  //Calculate S*c_x
  for(int x = 0; x < nroe/2; x++){
    pmv(S, &(MOs[x*nroao]), &(dummat[x*nroao]), nroao);
   }
  
  //ortho
  //Calculate c_y*S*c_x  (should be delta_xy)
  for(int x = 0; x < nroe/2; x++){
    for(int y = 0; y < nroe/2; y++){
      Complex n = 0.;
      for(int z = 0; z < nroao; z++)
        n += dummat[x*nroao+z]*conj(MOs[y*nroao+z]);

      
      //check norm of MOs
      if(x==y && norm(n-1.) >  norms[2]){
	norms[2]  =  norm(n-1.);
	norms[3] +=  norm(n-1.);
	//Claculate total norm
	tot_norm *= n;
      }
      
      //check orthogonality of Mos
      if(x!=y && norm(n)    > norms[4]){
	norms[4]  =  norm(n);
	norms[5] +=  norm(n);
      }

    }
  }

  norms[0] = real(tot_norm);
  norms[1] = imag(tot_norm);
  norms[3] /= nroe/2.;
  norms[5] /= nroe/2.;

}

/*******************************************************************************
 * Calc_pops in canonical basis                                                *
 *                                                                             *
 ******************************************************************************/

void calc_pops(double* S, double* initMOs, Complex* MOs, Complex*  dumvec, int nroao, int nroe,  double* pops){
  
  for(int x = 0; x < nroao; x++) pops[x] = 0.;
  
  for(int x = 0; x < nroe/2; x++){
    pmv(S, &(MOs[x*nroao]), dumvec, nroao);
    //loop over inital MOs
    for(int y = 0; y < nroao; y++){
      Complex coeff = 0.;
      for(int z = 0; z < nroao; z++)
	coeff += initMOs[y*nroao+z]*dumvec[z];
      pops[y] += norm(coeff);
    }
  }
}
