/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_mat.cc                                                             *
 *                                                                              *
 * contains functions for matrix operations                                     * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <complex>

using namespace std;

#define Complex complex<double>

//Functions
void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
void symortho_mat(int nroao, Complex* mat,  double *tmat, Complex* dummat);
void symortho_MOs(int nroao,int nroe, Complex* MOs, double* tmat, Complex* dumvec);
void pmv(double* mat, Complex* vi, Complex* vo, int nroao);

//Extern Functions
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA, 
                        double* W, double* WORK, int* LWORK, int*  INFO );

/*******************************************************************************
 * matrix transformation of mat with tmat (Som12)                              *
 *                                                                             *
 *  matneu = tmat mat tmat                                                     *
 *******************************************************************************/

void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat){
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      dummat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
	dummat[x*nroao+y] += tmat[x*nroao+z] * mat[z*nroao+y];
    }
  }
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      mat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        mat[x*nroao+y] += tmat[z*nroao+y] * dummat[x*nroao+z];
    }
  }
}

/*******************************************************************************
 * Transformation of  MO's back to AO-Basis                                    *
 *                                                                             *
 ******************************************************************************/

void transform_MOs(int nroao, double *MOs, double* tmat, double* tmpvec){

  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao; y++){
      tmpvec[y] = 0.;
      for(int z = 0; z < nroao; z++){
        tmpvec[y] += tmat[y*nroao+z]*MOs[x*nroao+z];
      }
    }
    for(int y = 0; y < nroao; y++)
      MOs[x*nroao+y] = tmpvec[y];
  }  
}

/*******************************************************************************
 * transformations for symmtric orhtogonalization         (td hf)             *
 *                                                                             *
 ******************************************************************************/

void   symortho_mat(int nroao, Complex* mat,  double *tmat, Complex* dummat){
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      dummat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        dummat[x*nroao+y] += tmat[x*nroao+z] * mat[z*nroao+y];
    }
  }
  for(int x = 0; x  < nroao; x++){
    for(int y = 0; y < nroao;y++){
      mat[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++)
        mat[x*nroao+y] += tmat[z*nroao+y] * dummat[x*nroao+z];
    }
  }
}

void symortho_MOs(int nroao,int nroe, Complex* MOs, double* tmat, Complex* dumvec){
  for(int x = 0; x < nroe/2;x++){
    for(int y = 0; y < nroao; y++){
      dumvec[y] = 0.;
      for(int z = 0; z < nroao; z++){
        dumvec[y] += tmat[y*nroao+z]*MOs[x*nroao+z];
      }
    }
    for(int y = 0; y < nroao; y++)
      MOs[x*nroao+y] = dumvec[y];
  }  
}


/*******************************************************************************
 *Complex vectro matrix operations td                                          *
 *                                                                             *
 ******************************************************************************/


void pmv(double* mat, Complex* vi, Complex* vo, int nroao){
  for(int x = 0; x < nroao; x++){
    vo[x] = 0.;
    for(int y = 0; y < nroao; y++)
      vo[x] += mat[x*nroao+y]*vi[y];
  }
}

void pmv(double* mat, double* vi, double* vo, int nroao){
  for(int x = 0; x < nroao; x++){
    vo[x] = 0.;
    for(int y = 0; y < nroao; y++)
      vo[x] += mat[x*nroao+y]*vi[y];
  }
}

