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
#include <stdio.h>
#include <fstream>
#include <string.h>

using namespace std;

#define Complex complex<double>

//Functions
void diag_mat(int nroao, double* mat, double* vals, double* vecs);
void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mat_matl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f);
void mat_mat_Tl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f);
void mat_T_matl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f);
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir); //for large matrices
void symortho_mat(int nroao, Complex* mat,  double *tmat, Complex* dummat);
void symortho_MOs(int nroao,int nroe, Complex* MOs, double* tmat, Complex* dumvec);
void symortho_MOsx(int nroao, double* MOs, double* tmat);
void pmv(double* mat, Complex* vi, Complex* vo, int nroao);
void renorm_MOs(int nroao, double* MOs);
void renorm_MO(int nroao, double* MOs, int i);
void mat_vec(int nroao, double* mat, double* vec, double* resvec);
void transpose_mat(int nroao, double* mat, double* transpose);
void mat_mat_Tn(int np, double* mat_i1, double* mat_i2, double* mat_f, int n);

//Extern Functions
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA, 
                        double* W, double* WORK, int* LWORK, int*  INFO );
extern "C" void  zheev_(char* JOBZ, char* UPLO, int* N, Complex* A, int* LDA,
			double* W, Complex* WORK, int* LWORK, double* RWORK, int* INFO);

/*******************************************************************************
 * Matrix diagonalization                                                      *
 *                                                                             *
 ******************************************************************************/

void   diag_mat(int nroao, double* mat, double* vals, double* vecs){
  static int N = -1 ;
  static double* WORK;
  
  int LDA = nroao, LWORK=3*nroao-1, INFO;
  
  if(N != nroao){
    if(WORK != NULL){
      delete [] WORK; 
     }
    
    WORK = new double[LWORK];
    N = nroao;
  }

  char JOBZ, UPLO;
  double* A    = vecs;
  double* W    = vals;
  JOBZ = 'V';
  UPLO = 'U';  

  for(int x = 0; x <  nroao*nroao; x++)
    A[x] = mat[x];
  dsyev_( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
}

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
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     SIMPLE MATRIX ROUTINES                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}


void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}

//For large matrices
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}


void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}

void mat_matl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}


void mat_mat_Tl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_matl(long long int np, long double* mat_i1, long double* mat_i2, long double* mat_f){
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(long long int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     SIMPLE MATRIX ROUTINES   (COMPLEX)                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void mat_mat(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}


void mat_mat_T(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     BASIS TRANSFORMAITIONS                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// FOR ALL ROUTINES
// DIR = 1 means:   mat = trans * mat * trans^T
// else             mat = trans^T*mat*trans


void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat 
    mat_mat( np, trans, mat, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);
  
}

// For large matrices
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat 
    mat_mat( np, trans, mat, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);
  
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     BASIS TRANSFORMAITIONS (COMPLEX)                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// FOR ALL ROUTINES
// DIR = 1 means:   mat = trans * mat * trans^T
// else             mat = trans^T*mat*trans


void trans_mat(int np, Complex* mat, Complex* trans, Complex* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat 
    mat_mat( np, trans, mat, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);
  
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

void symortho_MOsx(int nroao, double* MOs, double* tmat){
    double* dumvec = new double[nroao];
      for(int x = 0; x < nroao; x++){
        for(int y = 0; y < nroao; y++){
	  dumvec[y] = 0.;
	  for(int z = 0; z < nroao; z++){
	    dumvec[y] += tmat[y*nroao+z]*MOs[x*nroao+z];
	  }
        }
	for(int y = 0; y < nroao; y++) MOs[x*nroao+y] = dumvec[y];
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

/*******************************************************************************
 * Hermitian Matrix diagonalization                                            *
 *                                                                             *
 ******************************************************************************/

void   diag_matH(int nroao, Complex* mat, double* vals, Complex* vecs){
   static int N = -1 ;
   static Complex* WORK;
   static double*  RWORK;
  
   int LDA = nroao, LWORK=3*nroao-1, INFO;
  
   if(N != nroao){
     if(WORK != NULL){
       delete [] WORK; 
       delete [] RWORK;
     }
    
     WORK = new Complex[LWORK];
     RWORK = new double[4*nroao+2];
     N = nroao;
   }

   char JOBZ, UPLO;
   double* W    = vals;
   JOBZ = 'V';
   UPLO = 'U';  

   for(int x = 0; x <  nroao*nroao; x++)
     vecs[x] = mat[x];
   zheev_(&JOBZ, &UPLO,  &N, vecs, &LDA, W, WORK, &LWORK, RWORK, &INFO);
}


/*******************************************************************************
 * renorm_MOs: renorms molecular orbital vectors to 1                          *
 * *****************************************************************************/
// all in one step
void renorm_MOs(int nroao, double* MOs){
 double abs;
 double normfactor;
 for(int i = 0; i < nroao; i++){
  abs = 0.;
  for(int j = 0; j < nroao; j++){
   abs += MOs[i*nroao+j]*MOs[i*nroao+j];
  }
  normfactor = 1./sqrt(abs);
  for(int j = 0; j < nroao; j++){
   MOs[i*nroao+j] *= normfactor;
  }
 }
}

// only one, needed for modified Gramm-Schmidt
void renorm_MO(int nroao, double* MOs, int i){
 double abs;
 double normfactor;
  abs = 0.;
  for(int j = 0; j < nroao; j++){
   abs += MOs[i*nroao+j]*MOs[i*nroao+j];
  }
  normfactor = 1./sqrt(abs);
  for(int j = 0; j < nroao; j++){
   MOs[i*nroao+j] *= normfactor;
  }
}

void mat_vec(int nroao, double* mat, double* vec, double* resvec){
    for(int j = 0; j < nroao; j++){
	for(int i = 0; j < nroao; j++){
           resvec[j] += mat[i*nroao+j] * vec[j];
	}
    }
}

void transpose_mat(int nroao, double* mat, double* transpose){
  for(int i = 0; i < nroao; i++){
   for(int j = 0; j < nroao; j++){
     transpose[i*nroao+j] = mat[j*nroao+i];
   }
  }
}

void mat_mat_Tn(int np, double* mat_i1, double* mat_i2, double* mat_f, int n){
   for(int x = 0; x < np; x++){
//      for(int y = 0; y < np; y++){
//	if(y == n){
         mat_f[n*np+x] = 0.;
         for(int z = 0; z < np; z++){
           mat_f[n*np+x] += mat_i1[x*np+z]*mat_i2[n*np+z];
         }
//       }else{
//	 mat_f[x*np+y] = mat_i2[x*np+y]; 
//       } 
//      }
   }
}

