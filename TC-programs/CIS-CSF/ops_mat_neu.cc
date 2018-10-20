#include <math.h>
#include <stdlib.h>
#include <iostream>

//Functions
void diag_mat(int nroao, double* mat, double* vals, double* vecs);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir); //for large matrices

//Extern Functions
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA, 
                        double* W, double* WORK, int* LWORK, int*  INFO );

//Blas
extern "C" void dsymm_(char* SIDE, char* UPLO, int* M, int* N, 
		       double* ALPHA, double* A, int* LDA, 
		       double* B, int* LDB, double* BETA, 
		       double* C, int* LDC);

extern "C" void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
		       double* ALPHA, double* A, int* LDA, double* B, int* LDB,
		       double* BETA, double* C, int* LDC);



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

/*************************************************+
  Blas 
************************************************/

void blas_smat_mat(int np, int side, double* A, double* B, double* C){
  char SIDE, UPLO;
  int M, N, LDA, LDB, LDC;
  double ALPHA, BETA;
  
  if(side == 0)
    SIDE = 'L';  //A*B
  else
    SIDE = 'R';  //B*A
  
  UPLO = 'U';
  
  M   = np;
  N   = np;
  LDA = np;
  LDB = np;
  LDC = np;
  
  ALPHA = 1.;
  BETA  = 0.;
  
  dsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}

void blas_mat_mat(int np, int transa, int transb, double* A, double* B, double* C){
  char TRANSA, TRANSB;
  int M,N,K,LDA,LDB,LDC;
  double ALPHA, BETA;
  
  if(transa == 0)
    TRANSA = 'T';
  else
    TRANSA = 'N';

  if(transb == 0)
    TRANSB = 'T';
  else
    TRANSB = 'N';

  M   = np;
  N   = np;
  K   = np;
  LDA = np;
  LDB = np;
  LDC = np;  

  ALPHA = 1.;
  BETA  = 0.;

  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA,  B, &LDB, &BETA, C, &LDC);
  
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     BASIS TRANSFORMAITIONS                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// FOR ALL ROUTINES
// DIR = 1 means:   mat = trans * mat * trans^T
// else             mat = trans^T*mat*trans


void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir){


  if(dir == 0) //CALC mat*trans
    //mat_mat( np,  mat, trans, tmpmat);
    blas_smat_mat( np, 1, mat, trans, tmpmat);
  else         //CALC trans*mat 
    //mat_mat( np, trans, mat, tmpmat);
    blas_smat_mat( np, 0, mat, trans, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    //mat_T_mat( np, trans, tmpmat, mat);
    blas_mat_mat(np, 1, 0, trans,  tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    //mat_mat_T( np, tmpmat, trans, mat);
    blas_mat_mat(np, 0, 1, tmpmat, trans, mat);
  
}

// For large matrices
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir){

  
  if(dir == 0) //CALC mat*trans
    //mat_mat( np,  mat, trans, tmpmat);
    blas_smat_mat( np, 1, mat, trans, tmpmat);
  else         //CALC trans*mat 
    //mat_mat( np, trans, mat, tmpmat);
    blas_smat_mat( np, 0, mat, trans, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    //mat_T_mat( np, trans, tmpmat, mat);
    blas_mat_mat(np, 1, 0, trans,  tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    //mat_mat_T( np, tmpmat, trans, mat);
    blas_mat_mat(np, 0, 1, tmpmat, trans, mat);
  
}


