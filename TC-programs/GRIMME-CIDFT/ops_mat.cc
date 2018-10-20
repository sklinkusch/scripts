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
#include </home/klinkusch/bin/Libraries/mpich-3.0.4/src/include/mpi.h>

using namespace std;

#define Complex complex<double>

//Functions
void diag_mat(int nroao, double* mat, double* vals, double* vecs);
//void diag_mat_mt(int argc, char** argv, int nroao, double* mat, double* vals, double* vecs);
//void find_nps(int np, int &nprow, int & npcol);
void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir); //for large matrices
void symortho_mat(int nroao, Complex* mat,  double *tmat, Complex* dummat);
void symortho_MOs(int nroao,int nroe, Complex* MOs, double* tmat, Complex* dumvec);
void pmv(double* mat, Complex* vi, Complex* vo, int nroao);

//Extern Functions
//extern "C" void  pdsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* IA, int* JA, int* DESCA, double* W, double* Z, int* IZ, int* JZ, int* DESCZ, double* WORK, int* LWORK, int* INFO);
//extern "C" int   numroc_(int* N, int* NB, int* IPROC, int* ISRCPROC, int* NPROCS);
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA, 
                        double* W, double* WORK, int* LWORK, int*  INFO );
extern "C" void  zheev_(char* JOBZ, char* UPLO, int* N, Complex* A, int* LDA,
			double* W, Complex* WORK, int* LWORK, double* RWORK, int* INFO);
extern "C" void  dsyevd_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
	                double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
//extern "C" void  Cblacs_get(int context, int request, int* value);
//extern "C" int   Cblacs_gridinit(int* context, char* order, int nprow, int npcol);
//extern "C" void  Cblacs_gridinfo(int context, int* nprow, int* npcol, int* myrow, int* mycol);

/*******************************************************************************
 * Matrix diagonalization                                                      *
 *                                                                             *
 ******************************************************************************/

void   diag_mat(int nroao, double* mat, double* vals, double* vecs){
  static int N = -1 ;
  static double* WORK;
  
  int LDA = nroao, INFO;
  int LWORK = 3*nroao -1;
  
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
 * Matrix diagonalization multithread and inproved (ScaLAPACK)                 *
 *                                                                             *
 *******************************************************************************/

/*void diag_mat_mt(int argc, char** argv, int nroao, double* mat, double* vals, double* vecs){

    MPI::Init();
     int rank = MPI::COMM_WORLD.Get_rank();
     int nprocs = MPI::COMM_WORLD.Get_size();
    MPI_Status Mstatus;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    static int N = -1;
    static double* WORK;

    int IA = 0;
    int JA = 0;
    int IZ = 0;
    int JZ = 0;
    int LWORK = 3*nroao-1;
    int INFO;

    if(N != nroao){
        if(WORK != NULL){
	     delete [] WORK;
        }
	WORK = new double[LWORK];
	N = nroao;
    }

    char JOBZ, UPLO;
    double* Z = vecs;
    double* W = vals;
    JOBZ = 'V';
    UPLO = 'U';
    int matrix_size = nroao;
    int block = 1; 
    int myrow = 0;
    int mycol = 0;
    int izero = 0;
    int nprow = 0;
    int npcol = 0;
    int ictxt;
    char dumtext[256];
    sprintf(dumtext, "Row");
    find_nps(nprocs, nprow,npcol); 
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, dumtext, nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    int locR = numroc_(&matrix_size, &block, &myrow, &izero, &nprow);
//    int locC = numroc_(&matrix_size, &block, &mycol, &izero, &npcol);
    int* DESCA = new int[9];
    int* DESCZ = new int[9];
    DESCA[0] = 1; DESCZ[0] = 1;                      // descriptor type
    DESCA[1] = ictxt; DESCZ[1] = ictxt;             
    DESCA[2] = nroao; DESCZ[2] = nroao;              // number of rows
    DESCA[3] = nroao; DESCZ[3] = nroao;              // number of columns
    DESCA[4] = 1; DESCZ[4] = block;                  // row block size
    DESCA[5] = 1; DESCZ[5] = block;                  // column block size (must be equal to row block size)
    DESCA[6] = 0; DESCZ[6] = 0;                      // initial process row (0 by default)
    DESCA[7] = 0; DESCZ[7] = 0;                      // initial process column ( 0 by default)
    DESCA[8] = locR; DESCZ[8] = locR;                // leading dimension of local array
    double* A = new double[nroao*nroao];

    for(int x = 0; x < nroao*nroao; x++) A[x] = mat[x];
    pdsyev_(&JOBZ, &UPLO, &N, A, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ, WORK, &LWORK, &INFO);
}
*/
/*******************************************************************************
 * Find out values for nprow and npcol                                         *
 *******************************************************************************/
/*
void find_nps(int np, int &nprow, int & npcol) {

    int min_nprow=100000;
    int min_npcol=100000;


    nprow = np;
    npcol = np;


    while(1) {
     npcol--;
     if(np%2 == 0){
      if(npcol ==1){
       nprow --;
       npcol = nprow;
      }
     }else{
      if(npcol ==0){
       nprow --;
       npcol = nprow;
      }
     }

     if(nprow*npcol == np){
      min_npcol = npcol;
      if(nprow < min_nprow)    min_nprow = nprow;
     }
     if(nprow ==1 ) break;
    }
    nprow = min_nprow;
    npcol = min_npcol;
}
*/
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



