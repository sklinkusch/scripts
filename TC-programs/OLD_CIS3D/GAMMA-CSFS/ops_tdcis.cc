/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_tdcis.cc                                                           *
 *                                                                              *
 * contains functions for  time depentend   operations  needed for tdcis        *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   * 
 *                                                    Pascal Krause      2007   *
 ********************************************************************************/ 
#include <complex>
#include <math.h>

using namespace std;
#define Complex complex<double>

//Functions
void trans_fore(double* mat, Complex* vec_i, Complex* vec_o, int cis_size); //small cis_mat
void trans_back(double* mat, Complex* vec_i, Complex* vec_o, int cis_size); //small cis_mat
void trans_fore(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size); //big cis_mat
void trans_back(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size); //big cis_mat
double opval(double* vals, Complex* vec, int cis_size);

//Extern Functions


/*******************************************************************************
 * Vector transformations for propagation                                      *
 *                                                                             *
 ******************************************************************************/
//small mats
void trans_fore(double* mat, Complex* vec_i, Complex* vec_o, int cis_size){
  int Cis_size = cis_size;
  for(int x = 0; x < Cis_size; x++){
    vec_o[x] = 0.;
    for(int y = 0; y < Cis_size;  y++){
      vec_o[x] += mat[x*Cis_size+y]*vec_i[y];
    }
  }
}

void trans_back(double* mat, Complex* vec_i, Complex* vec_o, int cis_size){
  int Cis_size = cis_size;
  for(int x = 0; x < Cis_size; x++){
    vec_o[x] = 0.;
    for(int y = 0; y < Cis_size;  y++){
      vec_o[x] += mat[y*Cis_size+x]*vec_i[y];
    }
  }
}

//big mats
void trans_fore(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size){
  long long int Cis_size = cis_size;
  for(long long int x = 0; x < Cis_size; x++){
    vec_o[x] = 0.;
    for(long long int y = 0; y < Cis_size;  y++){
      vec_o[x] += mat[x*Cis_size+y]*vec_i[y];
    }
  }
}

void trans_back(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size){
  long long int Cis_size = cis_size;
  for(long long int x = 0; x < Cis_size; x++){
    vec_o[x] = 0.;
    for(long long int y = 0; y < Cis_size;  y++){
      vec_o[x] += mat[y*Cis_size+x]*vec_i[y];
    }
  }
}


/*******************************************************************************
 * calc expetation value of cherm op in its eigenspace                         *
 *                                                                             *
 ******************************************************************************/


double opval(double* vals, Complex* vec, int cis_size){
  double val = 0.;
  for(int x = 0; x < cis_size; x++)
    val += real(vec[x]*conj(vec[x]))*vals[x];
  return(val);
}  

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      blas_m_vec                                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
extern "C" void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, 
                       double* X, int*  INCX, double* BETA, double* Y, int* INCY);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      TRANS_FORE                                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void trans_f(int np, double* U, Complex* vi, Complex* vo, double* tmpvec){
  static char TRANS = 'T';

  static int  M     = np;
  static int  N     = np;
  static int  LDA   = np;

  static int  INCX  =  1;
  static int  INCY  =  1;
  
  static double ALPHA = 1.;
  static double BETA  = 0.;
  
  for(int x = 0; x < np; x++) vo[x] = 0.;
  for(int x = 0; x < np; x++) tmpvec[x] = real(vi[x]);
  
  dgemv_(&TRANS, &M, &N, &ALPHA, U, &LDA, tmpvec, &INCX, &BETA, &(tmpvec[np]), &INCY);

  for(int x = 0; x < np; x++) vo[x] = tmpvec[np+x]; 
  for(int x = 0; x < np; x++) tmpvec[x] = imag(vi[x]);
  
  dgemv_(&TRANS, &M, &N, &ALPHA, U, &LDA, tmpvec, &INCX, &BETA, &(tmpvec[np]), &INCY);

  for(int x = 0; x < np; x++) vo[x] += Complex(0.,1.)*tmpvec[np+x]; 
  
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      TRANS_BACK                                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void trans_b(int np, double* U, Complex* vi, Complex* vo, double* tmpvec){
  static char TRANS = 'N';
  
  static int  M     = np;
  static int  N     = np;
  static int  LDA   = np;

  static int  INCX  =  1;
  static int  INCY  =  1;
  
  static double ALPHA = 1.;
  static double BETA  = 0.;
  
  for(int x = 0; x < np; x++) vo[x] = 0.;
  for(int x = 0; x < np; x++) tmpvec[x] = real(vi[x]);
  
  dgemv_(&TRANS, &M, &N, &ALPHA, U, &LDA, tmpvec, &INCX, &BETA, &(tmpvec[np]), &INCY);

  for(int x = 0; x < np; x++) vo[x] = tmpvec[np+x]; 
  for(int x = 0; x < np; x++) tmpvec[x] = imag(vi[x]);
  
  dgemv_(&TRANS, &M, &N, &ALPHA, U, &LDA, tmpvec, &INCX, &BETA, &(tmpvec[np]), &INCY);

  for(int x = 0; x < np; x++) vo[x] += Complex(0.,1.)*tmpvec[np+x]; 
}
