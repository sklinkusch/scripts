/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_tdcis.cc                                                           *
 *                                                                              *
 * contains functions for  time depentend   operations  needed for tdcis        *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
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

