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
#include <fstream>

using namespace std;

#define Complex complex<double>

//Functions
void diag_mat(int nroao, double* mat, double* vals, double* vecs);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mt_mt(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_mt_mult_a(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mt_mt_mult_b(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mt_mt_t(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_mt_t_mult_a(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mt_t_mt(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_t_mt_mult_b(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
void mt_mt(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_mt_mult_a(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mt_mt_mult_b(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mt_mt_t(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_mt_t_mult_a(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
void mt_t_mt(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf);
void mt_t_mt_mult_b(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf);
void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void trans_mt(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf);
void trans_mt_mult(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* matx, double* maty, double* matz, double* trans, double* tmpmatx, double* tmpmaty,
	double* tmpmatz, int dir, ofstream* outf);
void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir); //for large matrices
void trans_mt(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf);
void trans_mt_mult(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* matx, double* maty, double* matz, double* trans, 
	double* tmpmatx, double* tmpmaty, double* tmpmatz, int dir, ofstream* outf);
//Extern Functions
extern "C" void  dsyev_(char* JOBZ, char*  UPLO,int* N, double* A, int* LDA, 
                        double* W, double* WORK, int* LWORK, int*  INFO );

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
    double dumc = 0.;
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_mat for non-quadratic matrices (na x nc  times nc x nb)
void mt_mt(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt\n";
	exit(201);
    }else{
	ne = nc_m1;
    }
    double dumc = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction (+:dumc)
	    for(int z = 0; z < ne; z++) dumc += mat_i1[x*ne+z]*mat_i2[z*nc_m2+y];
	    mat_f[x*nc_m2+y] = dumc;
	}
    }
}

void mt_mt_mult_a(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt\n";
	exit(201);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1a[x*ne+z]*mat_i2[z*nc_m2+y];
	       dumc_b += mat_i1b[x*ne+z]*mat_i2[z*nc_m2+y];
	       dumc_c += mat_i1c[x*ne+z]*mat_i2[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}

void mt_mt_mult_b(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt\n";
	exit(201);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1[x*ne+z]*mat_i2a[z*nc_m2+y];
	       dumc_b += mat_i1[x*ne+z]*mat_i2b[z*nc_m2+y];
	       dumc_c += mat_i1[x*ne+z]*mat_i2c[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}

void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f){
    double dumc = 0.;
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[y*np+z];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_mat_t for non-quadratic matrices (na x nb times (nc x nb)^T)
void mt_mt_t(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    int ne = 0;
    if(nc_m1 != nc_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt_t\n";
	exit(202);
    }else{
	ne = nc_m1;
    }
    double dumc = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nr_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction (+:dumc)
	    for(int z = 0; z < ne; z++) dumc += mat_i1[x*ne+z]*mat_i2[y*ne+z];
	    mat_f[x*nr_m2+y] = dumc;
	}
    }
}

void mt_mt_t_mult_a(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    int ne = 0;
    if(nc_m1 != nc_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt_t\n";
	exit(202);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nr_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
		dumc_a += mat_i1a[x*ne+z]*mat_i2[y*ne+z];
		dumc_b += mat_i1b[x*ne+z]*mat_i2[y*ne+z];
		dumc_c += mat_i1c[x*ne+z]*mat_i2[y*ne+z];
	    }
	    mat_fa[x*nr_m2+y] = dumc_a;
	    mat_fb[x*nr_m2+y] = dumc_b;
	    mat_fc[x*nr_m2+y] = dumc_c;
	}
    }
}


void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
    double dumc = 0.;
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[z*np+x]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_t_mat for non-quadratic matrices ((na x nc)^T times na x nb)
void mt_t_mt(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    int ne = 0;
    if(nr_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_t_mt\n";
	exit(203);
    }else{
	ne = nr_m1;
    }
    double dumc = 0.;
    for(int x = 0; x < nc_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction( +:dumc)
	    for(int z = 0; z < ne; z++) dumc += mat_i1[z*nc_m1+x]*mat_i2[z*nc_m2+y];
	    mat_f[x*nc_m2+y] = dumc;
	}
    }
}

void mt_t_mt_mult_b(int nr_m1, int nc_m1, int nr_m2, int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    int ne = 0;
    if(nr_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_t_mt\n";
	exit(203);
    }else{
	ne = nr_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nc_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction( +:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1[z*nc_m1+x]*mat_i2a[z*nc_m2+y];
	       dumc_b += mat_i1[z*nc_m1+x]*mat_i2b[z*nc_m2+y];
	       dumc_c += mat_i1[z*nc_m1+x]*mat_i2c[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}

//For large matrices
void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
    double dumc = 0.;
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(long long int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

void mt_mt(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    long long int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt (large)\n";
	exit(204);
    }else{
	ne = nc_m1;
    }
    double dumc = 0.;
    for(long long int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(long long int y = 0; y < nc_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction (+:dumc)
	    for(long long int z = 0; z < ne; z++) dumc += mat_i1[x*ne+z]*mat_i2[z*nc_m2+y];
	    mat_f[x*nc_m2+y] = dumc;
	}
    }
}

void mt_mt_mult_a(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    long long int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt\n";
	exit(201);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1a[x*ne+z]*mat_i2[z*nc_m2+y];
	       dumc_b += mat_i1b[x*ne+z]*mat_i2[z*nc_m2+y];
	       dumc_c += mat_i1c[x*ne+z]*mat_i2[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}

void mt_mt_mult_b(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    long long int ne = 0;
    if(nc_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt\n";
	exit(201);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1[x*ne+z]*mat_i2a[z*nc_m2+y];
	       dumc_b += mat_i1[x*ne+z]*mat_i2b[z*nc_m2+y];
	       dumc_c += mat_i1[x*ne+z]*mat_i2c[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}


void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
    double dumc = 0.;
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(long long int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[y*np+z];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_mat_t for non-quadratic matrices (na x nb times (nc x nb)^T)
void mt_mt_t(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    long long int ne = 0;
    if(nc_m1 != nc_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt_t (large)\n";
	exit(205);
    }else{
	ne = nc_m1;
    }
    double dumc = 0.;
    for(long long int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(long long int y = 0; y < nr_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction (+:dumc)
	    for(long long int z = 0; z < ne; z++) dumc += mat_i1[x*ne+z]*mat_i2[y*ne+z];
	    mat_f[x*nr_m2+y] = dumc;
	}
    }
}

void mt_mt_t_mult_a(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1a, double* mat_i1b, double* mat_i1c, double* mat_i2, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    long long int ne = 0;
    if(nc_m1 != nc_m2){
	cerr << "Matrices cannot be multiplied in routine mt_mt_t\n";
	exit(202);
    }else{
	ne = nc_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nr_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nr_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction (+:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
		dumc_a += mat_i1a[x*ne+z]*mat_i2[y*ne+z];
		dumc_b += mat_i1b[x*ne+z]*mat_i2[y*ne+z];
		dumc_c += mat_i1c[x*ne+z]*mat_i2[y*ne+z];
	    }
	    mat_fa[x*nr_m2+y] = dumc_a;
	    mat_fb[x*nr_m2+y] = dumc_b;
	    mat_fc[x*nr_m2+y] = dumc_c;
	}
    }
}


void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f){
    double dumc = 0.;
  for(long long int x = 0; x < np; x++){
    for(long long int y = 0; y < np; y++){
	dumc = 0.;
//      mat_f[x*np+y] = 0.;
#pragma omp parallel for reduction (+:dumc)
      for(long long int z = 0; z < np; z++){
	  dumc += mat_i1[z*np+x]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_t_mat for non-quadratic matrices ((na x nc)^T times na x nb)
void mt_t_mt(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2, double* mat_f, ofstream* outf){
    long long int ne = 0;
    if(nr_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_t_mt (large)\n";
	exit(206);
    }else{
	ne = nr_m1;
    }
    double dumc = 0.;
    for(long long int x = 0; x < nc_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(long long int y = 0; y < nc_m2; y++){
	    dumc = 0.;
#pragma omp parallel for reduction( +:dumc)
	    for(long long int z = 0; z < ne; z++) dumc += mat_i1[z*nc_m1+x]*mat_i2[z*nc_m2+y];
	    mat_f[x*nc_m2+y] = dumc;
	}
    }
}

void mt_t_mt_mult_b(long long int nr_m1, long long int nc_m1, long long int nr_m2, long long int nc_m2, double* mat_i1, double* mat_i2a, double* mat_i2b, double* mat_i2c, 
	double* mat_fa, double* mat_fb, double* mat_fc, ofstream* outf){
    long long int ne = 0;
    if(nr_m1 != nr_m2){
	cerr << "Matrices cannot be multiplied in routine mt_t_mt\n";
	exit(203);
    }else{
	ne = nr_m1;
    }
    double dumc_a = 0.;
    double dumc_b = 0.;
    double dumc_c = 0.;
    for(int x = 0; x < nc_m1; x++){
	*outf << x << "\t";
	if(x%10 == 9) *outf << "\n";
        outf->flush();
	for(int y = 0; y < nc_m2; y++){
	    dumc_a = 0.;
	    dumc_b = 0.;
	    dumc_c = 0.;
#pragma omp parallel for reduction( +:dumc_a, dumc_b, dumc_c)
	    for(int z = 0; z < ne; z++){
	       dumc_a += mat_i1[z*nc_m1+x]*mat_i2a[z*nc_m2+y];
	       dumc_b += mat_i1[z*nc_m1+x]*mat_i2b[z*nc_m2+y];
	       dumc_c += mat_i1[z*nc_m1+x]*mat_i2c[z*nc_m2+y];
	    }
	    mat_fa[x*nc_m2+y] = dumc_a;
	    mat_fb[x*nc_m2+y] = dumc_b;
	    mat_fc[x*nc_m2+y] = dumc_c;
	}
    }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     SIMPLE MATRIX ROUTINES   (COMPLEX)                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void mat_mat(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = complex<double>(0.,0.);
//      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_mat for non-quadratic matrices (na x nc  times nc x nb)
void mt_mt(int na, int nb, int nc, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
    for(int x = 0; x < na; x++){
	for(int y = 0; y < nb; y++){
	    dumc = complex<double>(0.,0.);
	    for(int z = 0; z < nc; z++) dumc += mat_i1[x*nc+z]*mat_i2[z*nb+y];
	    mat_f[x*nb+y] = dumc;
	}
    }
}


void mat_mat_T(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = complex<double>(0.,0.);
//      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[x*np+z]*mat_i2[y*np+z];
//        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_mat_t for non-quadratic matrices (na x nb times (nc x nb)^T)
void mt_mt_t(int na, int nb, int nc, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
    for(int x = 0; x < na; x++){
	for(int y = 0; y < nc; y++){
	    dumc = complex<double>(0.,0.);
	    for(int z = 0; z < nb; z++) dumc += mat_i1[x*nb+z]*mat_i2[y*nb+z];
	    mat_f[x*nc+y] = dumc;
	}
    }
}


void mat_T_mat(int np, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
	dumc = complex<double>(0.,0.);
//      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
	  dumc += mat_i1[z*np+x]*mat_i2[z*np+y];
//        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
      mat_f[x*np+y] = dumc;
    }
  }
}

// mat_t_mat for non-quadratic matrices ((na x nc)^T times na x nb)
void mt_t_mt(int na, int nb, int nc, Complex* mat_i1, Complex* mat_i2, Complex* mat_f){
    Complex dumc = complex<double>(0.,0.);
    for(int x = 0; x < nc; x++){
	for(int y = 0; y < nb; y++){
	    dumc = complex<double>(0.,0.);
	    for(int z = 0; z < na; z++) dumc += mat_i1[z*nc+x]*mat_i2[z*nb+y];
	    mat_f[x*nb+y] = dumc;
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

void trans_mt(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf){

    if(dir == 0) // CALC mat*trans
	mt_mt(nr_mat, nc_mat, nr_trans, nc_trans, mat, trans, tmpmat, outf);
    else         // CALC trans*mat
	mt_mt(nr_trans, nc_trans, nr_mat, nc_mat, trans, mat, tmpmat, outf);
    *outf << "First multiplication finished, starting second multiplication\n";
    outf->flush();

    if(dir == 0) // CALC trans^T*(mat*trans)
	mt_t_mt(nr_trans, nc_trans, nr_mat, nc_trans, trans, tmpmat, mat, outf);
    else         // CALC (trans*mat)*trans^T
	mt_mt_t(nr_mat, nc_trans, nr_trans, nc_trans, tmpmat, trans, mat, outf);
    *outf << "Second multiplication finished\n";
    outf->flush();
}

void trans_mt_mult(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* matx, double* maty, double* matz, double* trans, double* tmpmatx,
	double* tmpmaty, double* tmpmatz, int dir, ofstream* outf){

    if(dir == 0) // CALC mat*trans
	mt_mt_mult_a(nr_mat, nc_mat, nr_trans, nc_trans, matx, maty, matz, trans, tmpmatx, tmpmaty, tmpmatz, outf);
    else         // CALC trans*mat
	mt_mt_mult_b(nr_trans, nc_trans, nr_mat, nc_mat, trans, matx, maty, matz, tmpmatx, tmpmaty, tmpmatz, outf);
    *outf << "First multiplication finished, starting second multiplication\n";
    outf->flush();

    if(dir == 0) // CALC trans^T*(mat*trans)
	mt_t_mt_mult_b(nr_trans, nc_trans, nr_mat, nc_trans, trans, tmpmatx, tmpmaty, tmpmatz, matx, maty, matz, outf);
    else         // CALC (trans*mat)*trans^T
	mt_mt_t_mult_a(nr_mat, nc_trans, nr_trans, nc_trans, tmpmatx, tmpmaty, tmpmatz, trans, matx, maty, matz, outf);
    *outf << "Second multiplication finished\n";
    outf->flush();
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

void trans_mt(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf){

    if(dir == 0) // CALC mat*trans
	mt_mt(nr_mat, nc_mat, nr_trans, nc_trans, mat, trans, tmpmat, outf);
    else         // CALC trans*mat
	mt_mt(nr_trans, nc_trans, nr_mat, nc_mat, trans, mat, tmpmat, outf);
    *outf << "First multiplication finished, starting second multiplication\n";
    outf->flush();

    if(dir == 0) // CALC trans^T*(mat*trans)
	mt_t_mt(nr_trans, nc_trans, nr_mat, nc_trans, trans, tmpmat, mat, outf);
    else         // CALC (trans*mat)*trans^T
	mt_mt_t(nr_mat, nc_trans, nr_trans, nc_trans, tmpmat, trans, mat, outf);
    *outf << "Second multiplication finished\n";
    outf->flush();
}

void trans_mt_mult(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* matx, double* maty, double* matz, double* trans, 
	double* tmpmatx, double* tmpmaty, double* tmpmatz, int dir, ofstream* outf){

    if(dir == 0) // CALC mat*trans
	mt_mt_mult_a(nr_mat, nc_mat, nr_trans, nc_trans, matx, maty, matz, trans, tmpmatx, tmpmaty, tmpmatz, outf);
    else         // CALC trans*mat
	mt_mt_mult_b(nr_trans, nc_trans, nr_mat, nc_mat, trans, matx, maty, matz, tmpmatx, tmpmaty, tmpmatz, outf);
    *outf << "First multiplication finished, starting second multiplication\n";
    outf->flush();

    if(dir == 0) // CALC trans^T*(mat*trans)
	mt_t_mt_mult_b(nr_trans, nc_trans, nr_mat, nc_trans, trans, tmpmatx, tmpmaty, tmpmatz, matx, maty, matz, outf);
    else         // CALC (trans*mat)*trans^T
	mt_mt_t_mult_a(nr_mat, nc_trans, nr_trans, nc_trans, tmpmatx, tmpmaty, tmpmatz, trans, matx, maty, matz, outf);
    *outf << "Second multiplication finished\n";
    outf->flush();
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

void trans_mt(int na, int nb, int nc, Complex* mat, Complex* trans, Complex* tmpmat, int dir){

    if(dir == 0) // CALC mat*trans
	mt_mt(na, nb, nc, mat, trans, tmpmat);
    else         // CALC trans*mat
	mt_mt(na, nb, nc, trans, mat, tmpmat);

    if(dir == 0) // CALC trans^T*(mat*trans)
	mt_t_mt(na, nb, nc, trans, tmpmat, mat);
    else         // CALC (trans*mat)*trans^T
	mt_mt_t(na, nb, nc, tmpmat, trans, mat);
}

