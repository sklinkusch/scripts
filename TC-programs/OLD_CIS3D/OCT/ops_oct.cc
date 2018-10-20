/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_oct.cc                                                             *
 *                                                                              *
 * contains functions for oct propagations                                      * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <iostream>
#include <sstream>
#include <complex>

using namespace std;

#define Complex complex<double>

//Functions
void single_prop(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
		 double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
		 double* cisvecs,  double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
		 int cis_size, int nrots, double dt, double* field_h, char* pref, 
		 int mw_exp, int mw_pop, int mw_wav, ofstream *outf, int extra_ts, double* tmpvec);

void prop_fore(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
	       double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
	       double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
	       int cis_size, int ts, double dt, int nrots, double* field_h, double* tmpvec);

void prop_fore_mt(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
		  double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
		  double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
		  Complex* ChiE, Complex* ChiM_x, Complex* ChiM_y, Complex* ChiM_z,
		  double* MTcis_vals, double* MTmu_vals_x,  double* MTmu_vals_y, double* MTmu_vals_z,
		  double* MTmu_vecs_x,  double* MTmu_vecs_y, double* MTmu_vecs_z,
		  int cis_size, int ts, double dt, int nrots, double* field_h, double* field_b);

void apply_op(double* target_op, Complex* vec_in, Complex* vec_out, int cis_size);
void apply_op(double* target_op, Complex* vec_in, Complex* vec_out, long long int cis_size);
double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, int cis_size);
double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, long long int cis_size);
double calc_field_int(double* field, int nrots, double dt, double* shape);
double calc_field(Complex* wav_m1, Complex* wav_m2, double* mv, int cis_size, double penal);

//Extern Functions
extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, int cis_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, int cis_size);
extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size);
extern void   trans_f(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);
extern void   trans_b(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);
extern double opval(double* vals, Complex* vec, int cis_size);

void single_prop(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
		 double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
		 double* cisvecs,  double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
		 int cis_size, int nrots, double dt, double* field_h, char* pref, 
		 int mw_exp, int mw_pop, int mw_wav, ofstream *outf, int extra_ts, double* tmpvec){
  char dumc[256];
  long long int Cis_size = cis_size;
  static double*  zero_field = new double[3]; 
  
  for(int y = 0; y < 3; y++) zero_field[y] = 0.;

  *outf << "-------------------------------------------------------------------------------\n";
  *outf << "Single propagation function:\n";
  
  sprintf(dumc,"%s.pop",pref);
  ofstream popf(dumc);
  *outf << "Writing time dependent populations to " << dumc << "\n";
  
  sprintf(dumc,"%s.fad",pref);
  ofstream fadf(dumc);
  fadf << "#FORMAT:    $1     $2    $3   $4   $5   $6   $7  \n"
       << "#          time   E_x   E_y  E_z  Mu_x Mu_y Mu_z  \n";
  *outf << "Writing time dependent field and dipole to " << dumc << "\n";
  
  sprintf(dumc,"%s.wav",pref);
  ofstream wavf(dumc);
  *outf << "Writing time dependent PsiE to " << dumc << "\n";
  
  *outf << "Transfroming to dipole spaces\n";
//  if(cis_size < 10000){
//    trans_fore(mu_vecs_x, PsiE, PsiM_x, Cis_size);
//    trans_fore(mu_vecs_y, PsiE, PsiM_y, Cis_size);
//    trans_fore(mu_vecs_z, PsiE, PsiM_z, Cis_size); 
//  }else{
//    trans_fore(mu_vecs_x, PsiE, PsiM_x, cis_size);
//    trans_fore(mu_vecs_y, PsiE, PsiM_y, cis_size);
//    trans_fore(mu_vecs_z, PsiE, PsiM_z, cis_size); 
//  }
   trans_f(cis_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
   trans_f(cis_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
   trans_f(cis_size, mu_vecs_z, PsiE, PsiM_z, tmpvec);

  *outf << "Initial analysis\n";
  *outf << "E mx my mz: " << opval(cis_vals, PsiE,  cis_size) << " " <<  opval(mu_vals_x, PsiM_x,  cis_size)
       << " " <<  opval(mu_vals_y, PsiM_y,  cis_size)<< " " <<  opval(mu_vals_z, PsiM_z,  cis_size) << "\n";
  
  *outf << "-------------------------------------------------------------------------------\n";
  *outf << "Starting propagation\n";
  outf->flush();
  for(int t = 0; t < nrots+extra_ts; t++){
    if(t%mw_exp == 0){
      
      double E_el =  opval(cis_vals, PsiE,  cis_size);
      
      *outf <<  pref << " pTime: " << t*dt << " " << E_el << "\n";
      outf->flush();
      
      //Fields and dipole
//      if(cis_size < 10000){
//	trans_fore(mu_vecs_x, PsiE, PsiM_x, cis_size);
//	trans_fore(mu_vecs_y, PsiE, PsiM_y, cis_size);
//	trans_fore(mu_vecs_z, PsiE, PsiM_z, cis_size);
//      }else{
//	trans_fore(mu_vecs_x, PsiE, PsiM_x, Cis_size);
//	trans_fore(mu_vecs_y, PsiE, PsiM_y, Cis_size);
//	trans_fore(mu_vecs_z, PsiE, PsiM_z, Cis_size);
//      }
      trans_f(cis_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
      trans_f(cis_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
      trans_f(cis_size, mu_vecs_z, PsiE, PsiM_z, tmpvec);

      
      if(t < nrots)
	fadf << t*dt << " " << field_h[0*nrots+t] << " " << field_h[1*nrots+t] << " " << field_h[2*nrots+t] << " ";
      else
	fadf << t*dt << " 0. 0. 0. ";
      
      fadf   <<  opval(mu_vals_x, PsiM_x,  cis_size)
	   << " " <<  opval(mu_vals_y, PsiM_y,  cis_size)<< " " <<  opval(mu_vals_z, PsiM_z,  cis_size) << "\n";
      fadf.flush();
      
    }
    
    if(t%mw_pop == 0){
      for(int x = 0; x < cis_size; x++){
	popf << t*dt << " " << cis_vals[x]  << " " << norm(PsiE[x]) << "\n";
      }
      popf << "\n";
      popf.flush();
    }
    
    double curr_time = t*dt;
    if(t%mw_wav == 0){
      wavf.write((char* ) &curr_time, sizeof(double));
      wavf.write((char* ) PsiE, sizeof(Complex)*cis_size);
      wavf.flush();
    }
    if(t < nrots)
      prop_fore(PsiE,  PsiM_x,  PsiM_y,  PsiM_z,
		cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z,
		mu_vecs_x,   mu_vecs_y,  mu_vecs_z,
		cis_size, t, dt, nrots, field_h, tmpvec);
    else
      prop_fore(PsiE,  PsiM_x,  PsiM_y,  PsiM_z,
		cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z,
		mu_vecs_x,   mu_vecs_y,  mu_vecs_z,
		cis_size, 0, dt, 0, zero_field, tmpvec);
  }
}

//forward propagation
void prop_fore(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
	       double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
	       double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
	       int cis_size, int ts, double dt, int nrots, double* field_h, double* tmpvec){
  
  long long int Cis_size = cis_size;
  
  static double Efield[3];
  
  //Apply H
  for(int x = 0; x < cis_size; x++)
    PsiE[x] *= exp(-Complex(0.,1.)*cis_vals[x]*dt);
  
  //get_Efield 
  for(int x = 0; x < 3; x++) Efield[x] = field_h[x*nrots+ts];
    
  
  //Apply dipole ops
  //X
  if(Efield[0] != 0.){
//    if(cis_size < 10000)
//      trans_fore(mu_vecs_x, PsiE, PsiM_x, cis_size);
//    else
//      trans_fore(mu_vecs_x, PsiE, PsiM_x, Cis_size);
    trans_f(cis_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
    for(int x = 0; x < cis_size; x++)
      PsiM_x[x] *= exp(Complex(0.,1.)*mu_vals_x[x]*Efield[0]*dt);
//    if(cis_size < 10000)
//      trans_back(mu_vecs_x, PsiM_x, PsiE, cis_size);
//    else
//      trans_back(mu_vecs_x, PsiM_x, PsiE, Cis_size);
    trans_b(cis_size, mu_vecs_x, PsiM_x, PsiE, tmpvec);
  }
    
    
  //Y
  if(Efield[1] != 0.){
//    if(cis_size < 10000)
//      trans_fore(mu_vecs_y, PsiE, PsiM_y, cis_size);
//    else
//      trans_fore(mu_vecs_y, PsiE, PsiM_y, Cis_size);
    trans_f(cis_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
    for(int x = 0; x < cis_size; x++)
      PsiM_y[x] *= exp(Complex(0.,1.)*mu_vals_y[x]*Efield[1]*dt);
    trans_b(cis_size, mu_vecs_y, PsiM_y, PsiE, tmpvec);
//    if(cis_size < 10000)
//      trans_back(mu_vecs_y, PsiM_y, PsiE, cis_size);
//    else
//      trans_back(mu_vecs_y, PsiM_y, PsiE, Cis_size);
  }
  
  //Z
  if(Efield[2] != 0.){
//    if(cis_size < 10000)
//      trans_fore(mu_vecs_z, PsiE, PsiM_z, cis_size);
//    else
//      trans_fore(mu_vecs_z, PsiE, PsiM_z, Cis_size);
    trans_f(cis_size, mu_vecs_z, PsiE, PsiM_z, tmpvec);
    for(int x = 0; x < cis_size; x++)
      PsiM_z[x] *= exp(Complex(0.,1.)*mu_vals_z[x]*Efield[2]*dt);
    trans_b(cis_size, mu_vecs_z, PsiM_z, PsiE, tmpvec);
//    if(cis_size < 10000)
//      trans_back(mu_vecs_z, PsiM_z, PsiE, cis_size);
//    else
//      trans_back(mu_vecs_z, PsiM_z, PsiE, Cis_size);
  }
  
}


//Apply target operator (small space)
void apply_op(double* target_op, Complex* vec_in, Complex* vec_out, int cis_size){
  for(int x = 0; x < cis_size; x++){
    vec_out[x] = 0.;
    for(int y = 0; y < cis_size; y++){
      vec_out[x] += target_op[x*cis_size+y]*vec_in[y];
    }
  }
}

//Apply target operator (large space)
void apply_op(double* target_op, Complex* vec_in, Complex* vec_out, long long int cis_size){
  for(long long int x = 0; x < cis_size; x++){
    vec_out[x] = 0.;
    for(long long int y = 0; y < cis_size; y++){
      vec_out[x] += target_op[x*cis_size+y]*vec_in[y];
    }
  }
}

//Calc expectation value of target operator and Apply it (small space)
double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, int cis_size){
  apply_op(target_op, vec_in, vec_dum, cis_size);
  double opval = 0.;
  for(int x = 0; x < cis_size; x++)
    opval += real(conj(vec_in[x])*vec_dum[x]);
  return(opval);
}

//Calc expectation value of target operator and Apply it (large space)
double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, long long int cis_size){
  apply_op(target_op, vec_in, vec_dum, cis_size);
  double opval = 0.;
  for(long long int x = 0; x < cis_size; x++)
    opval += real(conj(vec_in[x])*vec_dum[x]);
  return(opval);
}

//Calculates the field energy 
double calc_field_int(double* field, int nrots, double dt, double* shape){
  double fint = 0.;
  for(int x = 0; x < nrots*3; x++)
    fint += pow(field[x],2);
  
  return(fint*dt);
}

//Calculates the optimal field
double calc_field(Complex* wav_m1, Complex* wav_m2, double* mv, int cis_size, double penal){
  Complex mufac = 0.;
  for(int x = 0; x < cis_size; x++){
    mufac += conj(wav_m1[x])*mv[x]*wav_m2[x];
  }
  return(-1./penal*imag(mufac));
}


/*******************************************************************************
//MULTI-THREADED  PROPAGATION
*******************************************************************************/



typedef struct {
  double re, im;
} mycom;

extern "C" void  mk_thread_trans(int dir, double* mat1, double* mat2, mycom* vec_in1, mycom* vec_in2,
				 mycom* vec_out1, mycom* vec_out2, int* nrop1, int* nrop2);

extern "C" void lmk_thread_trans(int dir, double* mat1, double* mat2, mycom* vec_in1, mycom* vec_in2,
				 mycom* vec_out1, mycom* vec_out2, long long int* nrop1, long long int* nrop2);
// DIR=1 is forward !!

void prop_fore_mt(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
		  double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
		  double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
		  Complex* ChiE, Complex* ChiM_x, Complex* ChiM_y, Complex* ChiM_z,
		  double* MTcis_vals, double* MTmu_vals_x,  double* MTmu_vals_y, double* MTmu_vals_z,
		  double* MTmu_vecs_x,  double* MTmu_vecs_y, double* MTmu_vecs_z,
		  int cis_size, int ts, double dt, int nrots, double* field_h, double* field_b){
  
  long long int Cis_size = cis_size;
  long long int Cs2      = Cis_size;
  int           cs2      = cis_size;
  

  static double Efield_h[3];
  static double Efield_b[3];
  
  //Apply H
  for(int x = 0; x < cis_size; x++)
    PsiE[x] *= exp(-Complex(0.,1.)*cis_vals[x]*dt);


  //Apply H
  for(int x = 0; x < cis_size; x++)
    ChiE[x] *= exp(-Complex(0.,1.)*cis_vals[x]*dt);
  
  //get_Efield 
  for(int x = 0; x < 3; x++) Efield_h[x] = field_h[x*nrots+ts];
  for(int x = 0; x < 3; x++) Efield_b[x] = field_b[x*nrots+ts];
    
  
  //Apply dipole ops
  //X
  if(Efield_h[0] != 0.){
    if(cis_size < 10000)
      mk_thread_trans(1, mu_vecs_x, MTmu_vecs_x, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_x, (mycom* )ChiM_x,  &cis_size , &cs2);
    else
      lmk_thread_trans(1, mu_vecs_x, MTmu_vecs_x, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_x, (mycom* )ChiM_x,  &Cis_size , &Cs2);
    
    for(int x = 0; x < cis_size; x++)  PsiM_x[x] *= exp(Complex(0.,1.)*mu_vals_x[x]*Efield_h[0]*dt);
    for(int x = 0; x < cis_size; x++)  ChiM_x[x] *= exp(Complex(0.,1.)*mu_vals_x[x]*Efield_b[0]*dt);

    if(cis_size < 10000)
      mk_thread_trans(-1, mu_vecs_x, MTmu_vecs_x, (mycom* ) PsiM_x, (mycom* ) ChiM_x, (mycom* )PsiE, (mycom* )ChiE,  &cis_size , &cs2);
    else
      lmk_thread_trans(-1, mu_vecs_x, MTmu_vecs_x, (mycom* ) PsiM_x, (mycom* ) ChiM_x, (mycom* )PsiE, (mycom* )ChiE,  &Cis_size , &Cs2);
  }
    
    
  //Y
  if(Efield_h[1] != 0.){
    if(cis_size < 10000)
      mk_thread_trans(1, mu_vecs_y, MTmu_vecs_y, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_y, (mycom* )ChiM_y,  &cis_size , &cs2);
    else
      lmk_thread_trans(1, mu_vecs_y, MTmu_vecs_y, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_y, (mycom* )ChiM_y, &Cis_size , &Cs2);
    
    for(int x = 0; x < cis_size; x++)  PsiM_y[x] *= exp(Complex(0.,1.)*mu_vals_y[x]*Efield_h[1]*dt);
    for(int x = 0; x < cis_size; x++)  ChiM_y[x] *= exp(Complex(0.,1.)*mu_vals_y[x]*Efield_b[1]*dt);

    if(cis_size < 10000)
      mk_thread_trans(-1, mu_vecs_y, MTmu_vecs_y, (mycom* ) PsiM_y, (mycom* ) ChiM_y, (mycom* )PsiE, (mycom* )ChiE,  &cis_size , &cs2);
    else
      lmk_thread_trans(-1, mu_vecs_y, MTmu_vecs_y, (mycom* ) PsiM_y, (mycom* ) ChiM_y, (mycom* )PsiE, (mycom* )ChiE, &Cis_size , &Cs2);
  }
  
  //Z
  if(Efield_h[2] != 0.){
    if(cis_size < 10000)
      mk_thread_trans(1, mu_vecs_z, MTmu_vecs_z, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_z, (mycom* )ChiM_z,  &cis_size , &cs2);
    else
      lmk_thread_trans(1, mu_vecs_z, MTmu_vecs_z, (mycom* ) PsiE, (mycom* ) ChiE, (mycom* )PsiM_z, (mycom* )ChiM_z,  &Cis_size , &Cs2);
    
    for(int x = 0; x < cis_size; x++)  PsiM_z[x] *= exp(Complex(0.,1.)*mu_vals_z[x]*Efield_h[1]*dt);
    for(int x = 0; x < cis_size; x++)  ChiM_z[x] *= exp(Complex(0.,1.)*mu_vals_z[x]*Efield_b[1]*dt);

    if(cis_size < 10000)
      mk_thread_trans(-1, mu_vecs_z, MTmu_vecs_z, (mycom* ) PsiM_z, (mycom* ) ChiM_z, (mycom* )PsiE, (mycom* )ChiE,  &cis_size , &cs2);
    else
      lmk_thread_trans(-1, mu_vecs_z, MTmu_vecs_z, (mycom* ) PsiM_z, (mycom* ) ChiM_z, (mycom* )PsiE, (mycom* )ChiE,  &Cis_size , &Cs2);
  }
  
}
