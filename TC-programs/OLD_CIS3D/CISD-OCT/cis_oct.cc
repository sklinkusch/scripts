 /********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cis_oct.cc                                                             *
 *                                                                              *
 * oct propagations                                                   *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
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
#include <stdlib.h>

using namespace std;
#define Complex complex<double>


//Functions


//Extern Functions
extern int    rem_com(char* filename, char* streamstring, int string_length);
extern void   status(ofstream* outf);
extern void single_prop(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
			double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
			double* cisvecs,  double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
			int cis_size, int nrots, double dt, double* field_h, char* pref, 
			int mw_exp, int mw_pop, int mw_wav, ofstream *outf, int extra_ts, double* tmpvec);

extern double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, long long int cis_size);
extern double calc_op(double* target_op, Complex* vec_in, Complex* vec_dum, int cis_size);
extern double calc_field_int(double* field, int nrots, double dt, double* shape);
extern double calc_field(Complex* wav_m1, Complex* wav_m2, double* mv, int cis_size, double penal);

extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, int cis_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, int cis_size);
extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, long long int cis_size);
extern void   trans_f(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);
extern void   trans_b(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);

extern void prop_fore(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
		      double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
		      double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
		      int cis_size, int ts, double dt, int nrots, double* field_h, double* tmpvec);

extern void prop_fore_mt(Complex* PsiE, Complex* PsiM_x, Complex* PsiM_y, Complex* PsiM_z,
			 double* cis_vals, double* mu_vals_x,  double* mu_vals_y, double* mu_vals_z,
			 double* mu_vecs_x,  double* mu_vecs_y, double* mu_vecs_z,
			 Complex* ChiE, Complex* ChiM_x, Complex* ChiM_y, Complex* ChiM_z,
			 double* MTcis_vals, double* MTmu_vals_x,  double* MTmu_vals_y, double* MTmu_vals_z,
			 double* MTmu_vecs_x,  double* MTmu_vecs_y, double* MTmu_vecs_z,
			 int cis_size, int ts, double dt, int nrots, double* field_h, double* field_b);

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  double ddum;
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "CIS_OCT [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";
 
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

 //SYSTEM
  char    bcsfile[256];   //binary system file
  char    iwavfile[256];  //inital wave function file if nrots < 0;
  char    opfile[256];    //operator or target wave function file
  char    effile[256];
  int     use_d    = 0;   //if 1 use the 2nd PT corrected energies for propagations
  int     use_iwav = 0;   //if 1 inital wavefunction is read in
  int     op_mode  = 0;   //if 0 use target wave function, if 1 (real) target operator

  //Propagation
  int     nrots;          //nr of time steps( if < 0 initial wave function will be read in, otherwise HF-Ground state)
  double  dt;             //time step length
  int     nroit;          //nr of oct iterations
  int     mw_exp;         //analy intervall expectation values
  int     mw_pop;         //analy intervall for populations
  int     mw_wav;         //write out intervall wave function
  
  //OCT
  int     pfac[3];        //polarization vector
  double  penal;          //penalty factor

  //Multi theread
  int     use_mt   = 0;   //use 2 threads to propagate the wave function 0=false 1=true
  int     extra_ts = 0;   //extra propagation steps with zero field after the control time
  
  ist >> bcsfile >> use_d;
  outf << "\nSystem data \nsys-file: " << bcsfile << "\n";
  if(use_d == 1) outf << "Using doubles PT corrections!\n\n";
  
  ist >> use_iwav;
  if(use_iwav == 1){
    ist >> iwavfile; 
    outf << "Reading initial wavefunction from " << iwavfile << "\n";
  }else
    outf << "Using ground state as inital wave function\n";
  
  ist >> opfile >> op_mode;
  if(op_mode == 0) 
    outf << "\nTarget wave function is read from " << opfile << "\n";
  else
    outf << "\nTarget operator is read from " << opfile << "\n";
  
  ist >> effile;
  outf << "Initial field is read from " << effile << "\n";
  
  ist >> nroit;
  outf << "Number of oct iterations is: " << nroit << "\n";

  ist >> penal;
  outf << "Penalty factor is " << penal << "\n";

  ist >> mw_exp >> mw_pop >> mw_wav;
  outf << "\nWrite out for initial and final forward propagation:\n"
       << "\nExpectation values  written out every " << mw_exp << " time steps" 
       << "\nPopulations         written out every " << mw_pop << " time steps" 
       << "\nMOs                 written out every " << mw_wav << " time steps\n\n";
  
  ist >> use_mt;
  if(use_mt == 1)
    outf << "Using 2 threads to propagate the wave function\n";

  ist >> extra_ts;
  outf << "Using " << extra_ts << " time steps after the controll time\n";
  
  outf << "*******************************************************************************\n";
  //READ SYSTEM DATA
  //CIS READ IN
  ifstream datf(bcsfile);

  int nroao, nroe, llim, ulim;

  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));

  outf << "CIS data read: \n";
  int nrof = ulim - llim +1;
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;

   // Loops to determine cisd_size
  int sum2 = 0;
  for (int occ = 1; occ <= omo-1; occ++){
    sum2 += occ; 
  }

  int sum1 = 0;
  for (int virt = 1; virt <= umo-1; virt++){
    sum1 += virt; 
  }

  int cis_size = 1 + 2*omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;
  long long int Cis_size = cis_size;

  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  outf << "nr of CSF               : " << cis_size << "\n";


  //VARIABLES
  outf << "Allocating " << (long long int) cis_size * (long long int) cis_size * 5 *8  + 4* cis_size *8 << " bytes of memory for CIS operators\n";
  outf.flush();

  double* cis_vals  = new double[cis_size];        //eigenvlaues of CI-matrix
  double* mu_vals_x = new double[cis_size];        //eigenvalues of dipole operator in x
  double* mu_vals_y = new double[cis_size];        //eigenvalues of dipole operator in x
  double* mu_vals_z = new double[cis_size];        //eigenvalues of dipole operator in x
  double* corr_vals = new double[cis_size];        //PT corrections to CIS ex ens
  double* tmpvec    = new double[2*cis_size];

  double* cisvecs   = new double[(long long int) cis_size * (long long int) cis_size];        //eigenvectros of CI-matrix
  double* mu_vecs_x = new double[(long long int) cis_size * (long long int) cis_size];        //eigenvectors of dipole operator in x
  double* mu_vecs_y = new double[(long long int) cis_size * (long long int) cis_size];        //eigenvectors of dipole operator in y
  double* mu_vecs_z = new double[(long long int) cis_size * (long long int) cis_size];        //eigenvectors of dipole operator in z

  double* target_op = new double[(long long int) cis_size * (long long int) cis_size];        //target operator is CIS space

  Complex* Psi0;                              //Inital wave function
  Complex* PsiE;                              //Wave funcion in energy space
  Complex* PsiM_x;                            //Wave function in eigen space of mu_x
  Complex* PsiM_y;                            //Wave function in eigen space of mu_y
  Complex* PsiM_z;                            //Wave function in eigen space of mu_z
  Complex* ChiE;                              //Wave funcion in energy space
  Complex* ChiM_x;                            //Wave function in eigen space of mu_x
  Complex* ChiM_y;                            //Wave function in eigen space of mu_y
  Complex* ChiM_z;                            //Wave function in eigen space of mu_z

  //MULTI THREAD VARIABLES
  double* MTcis_vals = NULL;
  double* MTmu_vals_x= NULL; 
  double* MTmu_vals_y= NULL; 
  double* MTmu_vals_z= NULL; 
  double* MTmu_vecs_x= NULL; 
  double* MTmu_vecs_y= NULL;
  double* MTmu_vecs_z= NULL; 
  
  if(use_mt == 1){
    outf << "Allocating extra Memory for multi thread propagation!\n";
    MTcis_vals    = new double[cis_size]; 
    MTmu_vals_x   = new double[cis_size];
    MTmu_vals_y   = new double[cis_size];
    MTmu_vals_z   = new double[cis_size];
    MTmu_vecs_x   = new double[(long long int) cis_size * (long long int) cis_size];
    MTmu_vecs_y   = new double[(long long int) cis_size * (long long int) cis_size];
    MTmu_vecs_z   = new double[(long long int) cis_size * (long long int) cis_size];
  }
  
  int complmem = 9*cis_size; int inc = 0;
  Complex* dummc = new Complex[complmem];
  
  Psi0 = &(dummc[inc]); inc += cis_size; 
  PsiE = &(dummc[inc]); inc += cis_size; 
  PsiM_x = &(dummc[inc]); inc += cis_size; PsiM_y = &(dummc[inc]); inc += cis_size; PsiM_z = &(dummc[inc]); inc += cis_size; 
  ChiE = &(dummc[inc]); inc += cis_size; 
  ChiM_x = &(dummc[inc]); inc += cis_size; ChiM_y = &(dummc[inc]); inc += cis_size; ChiM_z = &(dummc[inc]); inc += cis_size; 
    
  outf << "Reading data:\n";
   
  datf.read((char *) cis_vals, cis_size * sizeof(double));
  datf.read((char *) cisvecs, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "First uncorrected excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
  outf << "Highest uncorrected excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";

  datf.read((char *) mu_vals_x, cis_size * sizeof(double));
  datf.read((char *) mu_vecs_x, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "Dipole eigenvalues in x range from " << mu_vals_x[0] << " to " << mu_vals_x[cis_size-1] << "\n";

  datf.read((char *) mu_vals_y, cis_size * sizeof(double));
  datf.read((char *) mu_vecs_y, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "Dipole eigenvalues in y range from " << mu_vals_y[0] << " to " << mu_vals_y[cis_size-1] << "\n";

  datf.read((char *) mu_vals_z, cis_size * sizeof(double));
  datf.read((char *) mu_vecs_z, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "Dipole eigenvalues in z range from " << mu_vals_z[0] << " to " << mu_vals_z[cis_size-1] << "\n";
  if(use_mt == 1){
    //COPY orgignal data to MT data
    memcpy((void *) MTcis_vals, cis_vals, cis_size*sizeof(double));
    memcpy((void *) MTmu_vals_x, mu_vals_x, cis_size*sizeof(double));
    memcpy((void *) MTmu_vals_y, mu_vals_y, cis_size*sizeof(double));
    memcpy((void *) MTmu_vals_z, mu_vals_z, cis_size*sizeof(double));
    
    memcpy((void *) MTmu_vecs_x, mu_vecs_x, Cis_size*Cis_size*sizeof(double));
    memcpy((void *) MTmu_vecs_y, mu_vecs_y, Cis_size*Cis_size*sizeof(double));
    memcpy((void *) MTmu_vecs_z, mu_vecs_z, Cis_size*Cis_size*sizeof(double));
  }
  outf << "-------------------------------------------------------------------------------\n";


  outf.flush();

  if(use_d == 1){
    outf << "Reading PT-correction\n";
    datf.read((char *) corr_vals, cis_size * sizeof(double));
    for(int x = 0; x < cis_size; x++)
      cis_vals[x] +=  corr_vals[x];
    outf << "First CORRECTED excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
    outf << "Highest CORRECTED excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";
    outf << "-------------------------------------------------------------------------------\n";
  }

  datf.close();
  
  outf << "*******************************************************************************\n";
  outf << "Setting up initial wave function\n";
  for(int x = 0; x < cis_size; x++) Psi0[x] = 0.; Psi0[0] = 1.;
  if(use_iwav == 1){
    ifstream iw(iwavfile);
    iw.read((char *) &ddum, sizeof(double));
    iw.read((char *) Psi0, cis_size*sizeof(Complex));
    outf << "Read in from " << iwavfile << "\n";
  }else{
    outf << "Starting from the ground state.\n";
  }
  outf << "*******************************************************************************\n";
  outf << "Setting up target operator\n";
  if(op_mode == 0){
    outf << "Reading target wave from " << opfile << "\n";
    ifstream opf(opfile);
    opf.read((char *) &ddum, sizeof(double));
    opf.read((char *) ChiE, cis_size*sizeof(Complex));
    double max_im = 0.;
    for(int x = 0; x < cis_size; x++){if(fabs(imag(ChiE[x]))> max_im) max_im = imag(ChiE[x]);}
    outf << "Largest imaginary part is :" << max_im << "\n";
    for(long long int x = 0; x < Cis_size; x++){
      for(long long int y = 0; y < Cis_size; y++)
	target_op[x*Cis_size+y] = real(ChiE[x])*real(ChiE[y]);
    }
  }else{
    ifstream opf(opfile);
    opf.read((char *) target_op, Cis_size*Cis_size*sizeof(double));
    outf << "Target operator read from " << opfile << "\n";
  }
  outf << "*******************************************************************************\n";
  outf << "Reading inital field\n";
  ifstream ief(effile);
  ief.read((char *) &nrots, sizeof(int));
  ief.read((char *) &dt, sizeof(double));
  ief.read((char *) pfac, 3*sizeof(int));
  outf << "Nr of time steps: " << nrots << "\n";
  outf << "dt:               " << dt << "\n";
  outf << "Polarization:     " << pfac[0] << " " << pfac[1] << " " << pfac[2] << "\n";
  double* shape       = new double[nrots];
  double* field_h     = new double[3*nrots];
  double* field_b     = new double[3*nrots];
  ief.read((char *) field_h, 3*nrots*sizeof(double));
  ief.read((char *) shape,     nrots*sizeof(double));
  outf << "*******************************************************************************\n";
  outf << "Inital forward propagation\n";
  outf.flush();
  sprintf(dumc,"%s_I",argv[2]);
  for(int x = 0; x< cis_size; x++) PsiE[x] = Psi0[x];
  if(nroit < 1)
    single_prop( PsiE,  PsiM_x,  PsiM_y,  PsiM_z,
		 cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z,
		 cisvecs,   mu_vecs_x,   mu_vecs_y,  mu_vecs_z,
		 cis_size,  nrots,  dt,  field_h, dumc, 
		 mw_exp,  mw_pop,  mw_wav, &outf, extra_ts, tmpvec);
  else
    single_prop( PsiE,  PsiM_x,  PsiM_y,  PsiM_z,
		 cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z,
		 cisvecs,   mu_vecs_x,   mu_vecs_y,  mu_vecs_z,
		 cis_size,  nrots,  dt,  field_h, dumc, 
		 mw_exp,  mw_pop,  mw_wav, &outf, 0, tmpvec);

  outf << "*******************************************************************************\n";
  outf << "Starting OCT iterations\n";
  sprintf(dumc,"%s.its",argv[2]);
  ofstream it(dumc);
  outf << "Writing target operator value and objective functional value to " << dumc << "\n";
  outf.flush();
  
  double opval;
  if(cis_size < 10000)
    opval = calc_op( target_op, PsiE, ChiE, cis_size);
  else
    opval = calc_op( target_op, PsiE, ChiE, Cis_size);
  
  outf << "Inital target operator value is " << opval << "\n";
  
  double fval  =  calc_field_int(field_h, nrots, dt, shape);
  outf << "Field integral is " << fval << "\n";
  outf << "Total objective functional is " << opval -penal*fval << "\n";
  it << "0 " << opval << " " << opval -penal*fval << " " << fval << "\n"; it.flush();
  outf << "-------------------------------------------------------------------------------\n";
  ofstream it_field;
  for(int cit = 1; cit <= nroit; cit++){
    outf << "Iteration " << cit << "\n"; outf.flush();
    //establish  chi in mu space
//    if(cis_size < 10000){
//      trans_fore(mu_vecs_x, ChiE, ChiM_x, cis_size);
//      trans_fore(mu_vecs_y, ChiE, ChiM_y, cis_size);
//      trans_fore(mu_vecs_z, ChiE, ChiM_z, cis_size);
//    }else{
//      trans_fore(mu_vecs_x, ChiE, ChiM_x, Cis_size);
//      trans_fore(mu_vecs_y, ChiE, ChiM_y, Cis_size);
//      trans_fore(mu_vecs_z, ChiE, ChiM_z, Cis_size);
//    }
   trans_f(cis_size, mu_vecs_x, ChiE, ChiM_x, tmpvec);
   trans_f(cis_size, mu_vecs_y, ChiE, ChiM_y, tmpvec);
   trans_f(cis_size, mu_vecs_z, ChiE, ChiM_z, tmpvec);
    //Backward prop
    for(int x = nrots-1; x >= 0; x--){
      //Calc field_b
      for(int y = 0; y < 3; y++) field_b[y*nrots+x] = 0.;
      if(pfac[0] == 1)  field_b[0*nrots+x] =  calc_field(ChiM_x, PsiM_x, mu_vals_x,  cis_size,  penal)*shape[x];
      if(pfac[1] == 1)  field_b[1*nrots+x] =  calc_field(ChiM_y, PsiM_y, mu_vals_y,  cis_size,  penal)*shape[x];
      if(pfac[2] == 1)  field_b[2*nrots+x] =  calc_field(ChiM_z, PsiM_z, mu_vals_z,  cis_size,  penal)*shape[x];
      
      if(use_mt != 1){
	//Backward one step ( -dt !!!!)
	prop_fore(PsiE,  PsiM_x,  PsiM_y,  PsiM_z, cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z, mu_vecs_x,   mu_vecs_y,  mu_vecs_z,cis_size, x, -dt, nrots, field_h, tmpvec);
	prop_fore(ChiE,  ChiM_x,  ChiM_y,  ChiM_z, cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z, mu_vecs_x,   mu_vecs_y,  mu_vecs_z,cis_size, x, -dt, nrots, field_b, tmpvec);
      }else{
	//Backward mt (-dt)
	prop_fore_mt(PsiE, PsiM_x, PsiM_y, PsiM_z,   cis_vals,   mu_vals_x,   mu_vals_y,  mu_vals_z,    mu_vecs_x,   mu_vecs_y,   mu_vecs_z,
		     ChiE, ChiM_x, ChiM_y, ChiM_z, MTcis_vals, MTmu_vals_x, MTmu_vals_y, MTmu_vals_z, MTmu_vecs_x, MTmu_vecs_y, MTmu_vecs_z,
		     cis_size, x, -dt, nrots,  field_h, field_b);
      }

    }
    
    //REINIT PSI
    for(int x = 0; x< cis_size; x++) PsiE[x] = Psi0[x];

    //establish  psi in mu space
//    if(cis_size < 10000){
//      trans_fore(mu_vecs_x, PsiE, PsiM_x, cis_size);
//      trans_fore(mu_vecs_y, PsiE, PsiM_y, cis_size);
//      trans_fore(mu_vecs_z, PsiE, PsiM_z, cis_size);
//    }else{
//	trans_fore(mu_vecs_x, PsiE, PsiM_x, Cis_size);
//	trans_fore(mu_vecs_y, PsiE, PsiM_y, Cis_size);
//	trans_fore(mu_vecs_z, PsiE, PsiM_z, Cis_size);
//    }
   trans_f(cis_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
   trans_f(cis_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
   trans_f(cis_size, mu_vecs_z, PsiE, PsiM_z, tmpvec);
    //Forward prob
    for(int x = 0; x < nrots; x++){
      //Calc field_h
      for(int y = 0; y < 3; y++) field_h[y*nrots+x] = 0.;
      if(pfac[0] == 1)  field_h[0*nrots+x] =  calc_field(ChiM_x, PsiM_x, mu_vals_x,  cis_size,  penal)*shape[x];
      if(pfac[1] == 1)  field_h[1*nrots+x] =  calc_field(ChiM_y, PsiM_y, mu_vals_y,  cis_size,  penal)*shape[x];
      if(pfac[2] == 1)  field_h[2*nrots+x] =  calc_field(ChiM_z, PsiM_z, mu_vals_z,  cis_size,  penal)*shape[x];

      if(use_mt != 1){
	//foreward one step ( dt !!!!)
	prop_fore(PsiE,  PsiM_x,  PsiM_y,  PsiM_z, cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z, mu_vecs_x,   mu_vecs_y,  mu_vecs_z,cis_size, x, dt, nrots, field_h, tmpvec);
	prop_fore(ChiE,  ChiM_x,  ChiM_y,  ChiM_z, cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z, mu_vecs_x,   mu_vecs_y,  mu_vecs_z,cis_size, x, dt, nrots, field_b, tmpvec);
      }else{
	//forward mt (dt !!!)
	prop_fore_mt(PsiE, PsiM_x, PsiM_y, PsiM_z,   cis_vals,   mu_vals_x,   mu_vals_y,  mu_vals_z,    mu_vecs_x,   mu_vecs_y,   mu_vecs_z,
		     ChiE, ChiM_x, ChiM_y, ChiM_z, MTcis_vals, MTmu_vals_x, MTmu_vals_y, MTmu_vals_z, MTmu_vecs_x, MTmu_vecs_y, MTmu_vecs_z,
		     cis_size, x,  dt, nrots,  field_h, field_b);
      }
    }
    //calc objective func & APPLY OPERATOR
    if(cis_size < 10000)
	opval = calc_op( target_op, PsiE, ChiE, cis_size);
    else
      opval = calc_op( target_op,   PsiE, ChiE, Cis_size);
    
    outf << "Target operator value is " << opval << "\n";
  
    fval  =  calc_field_int(field_h, nrots, dt, shape);
    outf << "Field integral is " << fval << "\n";
    outf << "Total objective functional is " << opval -penal*fval << "\n";
    it << cit << " " << opval << " " << opval -penal*fval << " " << fval << "\n"; it.flush();
    
    sprintf(dumc,"%s_%i.ef",argv[2],cit);
    outf << "Writing current field to " << dumc << "\n";
    it_field.open(dumc);
    it_field.write((char *) &nrots,   sizeof(int));
    it_field.write((char *) &dt,      sizeof(double));
    it_field.write((char *) pfac,      3*sizeof(int));
    it_field.write((char *) field_h,  3*nrots*sizeof(double));
    it_field.write((char *) shape,    nrots*sizeof(double));
    it_field.flush();
    it_field.close();
    
    outf << "-------------------------------------------------------------------------------\n";
    outf.flush();
  }
  outf << "*******************************************************************************\n";
  if(nroit > 0){
    outf << "Final single forward propagation\n";
    sprintf(dumc,"%s_F",argv[2]);
    for(int x = 0; x< cis_size; x++) PsiE[x] = Psi0[x];
    single_prop( PsiE,  PsiM_x,  PsiM_y,  PsiM_z,
		 cis_vals,  mu_vals_x,   mu_vals_y,  mu_vals_z,
		 cisvecs,   mu_vecs_x,   mu_vecs_y,  mu_vecs_z,
		 cis_size,  nrots,  dt,  field_h, dumc, 
		 mw_exp,  mw_pop,  mw_wav, &outf, extra_ts, tmpvec);
  }
  status(&outf); 
}



