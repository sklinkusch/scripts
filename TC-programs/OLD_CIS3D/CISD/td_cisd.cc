/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: td_cis.cc                                                              *
 *                                                                              *
 * td_cis propagation  program                                                  *
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
extern void   calc_pol(int nrolp, double*  Ao_p, double* phase_p, 
		       int* pol_type, double* pol_vec);
extern void   wirte_las(char* dumc, int plas, double omega, double tpuls, double width, double* Ao, double* phase, 
			double* pol_vec, int pol_type);
extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, int cisd_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, int cisd_size);
extern void   trans_fore(double* mat, Complex* vec_i, Complex* vec_o, long long int cisd_size);
extern void   trans_back(double* mat, Complex* vec_i, Complex* vec_o, long long int cisd_size);
extern double opval(double* vals, Complex* vec, int cisd_size);
extern double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time);

extern void trans_f(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);
extern void trans_b(int np, double* U, Complex* vi, Complex* vo, double* tmpvec);


int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "TDCISD [CISD3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";


  sprintf(dumc,"%s.ten",argv[2]);

  ofstream tenf(dumc);
  tenf.precision(12);
  tenf << "#FORMAT:    $1     $2   \n"
       << "#          time   etot \n";
  outf << "Writing time dependent energy  to " << dumc << "\n";
  
  sprintf(dumc,"%s.fad",argv[2]);

  ofstream fadf(dumc);
  fadf << "#FORMAT:    $1     $2    $3   $4   $5   $6   $7  \n"
       << "#          time   E_x   E_y  E_z  Mu_x Mu_y Mu_z  \n";
  outf << "Writing time dependent field and dipole to " << dumc << "\n";
  
  sprintf(dumc,"%s.wav",argv[2]);
  ofstream wavf(dumc);
  outf << "Writing time dependent PsiE to " << dumc << "\n";
  
  sprintf(dumc,"%s.pop",argv[2]);
  ofstream popf(dumc);
  outf << "Writing time dependent populations to " << dumc << "\n";
  

  
  //INPUT VARIABLES
  
  //SYSTEM
  char    bcsfile[256];      //binary system file
  //  int     use_d = 0;         //if 1 use the 2nd PT corrected energies for propagations

  //Propagation
  int     nrots;          //nr of time steps( if < 0 initial wave function will be read in, otherwise HF-Ground state)
  char    iwavfile[256];  //inital wave function file if nrots < 0;
  int     use_iwav = 0;   //if 1 inital wavefunction is read in
  double  dt;             //time step length
  int     mw_exp;         //analy intervall expectation values
  int     mw_pop;         //analy intervall for populations
  int     mw_wav;         //write out intervall wave function
  
  //Laser
  int nrolp;
 
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

  
  ist >> bcsfile >> nrots;
  if(nrots < 0){ 
    nrots = -1*nrots; 
    ist >> iwavfile;  
    use_iwav = 1;
  }
  ist >> dt >>  mw_exp >> mw_pop >> mw_wav >> nrolp; 
  
  outf << "\nSystem data \nsys-file: " << bcsfile << "\n\n";
  outf << "TD-data: " << "\nNr of time steps: " << nrots << "\ndt: " << dt 
       << "\nExpectation values  written out every " << mw_exp << " time steps" 
       << "\nPopulations         written out every " << mw_pop << " time steps" 
       << "\nMOs                 written out every " << mw_wav << " time steps\n\n";
  if(use_iwav == 1){
    outf << "Reading inital wave functions from " << iwavfile << "\n";
  }
  //  if(use_d == 1)
  //    outf << "Will use the 2nd-PT corrected ex. energies\n";
  //  else
  //    outf << "No correction to CISD-energies used\n";
  outf << nrolp << " sets of laser pulse parameters will be read\n\n";
  outf << "-------------------------------------------------------------------------------\n";


  //CISD READ IN
  ifstream datf(bcsfile);

  int nroao, nroe, llim, ulim;

  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));

  outf << "CISD data read: \n";
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  int nrof = ulim - llim +1;
  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";

  // Loops to determine cisd_size
  int sum2 = 0;
  for (int occ = 1; occ <= omo-1; occ++){
    sum2 += occ;
  }

  int sum1 = 0;
  for (int virt = 1; virt <= umo-1; virt++){
    sum1 += virt;
  }

  int cisd_size = 1 + 2*omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;
  //  int cisd_size = omo*umo+1;
  outf << "nr of CSF               : " << cisd_size << "\n";


  //VARIABLES
  outf << "Allocating " << (long long int) cisd_size * (long long int) cisd_size * 4 *8  + 4* cisd_size *8 << " ??? bytes of memory for CISD operators\n";
  outf.flush();

  double* cisd_vals  = new double[cisd_size];        //eigenvlaues of CI-matrix
  double* mu_vals_x = new double[cisd_size];        //eigenvalues of dipole operator in x
  double* mu_vals_y = new double[cisd_size];        //eigenvalues of dipole operator in x
  double* mu_vals_z = new double[cisd_size];        //eigenvalues of dipole operator in x
  double* tmpvec    = new double[2*cisd_size];
  //  double* corr_vals = new double[cisd_size];        //PT corrections to CISD ex ens

  double* cisdvecs   = new double[(long long int) cisd_size * (long long int) cisd_size];        //eigenvectros of CI-matrix
  double* mu_vecs_x = new double[(long long int) cisd_size * (long long int) cisd_size];        //eigenvectors of dipole operator in x
  double* mu_vecs_y = new double[(long long int) cisd_size * (long long int) cisd_size];        //eigenvectors of dipole operator in y
  double* mu_vecs_z = new double[(long long int) cisd_size * (long long int) cisd_size];        //eigenvectors of dipole operator in z
    
  outf << "Reading data:\n";
   
  datf.read((char *) cisd_vals, cisd_size * sizeof(double));
  datf.read((char *) cisdvecs, (long long int) cisd_size * (long long int) cisd_size * sizeof(double));
  outf << "First uncorrected excitation energies are: " << cisd_vals[1] << " " << cisd_vals[2] << "  ....\n";
  outf << "Highest uncorrected excitation energy  is: " <<  cisd_vals[cisd_size-1] << "\n";

  datf.read((char *) mu_vals_x, cisd_size * sizeof(double));
  datf.read((char *) mu_vecs_x, (long long int) cisd_size * (long long int) cisd_size * sizeof(double));
  outf << "Dipole eigenvalues in x range from " << mu_vals_x[0] << " to " << mu_vals_x[cisd_size-1] << "\n";

  datf.read((char *) mu_vals_y, cisd_size * sizeof(double));
  datf.read((char *) mu_vecs_y, (long long int) cisd_size * (long long int) cisd_size * sizeof(double));
  outf << "Dipole eigenvalues in y range from " << mu_vals_y[0] << " to " << mu_vals_y[cisd_size-1] << "\n";

  datf.read((char *) mu_vals_z, cisd_size * sizeof(double));
  datf.read((char *) mu_vecs_z, (long long int) cisd_size * (long long int) cisd_size * sizeof(double));
  outf << "Dipole eigenvalues in z range from " << mu_vals_z[0] << " to " << mu_vals_z[cisd_size-1] << "\n";
  outf << "-------------------------------------------------------------------------------\n";


  outf.flush();

//   if(use_d == 1){
//     outf << "Reading PT-correction\n";
//     datf.read((char *) corr_vals, cisd_size * sizeof(double));
//     for(int x = 0; x < cisd_size; x++)
//       cisd_vals[x] +=  corr_vals[x];
//     outf << "First CORRECTED excitation energies are: " << cisd_vals[1] << " " << cisd_vals[2] << "  ....\n";
//     outf << "Highest CORRECTED excitation energy  is: " <<  cisd_vals[cisd_size-1] << "\n";
//     outf << "-------------------------------------------------------------------------------\n";
//   }

  datf.close();
  outf << "Reading laserpulses\n";

  //laser 
  double     *omega_p;                //laser frequencies                        nrolp
  double     *t_puls;                 //laser peak times                         nrolp
  double     *width_p;                //laser pulse width                        nrolp
  double     *Ao_p;                   //Amplidute vector                       3*nrolp
  double     *phase_p;                //Phase shifts in x y z                  3*nrolp
  double     *pol_vec;                //Polarization vectros                   3*nrolp  
  
  int doubmem   = 12*nrolp;
  double* dumd = new double[doubmem];int inc = 0;

  //laser
  omega_p = &(dumd[inc]); inc += nrolp; t_puls  = &(dumd[inc]); inc += nrolp; width_p = &(dumd[inc]); inc += nrolp;
  Ao_p    = &(dumd[inc]); inc += 3*nrolp; phase_p = &(dumd[inc]); inc += 3*nrolp; pol_vec = &(dumd[inc]); inc += 3*nrolp;

  int*    pol_type = new int[nrolp];
  
    for(int x = 0; x < nrolp; x++) {
    ist >> omega_p[x] >> t_puls[x] >> width_p[x]
        >> Ao_p[3*x+0] >> phase_p[3*x+0] >> Ao_p[3*x+1] >> phase_p[3*x+1] >> Ao_p[3*x+2] >> phase_p[3*x+2];
  }
  outf << "Nr   Omega   t_peak     width    X (Ao,phase),  Y (Ao,phase),  Z (Ao,phase)\n";
  outf << "..............................................................................\n";
    for(int x = 0; x < nrolp; x++){
    sprintf(dumc,"%i %2.4f %9.2f  %9.2f %2.4f %+2.4f, %2.4f %+2.4f, %2.4f %+2.4f\n",x,  
            omega_p[x],t_puls[x],width_p[x],Ao_p[3*x+0],phase_p[3*x+0],Ao_p[3*x+1],phase_p[3*x+1],Ao_p[3*x+2],phase_p[3*x+2]);
    outf << dumc;
  }
  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Analyzing pulses:\n";
  
  calc_pol( nrolp,  Ao_p, phase_p,  pol_type,  pol_vec);
  
  for(int x = 0 ; x < nrolp; x++){
    outf << "\nPulse " << x << " ";
    if(pol_type[x] == 0 ) 
      outf << " is linearly polarized in (x,y,z)\n\t ";
    else
      outf << " is eliptical polarized pointing along (x,y,z)\n\t " ;
    outf  << pol_vec[3*x+0] << " "  << pol_vec[3*x+1] << " "  << pol_vec[3*x+2] << "\n";  
    sprintf(dumc,"%s-lp%i.dat",argv[2],x);
    outf << "Writing  laser pulse data to " << dumc << "\n";
    wirte_las( dumc, (int) (10.* width_p[x]*omega_p[x]), omega_p[x], t_puls[x], width_p[x], &(Ao_p[3*x]), &(phase_p[3*x]), 
               &(pol_vec[3*x]),  pol_type[x]);
  }

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Init Propagation: \n";
  outf << "Init  Wavefunctions and operators\n";
  
  

  Complex* PsiE;                              //Wave funcion in energy space
  Complex* PsiM_x;                            //Wave function in eigen space of mu_x
  Complex* PsiM_y;                            //Wave function in eigen space of mu_y
  Complex* PsiM_z;                            //Wave function in eigen space of mu_z
  Complex* opH;                               //-i*H*dt in energy space
  
  int complmem = 5*cisd_size; inc = 0;
  Complex* dummc = new Complex[complmem];
  
  PsiE = &(dummc[inc]); inc += cisd_size; PsiM_x = &(dummc[inc]); inc += cisd_size; PsiM_y = &(dummc[inc]); inc += cisd_size; PsiM_z = &(dummc[inc]); inc += cisd_size; 
  opH = &(dummc[inc]); inc += cisd_size;
  
  for(int x = 0; x < cisd_size; x++){
    opH[x]   = exp(-Complex(0.,1.)*cisd_vals[x]*dt);
  }
  
  for(int x = 0; x < cisd_size; x++)
    PsiE[x] = 0.;
  PsiE[0] = 1.;
  double curr_time = 0.;


  if(use_iwav == 1){
    ifstream iw(iwavfile);
    iw.read((char *) &curr_time, sizeof(double));
    iw.read((char *) PsiE, cisd_size*sizeof(Complex));
  }
  
  outf << "Transfroming to dipole spaces\n";
  trans_f(cisd_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
  trans_f(cisd_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
  trans_f(cisd_size, mu_vecs_z, PsiE, PsiM_z, tmpvec); 

  outf << "Initial analysis\n";
  outf << "E mx my mz: " << opval(cisd_vals, PsiE,  cisd_size) << " " <<  opval(mu_vals_x, PsiM_x,  cisd_size)
       << " " <<  opval(mu_vals_y, PsiM_y,  cisd_size)<< " " <<  opval(mu_vals_z, PsiM_z,  cisd_size) << "\n";
  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Starting propagation\n";
  double Efield[3];
  for(int t = 0; t < nrots; t++){
    if(t%mw_exp == 0){
      
      double E_el =  opval(cisd_vals, PsiE,  cisd_size);
      
      outf <<  "pTime: " << curr_time << " " << E_el << "\n";
      outf.flush();
      
      //Fields and dipole
      trans_f(cisd_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
      trans_f(cisd_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
      trans_f(cisd_size, mu_vecs_z, PsiE, PsiM_z, tmpvec); 

      fadf << curr_time << " " << Efield[0] << " " << Efield[1] << " " << Efield[2] << " "
	   <<  opval(mu_vals_x, PsiM_x,  cisd_size)
	   << " " <<  opval(mu_vals_y, PsiM_y,  cisd_size)<< " " <<  opval(mu_vals_z, PsiM_z,  cisd_size) << "\n";
      fadf.flush();
      
      //Energies 
      tenf << curr_time << " " <<  E_el << "\n"; 

      tenf.flush();
    }

    if(t%mw_pop == 0){
      for(int x = 0; x < cisd_size; x++){
        popf << curr_time << " " << cisd_vals[x]  << " " << norm(PsiE[x]) << "\n";
      }
      popf << "\n";
      popf.flush();
    }
    
    if(t%mw_wav == 0){
      wavf.write((char* ) &curr_time, sizeof(double));
      wavf.write((char* ) PsiE, sizeof(Complex)*cisd_size);
      wavf.flush();
    }
    

    //Apply H
    for(int x = 0; x < cisd_size; x++)
      PsiE[x] *= opH[x];
    
    //Calc_Efield 
    for(int x = 0; x < 3; x++) Efield[x] = 0.;
    
    for(int l = 0; l < nrolp; l++){
      for(int x = 0; x < 3; x++){
	Efield[x] += calc_field( omega_p[l], t_puls[l], width_p[l], Ao_p[3*l+x],  phase_p[3*l+x], curr_time);
      }
    }

    //Apply dipole ops
    //X
    if(Efield[0] != 0.){
      trans_f(cisd_size, mu_vecs_x, PsiE, PsiM_x, tmpvec);
      for(int x = 0; x < cisd_size; x++)
	PsiM_x[x] *= exp(Complex(0.,1.)*mu_vals_x[x]*Efield[0]*dt);
      trans_b(cisd_size, mu_vecs_x, PsiM_x, PsiE, tmpvec);
    }
    
    
    //Y
    if(Efield[1] != 0.){
      trans_f(cisd_size, mu_vecs_y, PsiE, PsiM_y, tmpvec);
      for(int x = 0; x < cisd_size; x++)
	PsiM_y[x] *= exp(Complex(0.,1.)*mu_vals_y[x]*Efield[1]*dt);
      trans_b(cisd_size, mu_vecs_y, PsiM_y, PsiE, tmpvec);
    }
    
    //Z
    if(Efield[2] != 0.){
      trans_f(cisd_size, mu_vecs_z, PsiE, PsiM_z, tmpvec);
      for(int x = 0; x < cisd_size; x++)
	PsiM_z[x] *= exp(Complex(0.,1.)*mu_vals_z[x]*Efield[2]*dt);
      trans_b(cisd_size, mu_vecs_z, PsiM_z, PsiE, tmpvec);
    }
    
    curr_time += dt;
  }
    

}
