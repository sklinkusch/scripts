/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: td_rhf.cc                                                              *
 *                                                                              *
 * TD_RHF program                                                               *
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
void calc_ruku4_fields(int nrolp, double*  omega_p, double*  t_puls, double* width_p, double* Ao_p, double* phase_p,
		       double curr_time, double dt, double* Efields);

//Extern Functions
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void status(ofstream* outf);
extern void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
                     double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
                     double *Dz, long long int* sortcount, double* intval, 
                     unsigned short* intnums);
extern void calc_pol(int nrolp, double*  Ao_p, double* phase_p, 
		     int* pol_type, double* pol_vec);
extern void wirte_las(char* dumc, int plas, double omega, double tpuls, double width, double* Ao, double* phase, 
		      double* pol_vec, int pol_type);
extern double calc_ion_rep(int nroa, double* coord, double* charges);
extern void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
                         double* mu_core);
extern double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat, 
                       double* tmpvecs, double* tmpvals);
extern double calc_Sp12(int nroao, double* Smat, double* So12, double* tmpmat, 
			double* tmpvecs, double* tmpvals);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void build_Pmat(int nroao, int nroe, Complex* Pmat, Complex* MOs);
extern void build_Fmat(int nroao, Complex* Fmat, Complex* Pmat, double* Hmat,
		       double* intval, unsigned short* intnums,
		       long long int* sortcount, long long int nrofint);
extern Complex calc_e_el(int nroao, Complex* Fmat, double* Hmat, Complex* Pmat);
extern void    calc_mu(int nroao,  Complex* Pmat, 
		       Complex* mu, double *mu_x, double* mu_y, double* mu_z);
extern double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time);
extern void ruku4(int nroao, int nroe, double dt, double* Hmat,  Complex* Pmat,  Complex* Fmat, 
		  Complex* MOs, Complex* MOs_dt, Complex *alt_MOs, Complex *sum_MOs, 
		  Complex* dumvec, Complex* dummat,
		  double* Som12, double* So12, double* E_field, double *Dx, double* Dy, double* Dz,
		  double* intval, unsigned short* intnums,
		  long long int* sortcount, long long int nrofint, int nrofc);
extern Complex calc_one_p_op(int nroao,  Complex* Pmat, double *op);
extern void check_id_td(double* S, Complex* MOs, Complex*  dummat, int nroao, int nroe,  double* norms);
extern void calc_pops(double* S, double* initMOs, Complex* MOs, Complex*  dumvec, int nroao, int nroe,  double* pops);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  
  //SYSTEM
  int     nroe;           //Nr of electrons (if negative read in center of mass)
  int calc_COM = 0;       // 0 -> calc COM otherwise read in
  char sysfile[256];      //binary system file
  char wavfile[256];      //MO are read from

  //Propagation
  int     nrots;          //nr of time steps
  double  dt;             //time step length
  int     mw_exp;         //analy intervall expectation values
  int     mw_pop;         // "       "      populations
  int     mw_wav;         //write out intervall wave function

  int     nrofc=0;        //nr of frozen core orbitals (not changed in propagation)
  
  
  //Laser
  int nrolp;
  
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "TDRHF [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";

  sprintf(dumc,"%s.ten",argv[2]);
  ofstream tenf(dumc);
  tenf.precision(12);
  tenf << "#FORMAT:    $1     $2    $3  $4     $5          $6       $7            $8            $9              $10  \n"
       << "#          time   etot ekin epot real(norm) imag(norm) max_err_n(MO) av_err_n(MO) max_err_orhot av_err_ortho \n";
  outf << "Writing time dependent energies and norms to " << dumc << "\n";
  
  sprintf(dumc,"%s.fad",argv[2]);
  ofstream fadf(dumc);
  fadf << "#FORMAT:    $1     $2    $3   $4   $5   $6   $7  \n"
       << "#          time   E_x   E_y  E_z  Mu_x Mu_y Mu_z   (imaginary parts of mu)\n";
  outf << "Writing time dependent field and dipole to " << dumc << "\n";


  sprintf(dumc,"%s.pop",argv[2]);
  ofstream popf(dumc);
  popf << "#FORMAT:    $1     $2    $3    $nroao+1     \n"
       << "#          time   Pop_1  ...    sum of pops \n";
  outf << "Writing time dependent populations to " << dumc << "\n";
  
  sprintf(dumc,"%s.wav",argv[2]);
  ofstream wavf(dumc);
  outf << "Writing time dependent MOs to " << dumc << "\n";
  
  
  outf << "-------------------------------------------------------------------------------\n";
  outf.flush();
  
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

  
  ist >> sysfile >> nroe;
  
  double center_of_mass[3]; 

  if(nroe < 0){
    nroe *= -1;
    calc_COM = 1;
    for(int x = 0; x <3 ; x++) ist  >> center_of_mass[x];
    outf << "Center of mass read in!!!\n";
  }

  if(nroe%2 != 0){
    cerr << "Nr of electrons not even !!!\n";
    exit(2);
  }
  
  ist >> nrots >> dt >> wavfile >> mw_exp >> mw_pop >> mw_wav >> nrolp; 
  outf << "System data \nNr of electrons: " << nroe << "\nsys-file: " << sysfile << "\n\n";
  outf << "TD-data: " << "\nNr of time steps: " << nrots << "\ndt: " << dt 
       << "\nMOs read from " << wavfile 
       << "\nExpectation values  written out every " << mw_exp << " time steps" 
       << "\nPopulations         written out every " << mw_pop << " time steps" 
       << "\nMOs                 written out every " << mw_wav << " time steps\n\n";
  outf << nrolp << " sets of laser pulse parameters will be read\n\n";

  //HF VARIABLES 

  //system size
  int    nroao;                 //Nr of basis functions
  int    nroa;                  //Nr of atoms 
  long long  int nrofint;       //Nr of two electron Integrals
  
  get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

  outf << "System sizes read from " << sysfile << "\n";
  outf << "Nr of basis functions: " << nroao   << "\n";
  outf << "Nr of atoms: " << nroa << "\n";
  outf << "Nr of non zero 2el integrals " << nrofint << "\n\n";

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Allocating Memory\n";

  //double Memomry 
  double      *dumd;

  //atoms
  double      *coord;         //atomic coordinats               3*nroa
  double      *charges;       //atomic charges                    nroa
  double      *mass;          //atomic masses                     nroa

  //AO integrals
  double      *Smat;                  //overlap                                  nroao*nroao
  double      *Hmat;                  //core hamiltonian                         nroao*nroao
  double      *Tmat;                  //kinetic energy                           nroao*nroao
  double      *Som12;                 //S^-1/2                                   nroao*nroao
  double      *So12;                  //S^1/2                                    nroao*nroao

  double      *Dx, *Dy, *Dz;          //Dipole Matricies                       3*nroao*nroao  
  double      *svecs;                 //Eigen vectors of overlap matrix          nroao*nroao
  double      *initMOs;               //Read in MOs (real)                       nroao*nroao

  double      *dummat;                //temporary space                          nroao*nroao
  double      *dumvec;                // "                                       nroao*nroao

  double      *dumvals;               //eigen values                             nroao
  double      *initMOens;             //read in MOens (real)                     nroao

  //laser 
  double     *omega_p;                //laser frequencies                        nrolp
  double     *t_puls;                 //laser peak times                         nrolp
  double     *width_p;                //laser pulse width                        nrolp
  double     *Ao_p;                   //Amplidute vector                       3*nrolp
  double     *phase_p;                //Phase shifts in x y z                  3*nrolp
  double     *pol_vec;                //Polarization vectros                   3*nrolp  
  
  int doubmem   = 5*nroa+12*nroao*nroao+2*nroao+12*nrolp;
  outf << "\nNeed " << 8*doubmem << " bytes for AO integrals etc. pp \n";
  outf.flush();
  
  dumd  = new double[doubmem]; int inc = 0;

  //atoms
  coord = &(dumd[inc]); inc += 3*nroa; charges = &(dumd[inc]); inc += nroa; mass  = &(dumd[inc]); inc += nroa; 

  //AO integrals
  Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao; Tmat = &(dumd[inc]); inc += nroao*nroao;
  Som12 = &(dumd[inc]); inc += nroao*nroao; So12 = &(dumd[inc]); inc += nroao*nroao;

  Dx = &(dumd[inc]); inc += nroao*nroao; Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao; 
  svecs = &(dumd[inc]); inc += nroao*nroao; initMOs =  &(dumd[inc]); inc += nroao*nroao;

  dummat = &(dumd[inc]); inc += nroao*nroao; dumvec = &(dumd[inc]); inc += nroao*nroao;
  
  initMOens = &(dumd[inc]); inc += nroao; dumvals = &(dumd[inc]); inc += nroao;
  
  //laser
  omega_p = &(dumd[inc]); inc += nrolp; t_puls  = &(dumd[inc]); inc += nrolp; width_p = &(dumd[inc]); inc += nrolp;
  Ao_p    = &(dumd[inc]); inc += 3*nrolp; phase_p = &(dumd[inc]); inc += 3*nrolp; pol_vec = &(dumd[inc]); inc += 3*nrolp;

  int*    pol_type = new int[nrolp];

  //Complex memory 
  
  Complex     *dumC;

  Complex     *Pmat;                  //density matrix                           nroao*nroao
  Complex     *Fmat;                  //Fockmatrix                               nroao*nroao
  Complex     *cdummat;               //complex temp space                       nroao*nroao      
  Complex     *cdumvec;               // "        "   "                          nroao*nroao

  Complex     *MOs;                   //Fiele for  MO Coefficents                nroao*nroe/2
  Complex     *MOs_dt;                //time derivatives                         nroao*nroe/2
  Complex     *alt_MOs, *sum_MOs;     //neede for Ruku4 Integrator             2*nroao*nroe/2  

  Complex     *MOens;                 //MO Energies                              nroe/2 

  Complex     *mu;                    //Dipole Vector                            3 
  
  int ctotmem = 4*nroao*nroao+4*nroao*nroe/2 +nroe/2 +3;
  outf << "Need " << 16*ctotmem << " bytes for time dependent (complex) vectors \n";
  outf.flush();

  
  dumC = new Complex[ctotmem]; inc = 0;
  Pmat = &( dumC[inc]); inc += nroao*nroao;  Fmat = &( dumC[inc]); inc += nroao*nroao; cdummat = &( dumC[inc]); inc += nroao*nroao;
  cdumvec    = &( dumC[inc]); inc += nroao*nroao;

  MOs = &( dumC[inc]); inc += nroao*nroe/2; MOs_dt = &( dumC[inc]); inc += nroao*nroe/2;
  alt_MOs = &( dumC[inc]); inc += nroao*nroe/2; sum_MOs = &( dumC[inc]); inc += nroao*nroe/2;

  MOens = &( dumC[inc]); inc += nroao;       

  mu  = &( dumC[inc]); inc += 3;           

  //MEMORY ALLOCATION for two electron values 
  outf << "Need " << nrofint*(sizeof(double)+sizeof(unsigned short)*4) << " bytes for two electron data\n";
  outf.flush();
  
  double*         intval        = new double[nrofint];                     //two electron integrals
  unsigned short* intnums       = new unsigned short[nrofint*4];           //two electron indices
  long long int   sortcount[4];                                            //num of two electron Integrals in each perm. type
  
  outf << "-------------------------------------------------------------------------------\n";
  
  outf << "Reading sysfile: " <<  sysfile << "\n\n";
  outf.flush();

  //Read system data 
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);
  
  outf << "Reading laser pulse parameters\n";
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
  outf << "Initializing td calculation:\n";
  outf << "...............................................................................\n";

  ist >> nrofc;
  outf << "Nr of frozen core orbital in propagation is: " << nrofc << "\n";
  
  outf << "Ion Cores: coordinates mass charge\n";
  int count = 0;
  for(int x = 0; x < nroa; x++){
    outf << x << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%.5f", mass[x]);           outf  << dumc << "\t";
    sprintf(dumc,"%.0f", charges[x]);        outf  << dumc << "\n";     
  }
  
  outf.precision(10);


  double mu_core[3];
  double ion_rep =  calc_ion_rep( nroa, coord, charges);
  
  if(calc_COM == 0)
    calc_center_of_mass( nroa, coord,  mass,  center_of_mass);
  calc_mu_core( nroa, coord, charges, center_of_mass, mu_core);
  
  outf << "\n\nIon repulsion is: " << ion_rep  << "\n";
  outf << "Center of mass (x,y,z): " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << "\n";
  outf << "Core dipole moment (at C.o.M.): " << mu_core[0] << " " << mu_core[1] << " " << mu_core[2] << "\n";
  outf << "...............................................................................\n";

  
  outf << "Init operators:\n";
  
  outf << "Calculating S^-1/2\n";
  outf << "Minimal eigenvalue of S is " <<  calc_S12( nroao,  Smat,  Som12, dummat,  dumvec, dumvals);
  outf << "\n\nCalculating S^+1/2\n";
  outf << "Minimal eigenvalue of S is " <<  calc_Sp12( nroao,  Smat,  So12, dummat,  dumvec, dumvals) << "\n\n";
  
  read_wav_HF(wavfile, nroao, initMOens, initMOs);
  outf << "Inital wave function read from " << wavfile << "\n";

  outf << "Setting up initial wave function\n";
  
  for(int x = 0; x < nroe/2*nroao;x++) MOs[x] = initMOs[x];
  
  outf << "Bulfing inital P and F matrix\n";
  
  build_Pmat( nroao,  nroe, Pmat,  MOs);
  build_Fmat( nroao,  Fmat, Pmat,  Hmat, intval, intnums, sortcount, nrofint);
  
  
  Complex E_el =  calc_e_el( nroao, Fmat, Hmat, Pmat);
  outf << "Initial energy: " << real(E_el)+ion_rep << "  (imag: " << imag(E_el) << ")\n";
  calc_mu(nroao,  Pmat,  mu, Dx, Dy, Dz);
  
  outf << "Initial electronic dipole: (" << -real(mu[0])+mu_core[0] << "," << -real(mu[1])+mu_core[1] << "," << -real(mu[2])+mu_core[2] << ")  [imag ( " 
       << imag(mu[0]) << "," << imag(mu[1]) << "," << imag(mu[2]) << ") ]\n";
  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Starting propagation\n";
  outf << "...............................................................................\n";
  
  outf << "...............................................................................\n";

  double Efields[9];
  double td_norms[6];
  
  double curr_time = 0.;
  //DUMMY Call for printing !!
  calc_ruku4_fields(nrolp, omega_p,  t_puls,  width_p,  Ao_p,  phase_p, curr_time,  dt,  Efields);

  for(int x = 0; x < nrots; x++){
    
    if(x%mw_exp == 0){
      check_id_td(Smat,  MOs, cdummat, nroao, nroe,  td_norms);
      

      build_Pmat( nroao,  nroe, Pmat,  MOs);
      build_Fmat( nroao,  Fmat, Pmat,  Hmat, intval, intnums, sortcount, nrofint);
      E_el =  calc_e_el( nroao, Fmat, Hmat, Pmat);
      
      outf <<  "pTime: " << curr_time << " " << real(E_el)+ion_rep << "  (imag: " << imag(E_el) << ") " 
	   <<   "  err_norm " << td_norms[0] -1. << "   max_err_orotho " << td_norms[5] << "\n";
      outf.flush();
      
      //Fields and dipole
      calc_mu(nroao,  Pmat,  mu, Dx, Dy, Dz);
      fadf << curr_time << " " << Efields[6] << " " << Efields[7] << " " << Efields[8] << " "
	   << -real(mu[0])+mu_core[0] << " " << -real(mu[1])+mu_core[1] << " " << -real(mu[2])+mu_core[2] << " " 
	   << imag(mu[0]) << " " << imag(mu[1]) << " " << imag(mu[2]) << "\n";
      fadf.flush();
      
      //Energies and norms
      Complex ekin = calc_one_p_op( nroao, Pmat, Tmat);
      tenf << curr_time << " " <<  real(E_el)+ion_rep << " " << real(ekin) << " " << real(E_el)+ion_rep-real(ekin) << " " 
	   << td_norms[0] << " " << td_norms[1] << " " <<  td_norms[2] << " " <<  td_norms[3] << " " <<  td_norms[4] << " " 
	   << td_norms[5] << "\n"; 

      tenf.flush();
      
    }

    if(x%mw_pop == 0){
      calc_pops( Smat, initMOs, MOs, cdumvec,  nroao, nroe, dumvals); //pops dumvals
      for(int x = 0; x < nroao; x++){
	sprintf(dumc,"   orb%.3i\n",x);
	popf << curr_time << " " << initMOens[x] << " " << dumvals[x] << " " << x <<  dumc;
      }
      popf << "\n";
      popf.flush();

    }
    

    if(x%mw_wav == 0){
      wavf.write((char* ) &curr_time, sizeof(double));
      wavf.write((char* ) MOs, sizeof(Complex)*nroao*nroe/2);
      wavf.flush();
    }

    calc_ruku4_fields(nrolp, omega_p,  t_puls,  width_p,  Ao_p,  phase_p, curr_time,  dt,  Efields);
    ruku4(nroao, nroe, dt, Hmat, Pmat,  Fmat, 
	  MOs, MOs_dt, alt_MOs, sum_MOs, 
	  cdumvec,  cdummat,
	  Som12, So12, Efields, Dx, Dy, Dz,
	  intval, intnums,
	  sortcount,  nrofint, nrofc);
    
    curr_time += dt;
  }
}




void calc_ruku4_fields(int nrolp, double*  omega_p, double*  t_puls, double* width_p, double* Ao_p, double* phase_p,
		       double curr_time, double dt, double* Efields){
  for(int x = 0; x < 9; x++) Efields[x] = 0.;

  for(int l = 0; l < nrolp; l++){
    double dt2 = dt/2.;
    for(int x = 0; x < 3; x++){
      Efields[x*3+0] += calc_field( omega_p[l], t_puls[l], width_p[l], Ao_p[3*l+0],  phase_p[3*l+0], curr_time+x*dt2);
      Efields[x*3+1] += calc_field( omega_p[l], t_puls[l], width_p[l], Ao_p[3*l+1],  phase_p[3*l+1], curr_time+x*dt2);
      Efields[x*3+2] += calc_field( omega_p[l], t_puls[l], width_p[l], Ao_p[3*l+2],  phase_p[3*l+2], curr_time+x*dt2);
    }
  }
}
  

 

