/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: rhf.cc                                                                 *
 *                                                                              *
 * RHF program                                                                  *
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
#include <stdlib.h>
#include <complex>

using namespace std;

#define Complex complex<double>

//Functions


//Extern Functions
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void status(ofstream* outf);
extern void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
		     double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		     double *Dz, long long int* sortcount, double* intval, 
		     unsigned short* intnums);
extern double calc_ion_rep(int nroa, double* coord, double* charges);
extern void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
			 double* mu_core);
extern double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat, 
		       double* tmpvecs, double* tmpvals);
extern double calc_Sp12(int nroao, double* Smat, double* Sop12, double* tmpmat,
	               double* tmpvecs, double* tmpvals);
extern void diag_Fmat(int nroao, double* Fmat, double* MOs, 
		      double* MOens, double* Som12, double* tmpmat);
extern double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old, 
			      double* MOs, double damp);
extern void build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,
		       double* intval, unsigned short* intnums,
		       long long int* sortcount, long long int nrofint);
extern double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat);
extern double calc_op_1el(int nroao, double* opmat, double* Pmat);
extern void write_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void renorm_MOs(int nroao, double* MOs);
extern void renorm_MO(int nroao, double* MOs, int monr);
extern void symortho_MOsx(int nroao, double* MOs, double* tmat);
extern void mat_vec(int nroao, double* mat, double* vec, double* resvec);
extern void mat_mat_Tn(int np, double* mat_i1, double* mat_i2, double* mat_f, int n);
extern void transpose_mat(int nroao, double* mat, double* transpose);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  int     nroe;           //Nr of electrons (if negative, read in center of mass)
  int     max_cyc;        //maximum nr of SCF cycles 
  double  damp;           //scf damping
  double  tol_E;          //tolerance for total energy
  double  tol_eps;        //tolerance for MO energies
  double  tol_dens;       //tolerance for density
  char sysfile[256];      //binary system file
  int guess = 0;          //guess: 0 -> core guess, 1 -> read guess from file
  char guessfile[256];    //guessfile if guess == 1
  int calc_COM = 0;       // 0 -> calc COM otherwise read in

  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "RHF [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";
  
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
  
  ist >> max_cyc >> damp >> tol_E >> tol_eps >> tol_dens >> guess; 
  outf << "System data \nNr of electrons: " << nroe << "\nsys-file: " << sysfile << "\n\n";
  outf << "SCF data\nMax. nr of cylces: " << max_cyc << "\nDamping: " << damp << "\nTolerance total energy: " << tol_E 
       << "\nTolerance MO-energies: " << tol_eps << "\nTolerance density: " << tol_dens << "\n";
  if(guess == 0) outf << "Core Guess\n";
  if(guess == 1){
    ist >> guessfile;
    outf << "Reading guess from " << guessfile << "\n";
  }
  outf << "-------------------------------------------------------------------------------\n";

  //SCF VARIABLES 

  //system size
  int    nroao;                 //Nr of basis functions
  int    nroa;                  //Nr of atoms 
  long long  int nrofint;      //Nr of two electron Integrals
  
  get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

  outf << "System sizes read from " << sysfile << "\n";
  outf << "Nr of basis functions: " << nroao   << "\n";
  outf << "Nr of atoms: " << nroa << "\n";
  outf << "Nr of non zero 2el integrals " << nrofint << "\n";
  outf << "Allocating Memory\n";

  //atoms
  double*        coord;         //atomic coordinats               3*nroa
  double*        charges;       //atomic charges                    nroa
  double*        mass;          //atomic masses                     nroa

  //one electron mat&vecs
  double*        Smat;          //Overlap matrix S                  nroao*nroao
  double*        Hmat;          //one electron Hamiltionian         nroao*nroao
  double*        Tmat;          //Kinetic energy operator           nroao*nroao
  
  double*        Som12;         //S^-1/2                            nroao*nroao
  double*        Pmat;          //density matrix                    nroao*nroao
  double*        Pmat_old;      //Old density matrix                nroao*nroao

  double*        Fmat;          //Fock matrix                       nroao*nroao
  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
  double*        MOens_old;     //Old MO energies;                  nroao


  //Temporary memory spaces
  double*  tmpmat1;                //                               nroao*nroao
  double*  tmpmat2;                //                               nroao*nroao
  double*  tmpvecs;                //                               nroao*nroao

  double*  tmpvals;                //                               nroao
 
  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+14*nroao*nroao+3*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]); inc += 3*nroa; charges = &(dumd[inc]); inc += nroa; mass  = &(dumd[inc]); inc += nroa;
  Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao; Tmat = &(dumd[inc]); inc += nroao*nroao;
  Som12 = &(dumd[inc]); inc += nroao*nroao; Pmat = &(dumd[inc]); inc += nroao*nroao;Pmat_old = &(dumd[inc]); inc += nroao*nroao;
  Fmat = &(dumd[inc]); inc += nroao*nroao; MOs = &(dumd[inc]); inc += nroao*nroao; Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao; 
  
  MOens = &(dumd[inc]); inc += nroao; MOens_old = &(dumd[inc]); inc += nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; tmpmat2 = &(dumd[inc]); inc += nroao*nroao; tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
  tmpvals = &(dumd[inc]); inc += nroao;

  //MEMORY ALLOCATION for two electron values 
  outf << "Need " << nrofint*(sizeof(double)+sizeof(unsigned short)*4) << " bytes for two electron data\n";
  outf.flush();
  
  double*         intval        = new double[nrofint];                     //two electron integrals
  unsigned short* intnums       = new unsigned short[nrofint*4];           //two electron indices
  long long int   sortcount[4];                                            //num of two electron Integrals in each perm. type

  
  //Read system data 
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);
  

  //INIT SCF
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Iinit. SCF:\n";
  outf << "...............................................................................\n";
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
  outf << "Minimal eigenvalue of S is " << calc_S12( nroao, Smat, Som12, tmpmat1, tmpvecs, tmpvals) <<"\n";
  
  if(guess == 0){
    outf << "Performing core guess\n";
    for(int x = 0; x < nroao*nroao; x++)
      Fmat[x] = Hmat[x];
    diag_Fmat( nroao, Fmat, MOs, MOens, Som12, tmpmat1);
  }
  if(guess == 1){
    read_wav_HF(guessfile, nroao, MOens, MOs);
    outf << "Inital guess read from " << guessfile << "\n";
  }

  outf << "Building initial density and Fock matrix\n";
  
  build_Pmat_dscf( nroao, nroe,  Pmat,  Pmat_old,  MOs, 0.);
  build_Fmat(nroao,  Fmat, Pmat, Hmat, intval, intnums, sortcount, nrofint);


  outf << "Guess Energy: " << Calc_e_el( nroao, Fmat, Pmat, Hmat) + ion_rep << "\n";
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Entering SCF:\n";
  outf << "cyc         E             dE             dD               deps\n";
  outf.flush();

  double currE = 0. , maxd_dens, maxd_eps =0., oldE, dE;
  
  oldE = Calc_e_el( nroao, Fmat, Pmat, Hmat);
  for(int x = 0; x < nroao; x++) MOens_old[x] = MOens[x];

  int conv = 0;

  for(int c = 0; c < max_cyc; c++){
    
    diag_Fmat( nroao, Fmat, MOs,  MOens,  Som12, tmpmat1);
    maxd_dens = build_Pmat_dscf( nroao, nroe,  Pmat,  Pmat_old,  MOs, damp);
    build_Fmat(nroao,  Fmat, Pmat, Hmat, intval, intnums, sortcount, nrofint);

    currE =  Calc_e_el( nroao, Fmat, Pmat, Hmat);
    dE    =  currE-oldE;

    for(int x = 0; x < nroao; x++){
      if(fabs(MOens[x]-MOens_old[x]) >  fabs(maxd_eps))maxd_eps = MOens[x]-MOens_old[x];
      MOens_old[x] = MOens[x];
    }
    sprintf(dumc,"%.12f %.12f %.12f %.12f\n",currE+ion_rep, dE, maxd_dens, maxd_eps);
    outf << c << " " << dumc;
    outf.flush();

    if(fabs(maxd_eps) < tol_eps && fabs(dE) < tol_E && fabs(maxd_dens) < tol_dens){
      c = max_cyc+10;
      conv = 1;
    }
    oldE = currE;
    maxd_eps = 0.; 
  }
    
  if(conv == 1){
    double ekin = calc_op_1el( nroao, Tmat, Pmat);
    outf << "SCF Converged!!!!\n";
    outf << "-------------------------------------------------------------------------------\n";
    outf << "Final electronic energy: " << currE << "\n";
    outf << "Kinetic energy: " << ekin  << "\n";
    outf << "Total energy: " << currE+ion_rep << "\n\n";
    outf << "-V/T : " << -(currE+ion_rep-ekin)/ekin << "\n\n";

    

    outf << "MO-energies:\n";
    for(int x = 0; x < nroao; x++){
      sprintf(dumc,"%+.5f ",MOens[x]);
      outf << dumc << "\t";
      if((x+1)%5==0) outf << "\n";
    }
    
    outf << "\n\nDipole moment (AU): " 
	 <<  -calc_op_1el( nroao, Dx, Pmat) +mu_core[0] << " " 
	 <<  -calc_op_1el( nroao, Dy, Pmat) +mu_core[1] << " " 
	 <<  -calc_op_1el( nroao, Dz, Pmat) +mu_core[2] << "\n";
    outf << "Dipole moment (Debye): " 
	 <<  -(calc_op_1el( nroao, Dx, Pmat) -mu_core[0])*2.5417462 << " " 
	 <<  -(calc_op_1el( nroao, Dy, Pmat) -mu_core[1])*2.5417462 << " " 
	 <<  -(calc_op_1el( nroao, Dz, Pmat) -mu_core[2])*2.5417462 << "\n";
  }
  else{
    outf << "SCF unconverged !! To many cycles !!!!!\n";
  }
  sprintf(dumc,"%s.hfw",argv[2]);
  outf << "\n\nWriting final HF-wavefunction to " << dumc << "\n";

  // symmortho_MO with S^1/2
  double* Sop12 = new double[nroao*nroao];
  double Sop_minval = calc_Sp12(nroao, Smat, Sop12, tmpmat1, tmpvecs, tmpvals); 
  double* Zvec = new double[nroao*nroao];
  for(int i = 0; i < nroao*nroao; i++) Zvec[i] = MOs[i];
  symortho_MOsx(nroao, Zvec, Sop12);
  // modified Gramm-Schmidt
  double sigma = 0.;
  for(int i = 0; i < nroao; i++){
   renorm_MO(nroao, Zvec,i);  
   mat_mat_Tn(nroao, Som12, Zvec, MOs, i);
   for(int j = i+1; j < nroao; j++){
    sigma = 0.;
    for(int k = 0; k < nroao; k++) sigma += Zvec[i*nroao+k] * Zvec[j*nroao+k];
    for(int k = 0; k < nroao; k++) Zvec[j*nroao+k] -= sigma * Zvec[i*nroao+k];
   }
  }
  write_wav_HF(dumc, nroao, MOens, MOs);
  char dumfile[512];
  sprintf(dumfile, "%s.olp", argv[2]);
  ofstream olpf;
  olpf.open(dumfile);
  for(int i = 0; i < nroao; i++){
   for(int j = 0; j < nroao; j++){
    olpf << Smat[i*nroao+j] << " ";
   }
   olpf << "\n";
   olpf.flush();
  }
  olpf.close();
  outf << "Execution ended on/at:\n";
  status(&outf);
}
