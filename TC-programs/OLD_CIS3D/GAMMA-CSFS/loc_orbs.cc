/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: loc_orbs.cc                                                            *
 *                                                                              *
 * localizes orbitals in z                                                      *
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
#include <stdlib.h>

using namespace std;

extern int    rem_com(char* filename, char* streamstring, int string_length);
extern void   status(ofstream* outf);
extern void   get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void   read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
			   double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
			   double *Dz);
extern void   read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void   pmv(double* mat, double* vi, double* vo, int nroao);
extern void   diag_mat(int nroao, double* mat, double* vals, double* vecs);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  int     llim;           //first occ  MO used for loc in z
  int     ulim;           //last  unocc MO used for loc in z
  
  char sysfile[256];      //binary system file
  char wavfile[256];      //HF-Wavefunction file
  
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "LOC_ORBS [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";

  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  
  ist >> sysfile;
  ist >> llim >> ulim >> wavfile;

  outf << "sys-file: " << sysfile << "\n\n";
  outf << "LOC data\nLower limit of MOs used for localization: " << llim << "\n";
  outf << "Upper limit of MOs used for localization: " << ulim << "\n";
  outf << "Reading HF-wave function from  " << wavfile << "\n";

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
  
  double*        Pmat;          //density matrix                    nroao*nroao
  double*        Pmat_Mul;      //Mulliken matrix                   nroao*nroao

  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
  double*        tdPOP_mo;      //td MO-populations                 nroao

  //Temporary memory spaces
  double*  tmpmat1;             //                                  nroao*nroao
  double*  tmpmat2;             //                                  nroao*nroao
  double*  tmpvecs;             //                                  nroao*nroao

  double*  tmpvals;             //                                  nroao

  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+12*nroao*nroao+3*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]); inc += 3*nroa;     charges = &(dumd[inc]); inc += nroa;      mass  = &(dumd[inc]); inc += nroa;
  Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao;  Tmat = &(dumd[inc]); inc += nroao*nroao;
  Pmat = &(dumd[inc]); inc += nroao*nroao; Pmat_Mul = &(dumd[inc]); inc += nroao*nroao;
  MOs = &(dumd[inc]); inc += nroao*nroao;  Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao;   Dz = &(dumd[inc]); inc += nroao*nroao; 
  
  MOens = &(dumd[inc]); inc += nroao;      tdPOP_mo = &(dumd[inc]); inc += nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; 
  tmpmat2 = &(dumd[inc]); inc += nroao*nroao; 
  tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
  tmpvals = &(dumd[inc]); inc += nroao;
 
  read_sys_1el(sysfile,  coord, charges, mass, Hmat,  Tmat,  Smat,   Dx,  Dy,  Dz);
  read_wav_HF(wavfile, nroao, MOens, MOs);

 
  outf << "HF reference data:\n";
  outf << "...............................................................................\n";
  outf << "Ion Cores: coordinates mass charge\n";
  int count = 0;
  for(int x = 0; x < nroa; x++){
    outf << x << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%.5f", mass[x]);            outf  << dumc << "\t";
    sprintf(dumc,"%.0f", charges[x]);         outf  << dumc << "\n";     
  }
  
  outf.precision(10);

  outf << "\nMO-energies:\n";
  for(int x = 0; x < nroao; x++){
    sprintf(dumc,"%+.5f ",MOens[x]);
    outf << dumc << "\t";
    if((x+1)%5==0) outf << "\n";
  }
  outf << "\n...............................................................................\n";
  outf << "Buliding <z> in orbital space between " << llim << ", " << ulim << "\n";
  
  int lsize = ulim-llim+1;
  
  for(int x = 0; x < lsize ; x++){ 
    pmv(Dz, &(MOs[(x+llim)*nroao]), tmpvecs, nroao);
    for(int y = 0; y < lsize; y++){
      tmpmat1[x*lsize+y] = 0.;
      for(int z = 0; z < nroao; z++)
	tmpmat1[x*lsize+y] +=  tmpvecs[z]*MOs[(y+llim)*nroao+z];
    }
  }
  diag_mat(lsize, tmpmat1,  tmpvals, tmpvecs);
  
  outf << "Eigenvalues  z, expectation values for epsilon:\n";
  for(int x = 0; x <  lsize ; x++){
    double eps = 0.;
    for(int y = 0; y < lsize; y++)
      eps += pow(tmpvecs[x*lsize+y],2)*MOens[y+llim];
    outf << x << " " << tmpvals[x] << " " << eps << "\n";
  }

  outf << "\n\nOutput of individual Orbitals:\n";
  for(int x = 0; x < lsize ; x++){
    outf << "Orbital " << x << ":\n";
    for(int y = 0; y < nroao; y++){
      double coeff = 0.;
      for(int z = 0; z < lsize; z++)
	coeff += MOs[(z+llim)*nroao+y]*tmpvecs[x*lsize+z];
      outf << y << "\t" << coeff << "\n";
    }
  }
}
