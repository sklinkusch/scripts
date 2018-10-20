/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: istochdft.cc                                                           *
 *                                                                              *
 * RHF program                                                                  *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 *                                                    Stefan Klinkusch 2015     *
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

//Functions

//Extern Functions
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
		     double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		     double *Dz); //, long long int* sortcount, double* intval, 
//		     unsigned short* intnums);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  int   nbasis;           //Nr of configurations
  int  nstates;           //Nr of electronic states
  int     nroe;           //Nr of electrons  (if negative, read in center of mass)
  int     llim;           //first MO used for correlation
  int     ulim;           //last  MO used for correlation
  char sysfile[256];      //binary system file
  char wavfile[256];      //HF-Wavefunction file
  char vecfile[256];      //Eigenvalue and Eigenvector file
  double damp = 1.;

  if(argc != 4){
    cerr << "Need input-file ip output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  cout << "DFT-CSF [CIS3(D) Suite]\n";
  cout << "Reading input from " << argv[1] << "\n\n";
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  ist >> sysfile >> nbasis >> nstates >> nroe;

  double center_of_mass[3]; 

  if(nroe < 0){
    nroe *= -1;
    for(int x = 0; x <3 ; x++) ist  >> center_of_mass[x];
    cout << "Center of mass read in!!!\n";
  }
  
  if(nroe%2 != 0){
    cerr << "Nr of electrons not even !!!\n";
    exit(2);
  }
  
  ist >> llim >> ulim >> wavfile >> vecfile;
  cout << "System data \nNr of electrons: " << nroe << "\n";
  cout << "Number of configurations: " << nbasis << "\n";
  cout << "Number of electronic states: " << nstates << "\n";
//  cout << "CIS data\nLower limit of MOs used for correlation: " << llim << "\n";
//  cout << "Upper limit of MOs used for correlation: " << ulim << "\n";
  cout << "Reading HF-wave function from  " << wavfile << "\n";
  cout << "Reading eigenvalues and vectors from " << vecfile << "\n";
  ifstream vecf(vecfile);

  //SCF VARIABLES 

  //system size
  int    nroao;                 //Nr of basis functions
  int    nroa;                  //Nr of atoms 
  long long  int nrofint;      //Nr of two electron Integrals
  
  get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

  cout << "System sizes read from " << sysfile << "\n";
  cout << "Nr of basis functions: " << nroao   << "\n";
  cout << "Nr of atoms: " << nroa << "\n";
  cout << "Nr of non zero 2el integrals " << nrofint << "\n";
  int omo = nroe/2 - llim;
  int umo = ulim + 1 - llim - omo;
  int cis_size = omo * umo + 1;
  if(cis_size != nbasis){
      cerr << "Wrong limits for the orbitals\n";
      exit(24);
  }
  long long int Cis_size = (long long int) cis_size;
  double* MOs = new double[nroao*nroao];
  double* MOens = new double[nroao];
  read_wav_HF(wavfile, nroao, MOens, MOs);
  cout << "HF-wave function  read from " << wavfile << "\n";

  double*  cisvals   = new double[cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors

  int real_nroao, nros, nroctot;
  vecf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao > nroao){
      cerr << "Wrong nr of AOs in vec file (" << real_nroao << " vs. " << nroao << ")!\n";
      exit(111);
  }
  vecf.read((char *) &nros, sizeof(int));
  if(nros != nstates){
      cerr << "Wrong nr of states in vec file!\n";
      exit(112);
  }
  vecf.read((char *) &nroctot, sizeof(int));
  if(nroctot != nbasis){
      cerr << "Wrong nr of configurations in vec file!\n";
      exit(113);
  }
  vecf.read((char *) cisvals, sizeof(double)*nstates);
  vecf.read((char *) cisvecs, sizeof(double)*nstates*Cis_size);
  cout << "Eigenvalues and vectors read\n";

   //##CALCULATION OF IONIZATION RATES / WRITE OUT TO IRX-FILE / WITHOUT DOUBLES CORRECTIONS
  double* gamma_csf = new double[cis_size];
  double* gamma = new double[nstates];
  double dumm = 0.;
  int MO_r;
  double ip = strtod(argv[2], NULL);
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      gamma_csf[x] = 0.;
      MO_r = (x-1)%umo+omo+llim;
      if(MOens[MO_r] > 0.) gamma_csf[x] = damp*sqrt(MOens[MO_r]);
  }
cout << "Ionization rates for CSFs calculated\n";
  for(int x = 0; x < nstates; x++){
      if(cisvals[x] > ip){
	  dumm = 0.;
	  for(int y = 0; y < nbasis; y++) dumm += pow(cisvecs[x*nbasis+y],2.)*gamma_csf[y];
	  gamma[x] = dumm;
      }else{
	  gamma[x] = 0.;
      }
  }
  sprintf(dumc, "%s.irx", argv[3]);
  ofstream irxf;
  irxf.open(dumc);
  for(int x = 0; x < nstates; x++) irxf << x << " " << gamma[x] << "\n";
  irxf.close();
  cout << "Ionization rates for eigenstates calculated and written\n";
}

