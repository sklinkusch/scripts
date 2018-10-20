/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: mk_init_wav.cc                                                         *
 *                                                                              *
 * computes initial wave function in                                            *
 *                                                                              *
 *        1.    Eigenstate basis                                                *
 *        2.    CSF basis                                                       *
 *        3.    Projected-MO  CSFs                                              *
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
#include <stdlib.h>

using namespace std;
#define Complex complex<double>

//Functions
void check_range(int tomuch, int current_number, ofstream *outf);


//Extern Functions
extern int    rem_com(char* filename, char* streamstring, int string_length);
extern void   status(ofstream* outf);
extern void   get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void   read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
			   double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
			   double *Dz);
extern void   read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern int    pos_i_a(int i, int a, int omo, int umo, int llim);
extern void pmv(double* mat, double* vi, double* vo, int nroao);

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "MK_INI_WAV [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";

  //INPUT VARIABLES
  //SYSTEM
  char    bcsfile[256];         //binary system file (CIS space)
  
  //INITAL WAVE FUNCTION
  int     iwav_mode;            //Mode for  inital wave function  creation
                                // 0 -> Eigenstate basis
                                // 1 -> CSF basis
                                // 2 -> Projected-MO CSFs
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  
  ist >> bcsfile;
  outf << "Binary system file (CIS space):  " << bcsfile << "\n";
  
  //CIS-VARIABLES
  outf << "\nREADING CIS-SPACE\n";
  int     nroao;
  int     nroe;                 //Nr of electrons
  int     llim;                 //first MO used for correlation
  int     ulim;                 //last  MO used for correlation
  
  ifstream datf(bcsfile);
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));

  outf << "CIS data read: \n";
  int nrof = ulim - llim +1;
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  int cis_size = omo*umo+1;
  long long int Cis_size =  cis_size;
  outf << "nr of CSF               : " << cis_size << "\n";
  outf << "Allocating Memory\n";

  double* cisvecs     = new double[(long long int) cis_size* (long long int) cis_size];
  double* cis_vals    = new double[cis_size];    //CIS (or CIS(D) eigenvalues)
  double* dumv        = new double[cis_size];    //temporary space

  Complex* wav        = new Complex[cis_size];   //Time dependent wave function (eigen state basis)
  Complex* SEwav      = new Complex[cis_size];   //Time dependent wave function (CSF basis)

  for(int x = 0; x < cis_size; x++)
    wav[x] = SEwav[x] = 0.;

  outf << "Reading data\n";    
  datf.read((char *) cis_vals, cis_size * sizeof(double));
  datf.read((char *) cisvecs, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "-------------------------------------------------------------------------------\n";
  datf.close();
  
  ist >> iwav_mode;
  outf << "iwav_mode is " << iwav_mode << "\n";
  outf << "-------------------------------------------------------------------------------\n";
  
  if(iwav_mode == 0){
    outf << "Constructing inital wave in eigenstate basis\n";
    //Eigenstate basis
    int nroin;                    //Number of initial states
     ist >> nroin;
   for(int j = 0; j < nroin; j++){
    int i_state;
    ist >> i_state;
    ist >> wav[i_state];
    outf << "Initial state is " << i_state << " with coefficient " << wav[i_state] << "\n";
    }
  }


  if(iwav_mode == 1){
    //CSF basis
    int nroin;                 //Number of initial states
     ist >> nroin;
    int i_MO, f_MO;     //Numbers of initial and final MOs
    outf << "Constructing inital wave in CSF basis\n";
    for(int j = 0; j < nroin; j++){
        ist >> i_MO >> f_MO;
        ist >> SEwav[pos_i_a(i_MO, f_MO, omo, umo, llim)];
    outf << "Inital MO is " << i_MO << " with coefficient " << SEwav[pos_i_a(i_MO, f_MO, omo, umo, llim)] << "\n";
    outf << "Final MO is " << f_MO << "\n";}
     for(long long int x = 0; x < Cis_size; x++){
      for(long long int y = 0; y < Cis_size; y++){
	wav[x] += cisvecs[x*Cis_size+y]*SEwav[y];
       }
      }
    
  }
    
  
  if(iwav_mode == 2){
    //Projected-MO CSFs 
    char    sysfile[256];
    char    hfwfile[256];
    int     nroa;
    ist >> sysfile >> hfwfile;
    outf << "Sysfile: " << sysfile << "\n";
    outf << "HF-wavefunction: " << hfwfile << "\n";

    double* i_MO  = new double[nroao];
    double* f_MO  = new double[nroao];
    double* pro_i = new double[nroao];   //Projections on HF orbitals
    double* pro_f = new double[nroao];   //Projections on HF orbitals
    
    for(int x = 0; x < nroao; x++) 
      i_MO[x] = f_MO[x] = pro_i[x] = pro_f[x] =  0.;

    long long  int nrofint;       //Nr of two electron Integrals (not used here)
  
    get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

    outf << "System sizes read from " << sysfile << "\n";
    outf << "Nr of basis functions: " << nroao   << "\n";
    outf << "Nr of atoms:           " << nroa << "\n";
    outf << "Allocating Memory\n";


    double* MOs    = new double[nroao*nroao];
    double* Smat   = new double[nroao*nroao];
    double* dumvec  = new double[nroao];
    
    
    outf.flush();
  
 
    read_sys_1el(sysfile,  MOs, MOs, MOs, MOs,  MOs,  Smat,   MOs,  MOs,  MOs);
    read_wav_HF(hfwfile, nroao, dumvec, MOs);
    
    int nc_i, nc_f;
    int nc;
 
    double pop;

    ist >> nc_i;
    outf << "Reading " << nc_i << " coeffs for the inital Orbital\n";
    for(int x = 0; x < nc_i; x ++){
      ist >> nc;
      ist >> i_MO[nc];
      outf << nc << " " << i_MO[nc] << "\n";
    }


    ist >> nc_f;
    outf << "Reading " << nc_f << " coeffs for the final Orbital\n";
    for(int x = 0; x < nc_f; x ++){
      ist >> nc;
      ist >> f_MO[nc];
      outf << nc << " " << f_MO[nc] << "\n";
    }
    
    outf << "\n\n Calculationg projections on HF orbitals\n";

    pmv(Smat, i_MO , dumvec, nroao);
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++)
	pro_i[x] += dumvec[y]*MOs[x*nroao+y];
    }

    pmv(Smat, f_MO , dumvec, nroao);
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++)
	pro_f[x] += dumvec[y]*MOs[x*nroao+y];
    }
    
    outf << "-------------------------------------------------------------------------------\n";
    for(int x = 0; x < nroao; x++)
      outf << x << "\t" << pro_i[x] << "\t" << pro_f[x] << "\n"; 
      
    pop = 0.; for(int x = 0; x < nroao; x++)  pop += pow( pro_i[x],2);
    outf << "\nTotal population of inital MO in  HF-MOs: " << pop << "\n";
    pop = 0.; for(int x = 0; x < nroe/2; x++)  pop += pow( pro_i[x],2);
    outf << "Population of inital MO in occ-HF-MOs: " << pop << "\n";
    pop = 0.; for(int x = llim; x < nroe/2; x++)  pop += pow( pro_i[x],2);
    outf << "Population of inital MO in active occ-HF-MOs: " << pop << "\n";

    pop = 0.; for(int x = 0; x < nroao; x++)  pop += pow( pro_f[x],2);
    outf << "\nTotal population of final MO in  HF-MOs: " << pop << "\n";
    pop = 0.; for(int x = nroe/2; x < nroao; x++)  pop += pow( pro_f[x],2);
    outf << "Population of final MO in virt-HF-MOs: " << pop << "\n";
    pop = 0.; for(int x = nroe/2; x <= ulim; x++)  pop += pow( pro_f[x],2);
    outf << "Population of inital MO in active virt-HF-MOs: " << pop << "\n";
    
    
    outf << "\nConstructing inital wave function\n";
    
    double wav_norm = 0.;

    for(int x = 1; x < cis_size; x++){
        int i1 = (x-1)/umo+llim;
	int f1 = (x-1)%umo+omo+llim;
	SEwav[x] = pro_i[i1]*pro_f[f1];
	wav_norm += norm(SEwav[x]);
    }
    
    outf << "Initial norm is: " << wav_norm << "\n";
    outf << "Renorming.\n";
    
    for(int x = 0; x < cis_size; x++)
      SEwav[x] /= sqrt(wav_norm);
    
    for(long long int x = 0; x < Cis_size; x++){
      for(long long int y = 0; y < Cis_size; y++){
	wav[x] += cisvecs[x*Cis_size+y]*SEwav[y];
      }
    }
  }

  outf << "-------------------------------------------------------------------------------\n";
  double   curr_time = 0.;
  sprintf(dumc,"%s.iwav",argv[2]);
  ofstream wavf(dumc);
  outf << "Writing initial wave function to " << dumc << "\n";
  wavf.write((char* ) &curr_time, sizeof(double));
  wavf.write((char* ) wav, sizeof(Complex)*cis_size);
  wavf.flush();
      
 
}


void check_range(int tomuch, int current_number, ofstream *outf){
  if( current_number  < 0 ||  current_number >= tomuch){
    *outf << "Invalid input: " << current_number  << "\n";
    outf->flush(); exit(1);
  }
}
