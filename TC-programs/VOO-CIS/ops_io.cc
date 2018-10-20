/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_io.cc                                                              *
 *                                                                              *
 * contains io routines                                                         * 
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
#include <stdlib.h>


using namespace std;


//Functions
void status(ofstream* outf);                                                                       // prints status(date, time, host)
int  rem_com(char* filename, char* streamstring, int string_length);                               // reads from input file, removes comments, and writes data to a string
//void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);                   // reads #AOs, #atoms, and #integrals from a sysfile generated by mk_in
/*void read_sys(char* sysfile, double* coord, double* charges, double* mass,                         // reads atom coordinates, core charges and masses, the Hamiltonian matrix, kinetic energy matrix, overlap matrix, dipole matrix elements, two-electron integrals and 
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy,                   // indices from the sysfile
	      double *Dz, long long int* sortcount, double* intval, 
	      unsigned short* intnums);*/
/*void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass,                     // another reader for the sysfile
		  double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		  double *Dz);*/
//void write_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                           // writes the #AOs, MO energies, and MOs to a wavfile
void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                            // reads the #AOs, MO energies, and MOs from a wavfile
//void output_matrix(double* mat,int nr_of_col, int nrop, ofstream *outf);                           // prints a matrix to an output filestream
//void read_matrix(char* matfile, int nros, double* cisvals, double* cisvecs);                       // reads the CI matrix from a binary file
//void read_bin(char* binfile, int nros, double& cisEgs, double* cisvals);                           // reads the #states, ground state energy and CIS excitation energies

//Extern Functions

//extern "C" void dafrd_(double* V, long long int* LEN, long long int* RECN, char* FMANE);           // reads Fortran data from the F10 binary file


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      STATUS                                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void status(ofstream* outf){
  char str1[1024];
  
  size_t len = 1024;
  time_t curr_time;
  time(&curr_time);
  
  gethostname(str1, len);
  
  *outf << "Host: " << str1 << ", Date: " << ctime(&curr_time);
  outf->flush();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      REM COM                                                  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int rem_com(char* filename, char* streamstring, int string_length){
  const char com_B = '#';
  const char com_E = '\n';
  
  int pos = 0;
  char cc;

  ifstream inf(filename);
  
  while(inf.get(cc)&& pos < string_length-1){
    if(cc != com_B) 
      streamstring[pos++] = cc;
    else{
      while(cc != com_E && inf.get(cc));
      streamstring[pos++] = com_E;
    }
  }
  streamstring[pos] = 0;
  if(pos == string_length-1){
    cerr << "Buffer size exceeded !\n"; exit(0);
  }
  return(strlen(streamstring));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                get_sys_size(char* sysfile, int* nroao, int* nroa,             */
/*                              long long int* nrofint)                          */
/*                                                                               */
/* Read the system sizes form sysfile.                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint){
  
  ifstream inf(sysfile);
  
  inf.read((char *) nroao,   sizeof(int));
  inf.read((char *) nroa,    sizeof(int));
  inf.read((char *) nrofint, sizeof(long long int));
  
  inf.close();
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys(........)                                             */
/*                                                                               */
/*                                                                               */
/* Read the system data                                                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
	      double *Dz, long long int* sortcount, double* intval, 
	      unsigned short* intnums){
  ifstream datf(sysfile);
  
  int nroao, nroa;
  long long int nrofint;

  //SYSTEM DATA
  datf.read((char *) &nroao , sizeof(int));
  datf.read((char *) &nroa  , sizeof(int));
  datf.read((char *) &nrofint,  sizeof(long long int));
  datf.read((char *) coord  , sizeof(double)*3*nroa);
  datf.read((char *) charges, sizeof(double)*nroa);
  datf.read((char *) mass,    sizeof(double)*nroa);
  
  //ONEL EL INTEGRAL DATA
  datf.read((char * ) Hmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Tmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Smat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dx    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dy    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dz    , sizeof(double)*nroao*nroao);
  
  //TWO EL INTEGRAL DATA
  datf.read((char *) sortcount, sizeof(long long int)*4);
  datf.read((char *) intval,    sizeof(double)*nrofint);
  datf.read((char *) intnums,   sizeof(unsigned short)*nrofint*4);

  datf.close();
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys_1el(........)                                         */
/*                                                                               */
/*                                                                               */
/* Read the system data, only 1 electron integrals                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
	      double *Dz){
  ifstream datf(sysfile);
  
  int nroao, nroa;
  long long int nrofint;

  //SYSTEM DATA
  datf.read((char *) &nroao , sizeof(int));
  datf.read((char *) &nroa  , sizeof(int));
  datf.read((char *) &nrofint,  sizeof(long long int));
  datf.read((char *) coord  , sizeof(double)*3*nroa);
  datf.read((char *) charges, sizeof(double)*nroa);
  datf.read((char *) mass,    sizeof(double)*nroa);
  
  //ONEL EL INTEGRAL DATA
  datf.read((char * ) Hmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Tmat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Smat  , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dx    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dy    , sizeof(double)*nroao*nroao);
  datf.read((char * ) Dz    , sizeof(double)*nroao*nroao);
  
  datf.close();
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               write_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* write  HF wave_function                                                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void write_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs){
  ofstream outf(wavfile);
  
  outf.write((char *) &nroao, sizeof(int));
  outf.write((char *) MOens,  sizeof(double)*nroao);
  outf.write((char *) MOs,    sizeof(double)*nroao*nroao);

  outf.close();
  
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* Read  HF wave_function                                                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs){
  ifstream inf(wavfile);
  int real_nroao = 0;
  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
    cerr << "Wrong HF wavefunction size in read_wav_HF!\n";
    exit(3);
  }
  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                output_matrix(...)                                             */
/*                                                                               */
/*                                                                               */
/* Output lower triangular of symmetric matrix                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/*void   output_matrix(double* mat,int nr_of_col, int nrop, ofstream *outf){
  char dumchar[64];
  int nr_of_blocks = (int) ceil((double)nrop/(double)nr_of_col);
  for(int l = 0; l <  nr_of_blocks; l++){
    int upper_bound = (l+1)*nr_of_col;
    if(upper_bound > nrop) upper_bound = nrop;
    for(int x = l*nr_of_col; x < upper_bound; x++)
      *outf << "\t\t" << x+1 ;
    *outf << "\n\n";
    //loops over centers
    for(int y = l*nr_of_col; y < nrop; y++){
      *outf << y+1 << "\t" ;
      for(int x = l*nr_of_col; x < upper_bound && x <= y; x++){
      sprintf(dumchar,"%+.6e",mat[x*nrop+y]);
           *outf << dumchar << "\t";
      }
      *outf << "\n";
    }
    *outf << "\n";
  }
  outf->flush();
}*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* read_matrix                                                                   */
/*                                                                               */
/* reads CI matrix and energies from a binary matrix file                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void read_matrix(char* matfile, int nros, double* cisvecs, int type){
 ifstream inf(matfile);
 if(type == 1){
  char tmp[128];
  inf.read(tmp,4);
 }
  int dimension = (nros - 1)*(nros - 1);
  double* tmpmat = new double[dimension];
  inf.read((char *) tmpmat, sizeof(double)*dimension);
 inf.close();
 for(int i = 0; i < nros; i++){
  for(int j = 0; j < nros; j++){
      if(i == 0 && j == 0) cisvecs[i*nros+j] = 1.;
      if(i == 0 && j != 0) cisvecs[i*nros+j] = 0.;
      if(i != 0 && j == 0) cisvecs[i*nros+j] = 0.;
      if(i != 0 && j != 0){
	  cisvecs[i*nros+j] = tmpmat[(i-1)*(nros-1)+(j-1)];
      }
   }
 }
}*/


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* read_bin                                                                      */
/*                                                                               */
/* reads RHF energy and CI excitation energies from a binary file                */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*void read_bin(char* binfile, int nros, double& cisEgs, double* cisvals){
 ifstream inf(binfile);
 int real_nros = 0;
 inf.read((char *) &real_nros, sizeof(int));
 if(real_nros != nros){
  cerr << "Number of states in input and binfile not equivalent!\n";
  exit(51);
 }
 inf.read((char *) &cisEgs, sizeof(double));
 inf.read((char *) cisvals, sizeof(double)*nros);
 inf.close();
}*/
