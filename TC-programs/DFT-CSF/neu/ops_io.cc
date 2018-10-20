/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_io.cc                                                              *
 *                                                                              *
 * contains io routines                                                         * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
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
void status(ofstream* outf);
int  rem_com(char* filename, char* streamstring, int string_length);
void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
	      double *Dz); 
void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
		  double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		  double *Dz);
void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);

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

void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint){
  
  ifstream inf(sysfile);
  
  inf.read((char *) nroao,   sizeof(int));
  inf.read((char *) nroa,    sizeof(int));
  inf.read((char *) nrofint, sizeof(long long int));
  
  inf.close();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys(........)                                             */
/*                                                                               */
/*                                                                               */
/* Read the system data                                                          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
	      double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
	      double *Dz){ //, long long int* sortcount, double* intval, 
//	      unsigned short* intnums){
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
//  datf.read((char *) sortcount, sizeof(long long int)*4);
//  datf.read((char *) intval,    sizeof(double)*nrofint);
//  datf.read((char *) intnums,   sizeof(unsigned short)*nrofint*4);

  datf.close();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_sys_1el(........)                                         */
/*                                                                               */
/*                                                                               */
/* Read the system data, only 1 electron integrals                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
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
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                read_wav_HF(........)                                          */
/*                                                                               */
/*                                                                               */
/* Read  HF wave_function                                                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs){
  ifstream inf(wavfile);
  int real_nroao;

  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao < nroao){
    cerr << "Wrong HF wavefunction size in read_wav_HF!\n"; exit(3);
  }

  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();
}

