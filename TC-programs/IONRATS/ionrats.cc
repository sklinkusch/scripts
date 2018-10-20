/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cis.cc                                                                 *
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
#include <algorithm>

using namespace std;

// Functions
double calc_ionrat(double en, double ip, double damp, int state, int nros, int* frommo, int* tomo, double* mocoeff, double* MOens);
void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);


int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  char datfile[256];      // file containing MO transitions
  char ensfile[256];      // file containing CIS and CIS(D) energies
  char wavfile[256];      // file containing MOs and MO energies
  int   calc_d=0;         // if 1 -> calculate ionization rates from CIS(D) energies, else -> calculate ionization rates from CIS energies
  double damp = 1.;

  if(argc != 9){
    cerr << "Usage: ./ionrats <calc_d> <nroao> <nros> <ip> <ensfile> <datfile> <wavfile> <output-prefix>\n";
    exit(1);
  }
  
  //SCF VARIABLES 

  //system size
  calc_d = atoi(argv[1]);
  if(calc_d != 0 && calc_d != 1) calc_d = 0;
  int    nroao = atoi(argv[2]);                 //Nr of basis functions
  if(nroao <= 0){
      cerr << "Number of atomic orbitals must be positive\n";
      exit(2);
  }
  int nros = atoi(argv[3]);
  if(nros <= 0){
      cerr << "Number of excited electronic states must be positive\n";
      exit(3);
  }
  double ip = strtod(argv[4],NULL);
  sprintf(ensfile, "%s", argv[5]);
  sprintf(datfile, "%s", argv[6]);
  sprintf(wavfile, "%s", argv[7]);
  int cissize = nros+1;
  ifstream ensf;
  int state;
  double cisenergy;
  double* cisvals = new double[cissize];
  double* cisd_corr = new double[cissize];
  cisvals[0] = 0.;
  cisd_corr[0] = 0.;
  ensf.open(ensfile);

  for(int i = 1; i < cissize; i++){
      ensf >> state >> cisenergy >> cisvals[i] >> cisd_corr[i];
  }
  ensf.close();
  double* stateens = new double[cissize];
  if(calc_d == 1){
      for(int x = 0; x < cissize; x++){
	  stateens[x] = cisvals[x] + cisd_corr[x];
      }
  }else{
      for(int x = 0; x < cissize; x++){
	  stateens[x] = cisvals[x];
      }
  }
  double* MOs = new double[nroao*nroao];
  double* MOens = new double[nroao];
  read_wav_HF(wavfile, nroao, MOens, MOs);
  delete [] MOs;

  int* frommo = new int[nros*nros];
  int* tomo = new int[nros*nros];
  double* mocoeff = new double[nros*nros];
  ifstream datf;
  datf.open(datfile);
  for(int s = 0; s < nros; s++){
      for(int motr = 0; motr < nros; motr++){
	  datf >> state >> frommo[s*nros+motr] >> tomo[s*nros+motr] >> mocoeff[s*nros+motr];
	  frommo[s*nros+motr] -= 1;
	  tomo[s*nros+motr] -= 1;
      }
  }
  
   //##CALCULATION OF IONIZATION RATES / WRITE OUT TO IRX-FILE

  double* ionrat = new double[cissize];
  double* lifetime = new double[cissize];
  ionrat[0] = 0.;
  lifetime[0] = -1.;
  for(int x = 1; x < cissize; x++){
      ionrat[x] = calc_ionrat(stateens[x], ip, damp, x, nros, frommo, tomo, mocoeff, MOens);
      if(ionrat[x] == 0.){
	  lifetime[x] = -1.;
      }else{
	  lifetime[x] = 1./ionrat[x];
      }
  }

  char dumc[256];
  sprintf(dumc,"%s.irx", argv[8]);
  ofstream irxf;
  irxf.open(dumc);
  char dumline[4096];
  for(int x = 0; x < cissize; x++){
      sprintf(dumline, "%4d   %12.8f   %12.8f", x, ionrat[x], lifetime[x]);
      irxf << dumline<< "\n";
      irxf.flush();
  }
  irxf.close();
}

  
// Functions

double calc_ionrat(double en, double ip, double damp, int state, int nros, int* frommo, int* tomo, double* mocoeff, double* MOens){
    double gamma;
    int MO_r;
    if((en - ip) < 0.){
	gamma = 0.;
    }else{
	gamma = 0.;
	for(int x = 0; x < nros; x++){
	    MO_r = tomo[(state-1)*nros+x];
	    if(MOens[MO_r] > 0.){
		gamma += pow(fabs(mocoeff[(state-1)*nros+x]),2)*damp*sqrt(MOens[MO_r]);
	    }else{
		gamma += 0.;
	    }
	}
    }
    return(gamma);
}

void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs){
  ifstream inf(wavfile);
  int real_nroao;

  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
    cerr << "Wrong HF wavefunction size in read_wav_HF!\n"; exit(3);
  }

  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();
}

