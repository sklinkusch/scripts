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

/************************************************************************************
 * readWAV                                                                          *
 * developped by Stefan Klinkusch, Free University of Berlin (2014)                 *
 * reads the MO energies and vectors from a binary .dat.hwf file and writes the     *
 * data to a human readable log file                                                *
 ************************************************************************************/

extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);

int main(int argc, char* argv[]){
 if(argc != 4){
  cerr << "Usage: ./readwav wavfile #MOs logfile\n";
  exit(1);
 }

 int nroao = atoi(argv[2]);
 double* MOens = new double[nroao];
 double* MOs = new double[nroao*nroao];
 char wavfile[256];
 sprintf(wavfile, "%s", argv[1]);
 clog << "WAVfile is " << wavfile << "\n";
 read_wav_HF(wavfile, nroao, MOens, MOs);
 clog << "WAVfile read\n";
 ofstream outf;
 outf.precision(12);
 outf.open(argv[3]);
 clog << "Outfile open\n";
 outf << "Number of MOs: " << nroao << "\n";
 outf << "-------------------------------------------\n";
 outf << "Eigenenergies: \n";
 for(int x = 0; x < nroao; x++){
    outf << x << "  " << MOens[x] << "\n";
 }
 cout << "Done with energy part\n";
 outf << "-------------------------------------------\n";
 outf << "Eigenvectors: \n";
 for(int x = 0; x < nroao; x++){
  outf << "VECTOR " << x << " : \n";
  for(int y = 0; y < nroao; y++){
   outf << MOs[x*nroao+y] << "  ";
   if(y%5 == 4) outf << "\n";
  }
  outf << "\n";
 }
 clog << "Done with vector part\n";
 outf.close();
}

