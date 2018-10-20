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

/****************************************************************************************
 * rgw (Read GAMESS Wavefunction)                                                       *
 * developped by Tillmann Klamroth, University of Potsdam (2004)                        *
 * modified by Stefan Klinkusch, Free University of Berlin (2014)                       *
 * reads MO energies and MO vectors from a modified GAMESS punch file and writes it to  *
 * a binary .dat.hwf file                                                               *
 ****************************************************************************************/

void read_line(ifstream* inf, char* line_buf, int max_count){                   // routine to read strings from a file (a line is saved into line_buf)
  int count = 0;
  char c;
  int end = 0;
  while(end == 0 && inf->get(c)){
    line_buf[count++] = c;
    if(c == '\n') end = 1;
    if(count == max_count){
      cerr << "Buffer overflow while parsing\n"; exit(1);
    }
  }
  line_buf[count-1] = 0;
}

int main(int argc, char* argv[]){
  if(argc != 4){
    cerr << "Need <nrods> <GAMESS log-file> <outfile>\n";
    exit(0);
  }
  int nrods = atoi(argv[1]);                                                   // number of AOs
  double* cisens = new double[nrods];                                          // array for CIS excitation energies
  double* oscstr = new double[nrods];                                          // array for CIS oscillator strengths
  char* state = new char[2];
  double hart, evolt, kcal, cm, nm;
  char* oscillator = new char[10];
  char* strength = new char[8];
  char equal;
  double osc;
 
  ifstream inf(argv[2]);                                                       // input stream from file
  clog << "Reading from " << argv[2] << " Nrods is " << nrods << "\n";
  char line_buf[64000];                                                        // string where a single line from the GAMESS punch file is saved
  line_buf[1] = 0;                                                             // is first set to zero
  clog << "Reading CIS energies\n";
  while(strcmp("                  CI-SINGLES EXCITATION ENERGIES", line_buf)!=0) read_line(&inf, line_buf,64000); // read until line with " $EIGENVALUES" is reached
  for(int dumint = 0; dumint < 2; dumint++) read_line(&inf, line_buf, 64000);
    for(int i = 0; i < nrods; i++){                                        // loop over all lines with MO energies
     read_line(&inf, line_buf,64000);                                          // read each single line and write it to line_buf
     istringstream ist(line_buf);                                              // input stream from the buffer line_buf, able to read doubles
     ist >> state >> hart >> evolt >> kcal >> cm >> nm; 
     cisens[i] = hart;
     clog << hart << "\n";
     }
   inf.close();
   ifstream iof(argv[2]); 
  clog << "Reading oscillator strengths\n";
  while(strcmp("                   EXPECTATION VALUES OF DIPOLE MOMENTS",line_buf)!=0)   read_line(&inf, line_buf,64000);     // read lines until " $VEC   " is reached
  for(int dumint = 0; dumint < 12; dumint++) read_line(&iof, line_buf, 64000);
  for(int i = 0; i < nrods; i++){
      read_line(&iof, line_buf, 64000);
      istringstream ist(line_buf);
      ist >> oscillator >> strength >> equal >> osc;
      oscstr[i] = osc;
      for(int dumint = 0; dumint < 10; dumint++) read_line(&iof, line_buf, 64000);
  }
  ofstream outf;
  outf.open(argv[3]);
  for(int i = 0; i < nrods; i++){
  outf << cisens[i] << "    " << oscstr[i] << "\n";
      outf.flush();
  }
  outf.close();
}

