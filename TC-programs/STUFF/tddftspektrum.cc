/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cispektrum.cc                                                          *
 *                                                                              *
 *                                                                              * 
 *                                                                              *
 *                                                    Stefan Klinkusch   2011   *
 ********************************************************************************/ 
#include <fstream>
#include <iostream> 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h> 
#include <sys/resource.h>

using namespace std;

// Functions
void status(ofstream* outf);
int rem_com(char* filename, char* streamstring, int string_length);

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need <input-file> <output-prefix>\n";
    exit(0);
  }

    
  char dumc[1024];
  sprintf(dumc, "%s.log", argv[2]);
  ofstream outf(dumc);
  status(&outf);
  outf << ".....................................................................................\n";
  outf << "Reading input from: " << argv[1] << "\n";
  int buff_length = 65536;
  char* file_buff = new char[buff_length];
  rem_com(argv[1], file_buff, buff_length);
  istringstream ist(file_buff);
  char eaffile[1024];
  int nros;
  double sigma, extra, incr, l_min;
  ist >> eaffile >> nros >> sigma >> l_min >> extra >> incr;
  outf << "File containing excitation energies and oscillator strengths: " << eaffile << "\n";
  outf << "Number of electronic states: " << nros << "\n";
  outf << "Smooth factor for Gaussian functions: " << sigma << "\n";
  outf << "Minimum wavelength: " << l_min << "\n";
  outf << "extra wavelength: " << extra << "\n";
  outf << "Increment: " << incr;
  outf.flush();
  
  double* cens  = new double[nros];
  double* ens   = new double[nros];
  double* exv   = new double[nros];
  double* ev    = new double[nros];
  double* ozstr = new double[nros];
  ifstream eaff(eaffile);
  for(int i = 0; i < nros; i++){
      eaff >> ens[i] >> ozstr[i];
  }
#pragma omp parallel for
  for(int x = 0; x < nros; x++){
   cens[x] = ens[x];                                    // no need to correct energies for TD-DFT
   exv[x] = 219474.63*cens[x];                          // conversion: energy in hartree -> wavelength in cm
   ev[x] = 27.211385*cens[x];                           // conversion: hartree -> eV
  }

  for(int x = 0; x < nros; x++){
   outf << x+1 << " " << cens[x] << " " << ozstr[x] << "\n";
   outf.flush();
  }

  // Constants
  double h = 2*M_PI;
  double c0 = 137.036;
  double l_max = 0.052917721*h*c0/cens[1] + extra;      // maximal wavelength (150 nm longer than the peak for the first excited state)
  int nrop = (int) ((l_max-l_min)/incr) + 1;            // number of points in the spectrum file
  double epsilon;                                       // extinction in L/(mol cm)
  double cl, cv;                                        // current wavelength in nm and cm
  double ch, cev;                                       // current energy in hartree and eV

  double dl = incr;
  double kappa = 4.3189984e-10;
  double nfac = sigma*sqrt(2*M_PI);
  double prefac = 0.1/(kappa*nfac);

  char specfile[1024];
  sprintf(specfile, "%s.spec", argv[2]);
  ofstream spcf;
  outf << "Writing data for spectrum to " << specfile << "\n";
  outf.flush();
  spcf.open(specfile);
  for(int i = 0; i < nrop; i++){
   epsilon = 0.;
   cl = l_min + i*dl;                                    // current wavelength in nm
   cv = 1./cl*10000000.;                                 // current wavelength in cm
   ch = cv/219474.63;                                    // current energy in hartree
   cev = 27.211385*ch;                                   // current energy in eV
#pragma omp parallel for reduction(+:epsilon)
   for(int s = 1; s < nros; s++){
    epsilon += ozstr[s]*prefac*exp(-0.5*pow((cv-exv[s])/sigma,2));
   }
   spcf << cl << " " << ch << " " << cev << " " << epsilon << "\n"; 
   spcf.flush();
  }
  spcf.close();
  status(&outf);
  outf.close();
 }

void status(ofstream* outf){
  char str1[1024];
  
  size_t len = 1024;
  time_t curr_time;
  time(&curr_time);
  
  gethostname(str1, len);
  
  *outf << "Host: " << str1 << ", Date: " << ctime(&curr_time);
  outf->flush();
}

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

