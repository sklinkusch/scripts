#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#define Complex complex<double>

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 6){
  cerr << "Need #timesteps, #CSFs, timestep for output, wav-file (input), iwav-file (output)!\n";
  exit(1);
 }

int nrots = atoi(argv[1]);
int cis_size = atoi(argv[2]);
int timestep = atoi(argv[3]);
double curr_time;
Complex PsiE[cis_size];

ifstream wavf(argv[4]);
ofstream iwav(argv[5]);

for(int x = 0; x < nrots; x++){
 wavf.read((char *) &curr_time, sizeof(double));
 wavf.read((char *) PsiE, sizeof(Complex)*cis_size);
 if(x == timestep){
  cout << "Timestep for output: " << curr_time << "\n";
  curr_time = 0.;
  iwav.write((char *) &curr_time, sizeof(double));
  iwav.write((char *) PsiE, sizeof(Complex)*cis_size);
  break;
 }
}
}


