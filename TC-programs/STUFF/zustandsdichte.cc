#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

//Programm zur Überprüfung, ob ein Kontinuum vorliegt
//Stefan Klinkusch, 2009

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 6){
  cerr << "Usage ./zustandsdichte <nros> <input-file> < 1 if eV, 0 if E_h> <constant width in previous unit> output-file!\n";
  exit(1);
 }

 ifstream inf;
 ofstream outf;

 int nros = atoi(argv[1]);                                  // Number of statesi
 cout << "Nr of states: " << nros << "\n";
 inf.open(argv[2]);
 cout << "Reading from file: " << argv[2] << "\n";
 double* energies = new double[nros];       // Array with energies
 int modus = atoi(argv[3]);
 if(modus == 1){
     cout << "Width is given in eV, converting to E_h\n";
 }else{
     cout << "Width is given in E_h\n";
 }
 double width = strtod(argv[4], NULL);
 cout << "width: " << width << "\n";
 if(modus == 1){
     width /= 27.211385;
     cout << "converted width: " << width << "\n";
 }
 double* rates    = new double[nros];       // Array with ionization rates
// double* dipolx   = new double[nros];       // Transition dipoles along x
// double* dipoly   = new double[nros];       // Transition dipoles along y
// double* dipolz   = new double[nros];       // Transition dipoles along z
 cout << "Reading energies and populations \n";
 for(int x = 0; x < nros; x++){
  inf >> energies[x] >> rates[x];
 }
 inf.close();
 outf.open(argv[5]);
 cout << "Spectrum is written out to " << argv[5] << "\n";
 outf.precision(12);
 double width_old = energies[nros-1];
 double width_new = width_old + 0.5;
 cout << "maximal frequency : " << width_new << "\n";
// const double ip  = 0.39079;

 double dE = 1e-04; 
 cout << "steplength: " << dE << "\n";
 long long int points = (long long int) (width_new / (dE)) + (long long int) 1;
 cout << "Nr of points: " << points << "\n";

 double sum1, sum2;
 double dosix = 0.;
 
 for(long long int x = 0; x <= points; x++){
  dosix = 0.;
#pragma omp parallel for reduction(+:dosix)
  for(int y = 0; y < nros; y++){
   sum1 = (x*dE)-energies[y];
   sum2 = width/2.;
   dosix += rates[y]*width/(2.*M_PI*pow(sum1,2.)+pow(sum2,2.));
//   if(fabs(sum1) < 5.e-5 && energies[y] < ip){
//    dos[x] = 500.;
//   }
  }
  outf << x*dE << "\t" << dosix << "\n";
 }
 outf.close();
}

