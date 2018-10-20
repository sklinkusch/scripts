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
 if(argc != 3){
  cerr << "Need input- and output-file!\n";
  exit(1);
 }

 ifstream inf;
 ofstream outf;

 inf.open(argv[1]);
 int nros;                                  // Number of states
 inf >> nros;
 double* energies = new double[nros];       // Array with energies
 double* rates    = new double[nros];       // Array with ionization rates
 double* dipoles  = new double[nros];       // Array with transition dipole moments
 double* dipolx   = new double[nros];       // Transition dipoles along x
 double* dipoly   = new double[nros];       // Transition dipoles along y
 double* dipolz   = new double[nros];       // Transition dipoles along z

 for(int x = 0; x < nros; x++){
  inf >> energies[x] >> rates[x] >> dipolx[x] >> dipoly[x] >> dipolz[x];
 }

 for(int x = 0; x < nros; x++){
  dipoles[x] = sqrt(pow(dipolx[x],2) + pow(dipoly[x],2) + pow(dipolz[x],2));
 }

 inf.close();
 outf.open(argv[2]);
 double width_old = energies[nros-1] - energies[0];
 double width_new = width_old + 0.5;

 int points = 1000.*nros;
 double dE = width_new/((double) points-1.);

 double* dos = new double[points];
 double sum1, sum2;
 
 for(int x = 0; x <= points; x++){
  dos[x] = 0.;
  for(int y = 0; y < nros; y++){
   sum1 = (x*dE)-energies[y];
   sum2 = rates[y]/2.;
   dos[x] += fabs(dipoles[y])*rates[y]/(2.*M_PI*pow(sum1,2.)+pow(sum2,2.));
  }
  outf << x*dE << "\t" << dos[x] << "\n";
 }
 outf.close();
}

