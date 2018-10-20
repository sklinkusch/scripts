#include<fstream>
#include<iostream>
#include<stdlib.h>

/* Program for separation of degenerated states' properties
Stefan Klinkusch, 2008 */

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 5 && argc != 6 && argc != 7){
  cerr << "Need in-index, out-index, #timesteps, #degenerated states, no. of this state in a row, (input-file, output-file)!\n";
  exit(1);
  }

ifstream inf;

int index = atoi(argv[1]);
int outdex = atoi(argv[2]);
int nrots = atoi(argv[3]);
int nrods = atoi(argv[4]);
int state = atoi(argv[5]);
int max = nrods*nrots;
double* time = new double[max];
double* energy = new double[max];
double* pop = new double[max];

if(index == 1){
 if(argc != 6 && argc != 7){
  cerr << "Need input-file!\n";
  exit(1);
 }
inf.open(argv[6]);
for(int x = 0; x < max; x++){
 inf >> time[x] >> energy[x] >> pop[x];
}
inf.close();
}else{
 if(argc != 5 && argc != 6){
  cerr << "Too many arguments!\n";
  exit(1);
 }
for(int x = 0; x < max; x++){
 cin >> time[x] >> energy[x] >> pop[x];
}
}

if(outdex == 1){
 if(argc != 6 && argc != 7){
  cerr << "Need output-file!\n";
  exit(1);
  }
ofstream outf;
outf.open(argv[6]);

for(int x = 0; x < max; x++){
 if(x%nrods == state){
  outf << time [x] << " " << energy[x] << " " << pop[x] << "\n";
 }else{
 continue;
 }
}
outf.close();
}else{
 if(argc != 5 && argc != 6){
  cerr << "Too many arguments!\n";
  exit(1);
 }
for(int x = 0; x < max; x++){
 if(x%nrods == state){
  cout << time[x] << " " << energy[x] << " " << pop[x] << "\n";
 }else{
  continue;
  }
 }
 }
}



