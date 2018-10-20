#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 5 && argc != 6){
  cerr << "Need out-index, Number of timesteps, CSFs, Input-File, (Output-File)!\n";
  exit(1);
 }
 int outdex = atoi(argv[1]);
 int dt = atoi(argv[2]);    // Number of time steps
 int N = atoi(argv[3]);    // Number of CSFs
 double* norm = new double[dt];  //Allocates an array
 long double pop;
 double* time = new double[dt];   
 double curr;
 double energ;
 ifstream inf;
 inf.open(argv[4]);
 
for(int i = 0; i < (dt); i++){
 time[i] = 0.;
 norm[i] = 0.;
 for(int x = 0; x < N; x++){ 
 inf >> curr >> energ >> pop;
 time[i] = curr;
 norm[i] += pop;
 }
}
 inf.close();
if(outdex == 1){
 if(argc != 6){
  cerr << "Need output-file!\n";
  exit(1);
 }
ofstream outf;
outf.precision(12);
outf.open(argv[5]);
for(int i = 0; i < (dt); i++){
 outf << time[i] << "\t" << norm[i] << "\n";
 }
 outf.close();
}else{
for(int i = 0; i < (dt); i++){
 cout << time[i] << "\t" << norm[i] << "\n";
 }
}
}


