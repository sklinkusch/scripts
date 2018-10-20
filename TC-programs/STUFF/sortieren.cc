#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 3){ 
  cerr << "Need input- and output file!\n";
  exit(1);
 }

  ifstream inf;
  ofstream outf;
  inf.open(argv[1]);
  int nrol;
  inf >> nrol;
  double* pops = new double[nrol];
  int* state = new int[nrol];
  double* energ = new double[nrol];
  double* tempval = new double[2];
  int tmpval;
  for(int x = 0; x < nrol; x++){
   inf >> pops[x] >> state[x] >> energ[x];
  }
  inf.close();
  int count;
  //Sortieralgorithmus
  do{
   count = 0;
   for(int x = 0; x < (nrol-1); x++){
    if(pops[x+1] > pops[x]){
     tempval[0] = pops[x+1];
     tmpval = state[x+1];
     tempval[1] = energ[x+1];
     pops[x+1] = pops[x];
     state[x+1] = state[x];
     energ[x+1] = energ[x];
     pops[x] = tempval[0];
     state[x] = tmpval;
     energ[x] = tempval[1];
     count++;
    }
   }
  }while(count != 0);

 outf.open(argv[2]);
 for(int x = 0; x < nrol; x++){
  outf << pops[x] << " " << state[x] << " " << energ[x] << "\n";
 }
 outf.close();
}

