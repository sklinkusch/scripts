#include<fstream>
#include<iostream>
#include<math.h>

using namespace std;


int main(int argc, char* argv[]){
 if(argc != 3){
 cerr << "Need input-file output-prefix!\n";
 exit(1);
 }
  ifstream inf(argv[1]);
  int nroao;


  inf.read((char *) &nroao, sizeof(int));

  double* MOens = new double[nroao];
  double* MOs = new double[nroao*nroao];

  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();

  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);

  outf << "Number of atomic orbitals: \t" << nroao << "\n";
  outf << "MO energies:\n";
  for(int x = 0; x < nroao; x++)
  outf << "MO No. " << x << ": " << MOens[x] << "\n";
  outf << "#-------------------------------------------------------#\n";
  outf << "MO coefficients:\n";
  for(int x = 0; x < nroao; x++){
   for(int y = 0; y < nroao; y++){
    outf << "(" << x << "," << y << ")\t" << MOs[y*nroao+x] << "\n";
   }
  }
}


