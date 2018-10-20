/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: read_bcs.cc                                                            *
 *                                                                              *
 * extrats human readable high precision data from an bcs file                  * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/ 
#include <fstream>
#include <iostream> 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);

using namespace std;

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Nedd <nr of print out states> <bcs-file>\n";
    exit(0);
  }

    
  int nros = atoi(argv[1]);

  cout << "Printing out information for " << nros << " states from " << argv[2] << "\n";

  char dumc[1024];
  int nroao, nroe, llim, ulim;
  ifstream datf(argv[2]);
  
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  
  //  int nrof               = ulim - llim +1;
  int omo                = nroe/2 - llim;
  int umo                = ulim - nroe/2+1;
  long long int cis_size = omo*umo+1;

  cout << "Nr of CIS states is " << cis_size << "\n";
  
  double* ens   = new double[cis_size];
  double* vecs  = new double[cis_size*cis_size];
  double* dmat  = new double[cis_size*cis_size];
  double* dvals = new double[cis_size];
  double* dummat= new double[cis_size*cis_size];
  double* corr  = new double[cis_size];
  
  datf.read((char *) ens, cis_size * sizeof(double));
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));
  
  cout << "===============================================================================\n";
  cout << "Uncorrected energies:\n";
  
  
  for(int x = 0; x<  nros; x++){
    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]);
    cout << dumc << "\n";
  }
 for(int x = 0; x < nros; x++){
  double max_coeff = 0.;
    int mi = 0, mf = 0;
    for(int y = 0; y < cis_size; y++){
      if( fabs(vecs[x*cis_size+y]) > fabs(max_coeff)){
        max_coeff = vecs[x*cis_size+y];
        mi = (y-1)/umo+llim;
        mf = (y-1)%umo+omo+llim;
      }
    }
  sprintf(dumc,"%.8f\t%.8f,\t %.3i->%.3i\t%+.8f\n",ens[x],ens[x]*27.211383,mi,mf,max_coeff/sqrt(2.));
  cout << dumc << "\n";
}
}

