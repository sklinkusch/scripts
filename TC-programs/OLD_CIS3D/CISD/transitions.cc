/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: read_cisd-bcs.cc                                                       *
 *                                                                              *
 * extrats human readable high precision data from an bcs file                  * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 *                                                    modified for CISD PK 2006 *
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
  //  long long int cis_size = omo*umo+1;
  long long int cisd_size;

  // Loops to determine cisd_size
  int sum2 = 0;
  for (int occ = 1; occ <= omo-1; occ++){
    sum2 += occ;
  }
  
  int sum1 = 0;
  for (int virt = 1; virt <= umo-1; virt++){
    sum1 += virt;
  }
  
  cisd_size = 1 + 2*omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;

  cout << "Nr of CISD states is " << cisd_size << "\n";
  
  double* ens   = new double[cisd_size];
  double* vecs  = new double[cisd_size*cisd_size];
  double* dmat  = new double[cisd_size*cisd_size];
  double* dvals = new double[cisd_size];
  double* dummat= new double[cisd_size*cisd_size];
  double* corr  = new double[cisd_size];
  
  datf.read((char *) ens, cisd_size * sizeof(double));
  datf.read((char *) vecs, cisd_size* cisd_size * sizeof(double));

for(int x = 1; x < nros; x++){
    double max_coeff = 0.;
    int mi = 0, mf = 0;
    for(int y = 0; y < cisd_size; y++){
      if( fabs(vecs[x*cisd_size+y]) > fabs(max_coeff)){
        max_coeff = vecs[x*cisd_size+y];
        mi = (y-1)/umo+llim;
        mf = (y-1)%umo+omo+llim;
      }
    }
    sprintf(dumc,"%.8f\t%.8f,\t %.3i->%.3i\t%+.8f\n",(ens[x]-ens[0]),(ens[x]-ens[0])*27.211383,mi,mf,max_coeff/sqrt(2.));
    cout << x << "          "<< dumc;
  }  

datf.close();
  
  

 }
