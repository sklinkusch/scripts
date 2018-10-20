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
  cout << "Nr of atomic orbitals: " << nroao << "\n";
  datf.read((char*) &nroe,  sizeof(int));
  cout << "Nr of electrons: " << nroe << "\n";
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  cout << "Limits for correlation: " << llim << " " << ulim << "\n";
  
  //  int nrof               = ulim - llim +1;
  int omo                = nroe/2 - llim;
  int umo                = ulim - nroe/2+1;
  long long int cis_size = omo*umo+1;
  cout << "MOs (occ./virt.), states: " << omo << " " << umo << " " << cis_size << "\n";
  
  double* ens   = new double[cis_size];
  double* vecs  = new double[cis_size*cis_size];
  double* dmat  = new double[cis_size*cis_size];
  double* dvals = new double[cis_size];
  double* dummat= new double[cis_size*cis_size];
  double* corr  = new double[cis_size];
  
  datf.read((char *) ens, cis_size * sizeof(double));
  cout << "Energies: " << "\n";
  for(int i = 0; i < cis_size; i++){
      cout << ens[i] << "\n";
  }
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));
  cout << "Vectors: " << "\n";
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++){
	  cout << vecs[i*cis_size+j] << "\n";
      }
  }
  datf.read((char *) dvals, cis_size * sizeof(double));
  cout << "Dx: " << "\n";
  for(int i = 0; i < cis_size; i++){
      cout << dvals[i] << "\n";
  }
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  cout << "Dxvecs: " << "\n";
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++){
	  cout << dmat[i*cis_size+j] << "\n";
      }
  }
  datf.read((char *) dvals, cis_size * sizeof(double));
  cout << "Dy: " << "\n";
  for(int i = 0; i < cis_size; i++){
      cout << dvals[i] << "\n";
  }
  
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  cout << "Dyvecs: " << "\n";
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++){
	  cout << dmat[i*cis_size+j] << "\n";
      }
  }
  datf.read((char *) dvals, cis_size * sizeof(double));
  cout << "Dz: " << "\n";
  for(int i = 0; i < cis_size; i++){
      cout << dvals[i] << "\n";
  }
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  cout << "Dzvecs: " << "\n";
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++){
	  cout << dmat[i*cis_size+j] << "\n";
      }
  }
  datf.read((char *) corr, cis_size * sizeof(double));
  cout << "Corr: " << "\n";
  for(int i = 0; i < cis_size; i++){
      cout << corr[i] << "\n";
  }
  datf.close();

 }
