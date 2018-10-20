/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: twostate.cc                                                            *
 *                                                                              *
 * transforms a N-state system to a two-level-system for pure 2-level dynamics  * 
 *                                                                              *
 *                                                     Stefan Klinkusch  2011   *
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
    cerr << "Need <bcs-file> <tls-prefix>\n";
    exit(0);
  }

  char dumc[1024];
  int nroao, nroe, llim, ulim;
  ifstream datf(argv[1]);
  sprintf(dumc, "%s.tls", argv[2]);
  ofstream tlsf(argv[2]);

// read and write unchanged data
  
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  tlsf.write((char*) &nroao, sizeof(int));
  tlsf.write((char*) &nroe,  sizeof(int));
  tlsf.write((char*) &llim,  sizeof(int));
  tlsf.write((char*) &ulim,  sizeof(int));
  
  int omo                = nroe/2 - llim;
  int umo                = ulim - nroe/2+1;
  int cis_size           = omo*umo+1;
  int nros               = cis_size;

  cout << "Nr of CIS states is " << cis_size << "\n";
  
  double* ens   = new double[cis_size];
  double* vecs  = new double[cis_size*cis_size];
  double* dmat  = new double[cis_size*cis_size];
  double* dvals = new double[cis_size];
  double* dummat= new double[cis_size*cis_size];
  double* corr  = new double[cis_size];
  
  datf.read((char *) ens, cis_size * sizeof(double));
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));
  
// Correction of Energies (from the 2nd excited state, + 10 hartrees)
    
  for(int x = 0; x < cis_size ; x++){
     if(x > 1){
      ens[x] += 10.;
     }
  }
  tlsf.write((char *) ens, cis_size * sizeof(double));
  tlsf.write((char *) vecs, cis_size* cis_size * sizeof(double));
  
//  Permanent dipoles in x
  
  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

  for(int x = 0; x<  nros; x++){
    for(int y = 0; y < cis_size; y++){
      if(x > 1){
       dmat[y*cis_size+x] = 0.;
      }
    }
  }
// Transition dipole moments in x
  for(int x = 0; x<  nros-1; x++){
      for(int y = x+1; y < nros; y++){
	for(int z = 0; z < cis_size; z++){
          if(x > 1){
           dmat[z*cis_size+x] = 0.;
          }
        }
      }
  }
  tlsf.write((char *) dvals, cis_size * sizeof(double));
  tlsf.write((char *) dmat, cis_size* cis_size * sizeof(double));

// Permanent dipoles in y

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];
  
  for(int x = 0; x<  nros; x++){
    for(int y = 0; y < cis_size; y++){
     if(x > 1){
      dmat[y*cis_size+x] = 0.;
     }
    }
  }
// Transition dipole moments in y
  for(int x = 0; x<  nros-1; x++){
      for(int y = x+1; y < nros; y++){
	for(int z = 0; z < cis_size; z++){
         if(x > 1){
          dmat[z*cis_size+x] = 0.;
         }
        }
      }
  }
  tlsf.write((char *) dvals, cis_size * sizeof(double));
  tlsf.write((char *) dmat, cis_size* cis_size * sizeof(double));

// Permanent dipole moments in z

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

  for(int x = 0; x<  nros; x++){
    for(int y = 0; y < cis_size; y++){
      if(x > 1){
       dmat[y*cis_size+x] = 0.;
      }
    }
  }
  for(int x = 0; x<  nros-1; x++){
// Transition dipole moments in z
      for(int y = x+1; y < nros; y++){
	for(int z = 0; z < cis_size; z++){
          if(x > 1){
           dmat[z*cis_size+x] = 0.;
          }
        }
      }
  }

  tlsf.write((char *) dvals, cis_size * sizeof(double));
  tlsf.write((char *) dmat, cis_size* cis_size * sizeof(double));

  datf.read((char *) corr, cis_size * sizeof(double));
  tlsf.write((char *) corr, cis_size * sizeof(double));
  
  datf.close();
  tlsf.close();
  

 }

