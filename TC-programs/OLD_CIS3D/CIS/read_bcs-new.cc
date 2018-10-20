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
  int nroao, nroe, llim, ulim, nstates;
  ifstream datf(argv[2]);
  
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  datf.read((char*) &nstates, sizeof(int));
  
  //  int nrof               = ulim - llim +1;
  int omo                = nroe/2 - llim;
  int umo                = ulim - nroe/2+1;
  long long int cis_size = omo*umo+1;

  cout << "Nr of CIS states is " << cis_size << "\n";
  
  double* ens   = new double[nstates];
  double* vecs  = new double[nstates*cis_size];
  double* dmat  = new double[nstates*nstates];
  double* dvals = new double[nstates];
  double* dummat= new double[nstates*nstates];
  double* corr  = new double[cis_size];
  cout << "Nr of included CIS states is " << nstates << "\n";
  
  datf.read((char *) ens, nstates* sizeof(double));
  cout << "Uncorrected energies read\n";
  datf.read((char *) vecs, nstates* cis_size * sizeof(double));
  cout << "Eigenvectors read\n";
  
  cout << "===============================================================================\n";
  cout << "Uncorrected energies:\n";
  
  
  for(int x = 0; x<  nros; x++){
    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]);
    cout << dumc << "\n";
  }
  
  cout << "===============================================================================\n";
  cout << "Permanent dipoles in x:\n";
  

  datf.read((char *) dvals, nstates * sizeof(double));
  datf.read((char *) dmat, nstates* nstates* sizeof(double));
  for(int x = 0; x < nstates; x++) dummat[x*nstates+x] = dvals[x];
//  cout << "Hier\n";
//  exit(1);


  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
  for(int x = 0; x<  nros; x++){
    double val = 0.;
    for(int y = 0; y < nstates; y++) val += dmat[y*nstates+x]*dvals[y]*dmat[y*nstates+x];
    sprintf(dumc,"State: %.5i   Perm. dip x : %4.14f",x,val);
    cout << dumc << "\n";
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    cout << "Transition moments from " << x << "\n";
    cout << "...............................................................................\n";
      for(int y = x+1; y < nros; y++){
	double val = 0.;
	for(int z = 0; z < nstates; z++) val += dmat[z*nstates+x]*dvals[z]*dmat[z*nstates+y];
	sprintf(dumc,"%.5i -->%.5i   P. dip x : %4.14f",x,y,val);
	cout << dumc << "\n";
      }
      cout << "...............................................................................\n";
  }

  cout << "===============================================================================\n";
  cout << "Permanent dipoles in y:\n";
  

  datf.read((char *) dvals, nstates * sizeof(double));
  datf.read((char *) dmat, nstates* nstates * sizeof(double));
  for(int x = 0; x < nstates*nstates; x++) dummat[x] = 0.;
  for(int x = 0; x < nstates; x++) dummat[x*nstates+x] = dvals[x];
  
  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
  for(int x = 0; x<  nros; x++){
    double val = 0.;
    for(int y = 0; y < nstates; y++) val += dmat[y*nstates+x]*dvals[y]*dmat[y*nstates+x];
    sprintf(dumc,"State: %.5i   Perm. dip y : %4.14f",x,val);
    cout << dumc << "\n";
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    cout << "Transition moments from " << x << "\n";
    cout << "...............................................................................\n";
      for(int y = x+1; y < nros; y++){
	double val = 0.;
	for(int z = 0; z < nstates; z++) val += dmat[z*nstates+x]*dvals[z]*dmat[z*nstates+y];
	sprintf(dumc,"%.5i -->%.5i   P. dip y : %4.14f",x,y,val);
	cout << dumc << "\n";
      }
      cout << "...............................................................................\n";
  }
  
  cout << "===============================================================================\n";
  cout << "Permanent dipoles in z:\n";
  

  datf.read((char *) dvals, nstates * sizeof(double));
  datf.read((char *) dmat, nstates* nstates * sizeof(double));
  for(int x = 0; x < nstates*nstates; x++) dummat[x] = 0.;
  for(int x = 0; x < nstates; x++) dummat[x*nstates+x] = dvals[x];

  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
  for(int x = 0; x<  nros; x++){
    double val = 0.;
    for(int y = 0; y < nstates; y++) val += dmat[y*nstates+x]*dvals[y]*dmat[y*nstates+x];
    sprintf(dumc,"State: %.5i   Perm. dip z : %4.14f",x,val);
    cout << dumc << "\n";
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    cout << "Transition moments from " << x << "\n";
    cout << "...............................................................................\n";
      for(int y = x+1; y < nros; y++){
	double val = 0.;
	for(int z = 0; z < nstates; z++) val += dmat[z*nstates+x]*dvals[z]*dmat[z*nstates+y];
	sprintf(dumc,"%.5i -->%.5i   P. dip z : %4.14f",x,y,val);
	cout << dumc << "\n";
      }
      cout << "...............................................................................\n";
  }
  cout << "===============================================================================\n";  

  datf.read((char *) corr, nstates * sizeof(double));
  cout << "Corrected energies:\n";
  
  
  for(int x = 0; x<  nros; x++){
    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]+corr[x]);
    cout << dumc << "\n";
  }
  cout << "===============================================================================\n";
  
 }
