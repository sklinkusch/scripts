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
  if(argc != 2){
    cerr << "Nedd <bcs-file>\n";
    exit(0);
  }

    
//  int nros = atoi(argv[1]);

//  cout << "Printing out information for " << nros << " states from " << argv[2] << "\n";

//  char dumc[1024];
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

//  cout << "Nr of CIS states is " << cis_size << "\n";
  
  double* ens   = new double[cis_size];
  double* vecs  = new double[cis_size*cis_size];
  double* dmat  = new double[cis_size*cis_size];
  double* dvals = new double[cis_size];
//  double* dummat= new double[cis_size*cis_size];
  double* corr  = new double[cis_size];
  
  datf.read((char *) ens, cis_size * sizeof(double));
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));
  
//  cout << "===============================================================================\n";
//  cout << "Uncorrected energies:\n";
  
  
//  for(int x = 0; x<  nros; x++){
//    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]);
//    cout << dumc << "\n";
//  }
  
//  cout << "===============================================================================\n";
//  cout << "Permanent dipoles in x:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
//  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];


  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
//  for(int x = 0; x<  nros; x++){
//    double val = 0.;
//    for(int y = 0; y < cis_size; y++) val += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
//    sprintf(dumc,"State: %.5i   Perm. dip x : %4.14f",x,val);
//    cout << dumc << "\n";
//  }
//  cout << "-------------------------------------------------------------------------------\n";  
//  for(int x = 0; x<  nros-1; x++){
//    cout << "Transition moments from " << x << "\n";
//    cout << "...............................................................................\n";
//      for(int y = x+1; y < nros; y++){
//	double val = 0.;
//	for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
//	sprintf(dumc,"%.5i -->%.5i   P. dip x : %4.14f",x,y,val);
//	cout << dumc << "\n";
//      }
//      cout << "...............................................................................\n";
//  }

//  cout << "===============================================================================\n";
//  cout << "Permanent dipoles in y:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
//  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
//  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];
  
  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
//  for(int x = 0; x<  nros; x++){
//    double val = 0.;
//    for(int y = 0; y < cis_size; y++) val += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
//    sprintf(dumc,"State: %.5i   Perm. dip y : %4.14f",x,val);
//    cout << dumc << "\n";
//  }
//  cout << "-------------------------------------------------------------------------------\n";  
//  for(int x = 0; x<  nros-1; x++){
//    cout << "Transition moments from " << x << "\n";
//    cout << "...............................................................................\n";
//      for(int y = x+1; y < nros; y++){
//	double val = 0.;
//	for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
//	sprintf(dumc,"%.5i -->%.5i   P. dip y : %4.14f",x,y,val);
//	cout << dumc << "\n";
//      }
//      cout << "...............................................................................\n";
//  }
  
//  cout << "===============================================================================\n";
//  cout << "Permanent dipoles in z:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
//  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
//  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

  //OLD SLOW VERSION  
//   if(cis_size < 10000)
//     trans_mat((int) cis_size, dummat, dmat, vecs, 0); //Forward transfornation
//   else
//     trans_mat(cis_size, dummat, dmat, vecs, 0); //Forward transfornation
  
//  for(int x = 0; x<  nros; x++){
//    double val = 0.;
//    for(int y = 0; y < cis_size; y++) val += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
//    sprintf(dumc,"State: %.5i   Perm. dip z : %4.14f",x,val);
//    cout << dumc << "\n";
//  }
//  cout << "-------------------------------------------------------------------------------\n";  
//  for(int x = 0; x<  nros-1; x++){
//    cout << "Transition moments from " << x << "\n";
//    cout << "...............................................................................\n";
//      for(int y = x+1; y < nros; y++){
//	double val = 0.;
//	for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
//	sprintf(dumc,"%.5i -->%.5i   P. dip z : %4.14f",x,y,val);
//	cout << dumc << "\n";
//      }
//      cout << "...............................................................................\n";
//  }
//  cout << "===============================================================================\n";  

  datf.read((char *) corr, cis_size * sizeof(double));
//  cout << "Corrected energies:\n";
  cout << corr[0] << "\n";
  
//  for(int x = 0; x<  nros; x++){
//    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]+corr[x]);
//    cout << dumc << "\n";
//  }
//  cout << "===============================================================================\n";
  
  datf.close();
  
  

 }

