/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cut_bcs.cc                                                             *
 *                                                                              *
 * cuts a selected nr of states from a bcs file                                 * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 
#include <fstream>
#include <iostream> 
#include <stdlib.h>
#include <stdio.h>

extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
using namespace std;

int main(int argc, char* argv[]){
  if(argc != 4){
    cerr << "Nedd <nr of print out states> <input-bcs-file> <output-bcs-file>\n";
    exit(0);
  }

  int nros = atoi(argv[1]);
  int nprint = 7;
  if(nros < nprint) nprint = nros;

  cout << "Cutting  out information for " << nros << " states from " << argv[2]
       << "  writing new bcs file to : " << argv[3] << "\n";

  char dumc[1024];
  int nroao, nroe, llim, ulim;
  ifstream datf(argv[2]);
  ofstream outf(argv[3]);
  
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  
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
  

  double* new_vals = new double[nros];
  double* new_vecs = new double[nros*nros];
  double* new_mat  = new double[nros*nros];

  datf.read((char *) ens, cis_size * sizeof(double));
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));

  cout << "===============================================================================\n";
  cout << "Uncorrected energies:\n";
  
  
  for(int x = 0; x<  nprint; x++){
    sprintf(dumc,"State: %.5i   Exc. Energy: %4.14f",x,ens[x]);
    cout << dumc << "\n";
  }
  
  // we take 2 electrons in nros MOs with llim 0 ulim nros-1 to get the correct cis_size
  outf.write((char* ) &nros, sizeof(int));   
  int dumn = 2; 
  outf.write((char* ) &dumn, sizeof(int));
  dumn = 0;
  outf.write((char* ) &dumn, sizeof(int));
  dumn = nros-1;
  outf.write((char* ) &dumn, sizeof(int));
  //Write the energies
  outf.write((char* ) ens, nros* sizeof(double));
  //write zero vectors
  for(int x = 0; x < nros*nros; x++) dummat[x] = 0.;
  outf.write((char* ) dummat, nros*nros*sizeof(double));
  
  cout.precision(12);

  cout << "===============================================================================\n";
  cout << "Permanent dipoles in x:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];
  
  for(int x = 0; x<  nros; x++){
    new_mat[x*nros+x] = 0.;
    for(int y = 0; y < cis_size; y++) new_mat[x*nros+x] += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
    if(x < nprint){
      sprintf(dumc,"State: %.5i   Perm. dip x : %4.14f",x,new_mat[x*nros+x]);
      cout << dumc << "\n";
    }
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    if(x < nprint){
      cout << "...............................................................................\n";
      cout << "Transition moments from " << x << "\n";
      cout << "...............................................................................\n";
    }
    for(int y = x+1; y < nros; y++){
      new_mat[x*nros+y] = 0.;
      for(int z = 0; z < cis_size; z++) new_mat[x*nros+y] += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
      if(y < nprint){
	sprintf(dumc,"%.5i -->%.5i   P. dip x : %4.14f",x,y,new_mat[x*nros+y]);
	cout << dumc << "\n";
      }
    }
  }
  cout << "...............................................................................\n";
  
  
  for(int x = 0; x < nros; x++){
    for(int y = x; y < nros; y++) 
      new_mat[y*nros+x] = new_mat[x*nros+y];
  }
  diag_mat(nros, new_mat, new_vals, new_vecs);
  
  cout << "Calculating dipole eigenvectors \n";
  cout << "Values ranging from " << new_vals[0] << " to " << new_vals[nros-1] << "\n";
  outf.write((char *) new_vals, nros*sizeof(double));
  outf.write((char *) new_vecs, nros*nros*(long long int) sizeof(double));
  outf.flush();


  cout << "===============================================================================\n";
  cout << "Permanent dipoles in y:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];
  
  for(int x = 0; x<  nros; x++){
    new_mat[x*nros+x] = 0.;
    for(int y = 0; y < cis_size; y++) new_mat[x*nros+x] += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
    if(x < nprint){
      sprintf(dumc,"State: %.5i   Perm. dip y : %4.14f",x,new_mat[x*nros+x]);
      cout << dumc << "\n";
    }
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    if(x < nprint){
      cout << "...............................................................................\n";
      cout << "Transition moments from " << x << "\n";
      cout << "...............................................................................\n";
    }
    for(int y = x+1; y < nros; y++){
      new_mat[x*nros+y] = 0.;
      for(int z = 0; z < cis_size; z++) new_mat[x*nros+y] += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
      if(y < nprint){
	sprintf(dumc,"%.5i -->%.5i   P. dip y : %4.14f",x,y,new_mat[x*nros+y]);
	cout << dumc << "\n";
      }
    }
  }
  cout << "...............................................................................\n";
  
  
  for(int x = 0; x < nros; x++){
    for(int y = x; y < nros; y++) 
      new_mat[y*nros+x] = new_mat[x*nros+y];
  }
  diag_mat(nros, new_mat, new_vals, new_vecs);
  
  cout << "Calculating dipole eigenvectors \n";
  cout << "Values ranging from " << new_vals[0] << " to " << new_vals[nros-1] << "\n";
  outf.write((char *) new_vals, nros*sizeof(double));
  outf.write((char *) new_vecs, nros*nros*(long long int) sizeof(double));
  outf.flush();


  cout << "===============================================================================\n";
  cout << "Permanent dipoles in z:\n";
  

  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size*cis_size; x++) dummat[x] = 0.;
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];
  
  for(int x = 0; x<  nros; x++){
    new_mat[x*nros+x] = 0.;
    for(int y = 0; y < cis_size; y++) new_mat[x*nros+x] += dmat[y*cis_size+x]*dvals[y]*dmat[y*cis_size+x];
    if(x < nprint){
      sprintf(dumc,"State: %.5i   Perm. dip z : %4.14f",x,new_mat[x*nros+x]);
      cout << dumc << "\n";
    }
  }
  cout << "-------------------------------------------------------------------------------\n";  
  for(int x = 0; x<  nros-1; x++){
    if(x < nprint){  
      cout << "...............................................................................\n";
      cout << "Transition moments from " << x << "\n";
      cout << "...............................................................................\n";
    }
    for(int y = x+1; y < nros; y++){
      new_mat[x*nros+y] = 0.;
      for(int z = 0; z < cis_size; z++) new_mat[x*nros+y] += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
      if(y < nprint){
	sprintf(dumc,"%.5i -->%.5i   P. dip z : %4.14f",x,y,new_mat[x*nros+y]);
	cout << dumc << "\n";
      }
    }
  }
  cout << "...............................................................................\n";


  
  for(int x = 0; x < nros; x++){
    for(int y = x; y < nros; y++) 
      new_mat[y*nros+x] = new_mat[x*nros+y];
  }
  diag_mat(nros, new_mat, new_vals, new_vecs);
  
  cout << "Calculating dipole eigenvectors \n";
  cout << "Values ranging from " << new_vals[0] << " to " << new_vals[nros-1] << "\n";
  outf.write((char *) new_vals, nros*sizeof(double));
  outf.write((char *) new_vecs, nros*nros*(long long int) sizeof(double));
  outf.flush();

  
  
  
//   datf.close();
  
}
