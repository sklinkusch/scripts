/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: calc_ct_op.cc                                                          *
 *                                                                              *
 * calculates a charge transfer target operator in one cartesion dimension      * 
 * for OCT                                                                      * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 



#include <fstream>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <iostream>
#include <sstream>
#include <complex>
#include <stdlib.h>

using namespace std;
#define Complex complex<double>

//Functions


//Extern Functions
extern int    rem_com(char* filename, char* streamstring, int string_length);
extern void   status(ofstream* outf);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void  diag_mat(int nroao, double* mat, double* vals, double* vecs);

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "TDCIS [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";

  //SYSTEM
  char    bcsfile[256];      //binary system file
  int     use_d = 0;         //if 1 use the 2nd PT corrected energies for propagations
  int     wdim;              //which cartisiean dimension to use 0=x,1=y,2=z
  double  fac;
  double  beta;              //beta in target operator := mu_dim*exp(-beta(H))+exp(-beta(H))mu_dim

  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  ist >> bcsfile >> use_d >> wdim >> fac >> beta;
  
  outf << "Reading system data from " << bcsfile << "\n";
  if(use_d == 1) 
    outf << "USING d-corrections in calculation of\n";
  else 
    outf << "NO   d-corrections int calculation of\n";
  outf << "\t O = fac*mu_dim*exp(-beta(H))+fac*exp(-beta(H))mu_dim\n";
  outf << "Dim is " << wdim << " ";
  if(wdim == 0) outf << "(x)\n";
  if(wdim == 1) outf << "(y)\n";
  if(wdim == 2) outf << "(z)\n";
  if(wdim < 0 || wdim > 2){ outf << "Invalid dimension spec.!"; outf.flush(); exit(1);}
  outf << "Fac is " << fac << "\n";
  outf << "Beta is " << beta << "\n";
  
  //CIS READ IN
  ifstream datf(bcsfile);

  int nroao, nroe, llim, ulim;

  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));

  outf << "CIS data read: \n";
  int nrof = ulim - llim +1;
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  int cis_size = omo*umo+1;
  long long int Cis_size = cis_size;
  outf << "nr of CSF               : " << cis_size << "\n";
  
  double* cis_vals     = new double[cis_size];        //eigenvlaues of CI-matrix
  double* mu_vals_dim  = new double[cis_size];        //eigenvalues of dipole operator in dim
  double* dumvals      = new double[cis_size];       //tmp vals
  double* corr_vals    = new double[cis_size];          //PT corrections to CIS ex ens

  double* dumvecs     = new double[(long long int) cis_size * (long long int) cis_size];        //tmp vecs
  double* mu_vecs_dim = new double[(long long int) cis_size * (long long int) cis_size];        //eigenvectors of dipole operator in x
  double* dummat      = new double[(long long int) cis_size * (long long int) cis_size];        //tmp mat
  
  datf.read((char *) cis_vals, cis_size * sizeof(double));
  datf.read((char *) dumvecs, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "...............................................................................\n";
  outf << "First uncorrected excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
  outf << "Highest uncorrected excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";
  
  for(int x = 0; x < 3; x++){
    if(x==wdim){
        datf.read((char *) mu_vals_dim, cis_size * sizeof(double));
	datf.read((char *) mu_vecs_dim, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "...............................................................................\n";
	outf << "Dipole eigenvalues in dim range from " << mu_vals_dim[0] << " to " << mu_vals_dim[cis_size-1] << "\n";
    }else{
      datf.read((char *) dumvals, cis_size * sizeof(double));
      datf.read((char *) dumvecs, (long long int) cis_size * (long long int) cis_size * sizeof(double));
    }
  }
  outf.flush();

  if(use_d == 1){
  outf << "...............................................................................\n";
    outf << "Reading PT-correction\n";
    datf.read((char *) corr_vals, cis_size * sizeof(double));
    for(int x = 0; x < cis_size; x++)
      cis_vals[x] +=  corr_vals[x];
    outf << "First CORRECTED excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
    outf << "Highest CORRECTED excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";
  }
  outf << "-------------------------------------------------------------------------------\n";

  datf.close();
 
  outf << "Transforming fac*(dipole matrix) to CIS-eigenspace\n"; outf.flush();
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = fac*mu_vals_dim[x];
  if(cis_size < 10000)
    trans_mat((int) cis_size, dummat, mu_vecs_dim, dumvecs, 0); //Forward transfornation
  else
    trans_mat(cis_size, dummat, mu_vecs_dim, dumvecs, 0); //Forward transfornation
  
  outf << "Done!\n\n"; outf.flush();
  
  
  outf << "First couple of moments\n";
  
  int nprint =  8; 
  if(nprint > cis_size) nprint = cis_size;
  
  for(int x = 0; x < nprint; x++) outf << "        " << x;
  outf << "\n";
  
  for(long long int x = 0; x < nprint ; x++ ){
    char tmpc[12];
    sprintf(dumc,"%i   ",(int) x);
    for(long long int y = 0; y < nprint; y++) {
      sprintf(tmpc,"%+.4f  ", dummat[x*cis_size+y]);
      strcat(dumc,tmpc);
    }
    outf << dumc << "\n";
  }
  outf << "\n\n";
  outf.flush();
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Building operator in CIS-Eigenspace\n";
  //Calc exp(-betaH)
  for(int x = 0; x < cis_size; x++){
    //Old version
    //dumvals[x] = exp(-beta*cis_vals[x]);
  
    //new version
    if(cis_vals[x] > beta) dumvals[x] = 0.;
    else dumvals[x] = 1.;
  }
   
  //build op in dumvecs !!
  for(int x = 0; x < cis_size; x++){
    for(int y = 0; y < cis_size; y++)
      dumvecs[x*Cis_size+y] = dummat[x*Cis_size+y]*(dumvals[y]+dumvals[x]);
  }
  
  outf << "Calculating eigenvalues of the operator\n";
  //eigenvalues in corr_vals eigenvectors in  dummat
  diag_mat(cis_size, dumvecs, corr_vals,  dummat);
  
  outf << "Lowest eigen value is " << corr_vals[0] << "\n";
  outf << "Highest eigen value is " << corr_vals[cis_size-1] << "\n";
  
  outf << "Shifting diagonal elements by " << corr_vals[0] << "\n";
  for(int x = 0; x < cis_size; x++)
    dumvecs[x*Cis_size+x] -= corr_vals[0];
  
  outf << "Calculating eigenvalues and eigenvectors  of the operator\n";
  //eigenvalues in corr_vals eigenvectors in  dummat
  diag_mat(cis_size, dumvecs, corr_vals,  dummat);
  outf << "Lowest eigen value is " << corr_vals[0] << "\n";
  outf << "-------------------------------------------------------------------------------\n";
  outf << "First 10 highest  eigen values:\n";
  outf << "...............................................................................\n";
  for(int v = 1; v <= 10; v++){
    outf <<  v << " highest  eigen values is " << corr_vals[cis_size-v] << "\n";
    double en = 0.;
    for(int x = 0; x < cis_size; x++)
      en += cis_vals[x]*pow(dummat[(cis_size-v)*Cis_size+x],2);
    outf << "Energy is " << en << "\n\n";
    outf << "Populations in CIS-eigen state basis over 1 %\n";
    for(int x = 0; x < cis_size; x++){
      if(pow(dummat[(cis_size-v)*Cis_size+x],2) > 0.01){
	outf << x << " " << pow(dummat[(cis_size-v)*Cis_size+x],2) << "\n";
      }
    }	
    outf << "...............................................................................\n";
  }
  outf << "-------------------------------------------------------------------------------\n";
  
  sprintf(dumc,"%s.top",argv[2]);
  ofstream topf(dumc);
  outf << "Writing target operator to " << dumc << "\n";
  topf.write((char *) dumvecs, Cis_size*Cis_size*sizeof(double));
  outf << "-------------------------------------------------------------------------------\n";
  status(&outf);
  outf.flush();
}

