/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cispektrum.cc                                                          *
 *                                                                              *
 *                                                                              * 
 *                                                                              *
 *                                                    Stefan Klinkusch   2011   *
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
    cerr << "Need <bcs-file> <sigma> \n";
    exit(0);
  }

    
  //char dumc[1024];
  int nroao, nroe, llim, ulim;
  ifstream datf(argv[1]);
  double sigma = strtod(argv[2],NULL);

  
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));
  
  int omo                = nroe/2 - llim;
  int umo                = ulim - nroe/2+1;
  long long int cis_size = omo*umo+1;
  int nros = (int) cis_size;


  double* ens   = new double[cis_size];
  double* vecs  = new double[cis_size*cis_size];
  double* dmat  = new double[cis_size*cis_size];
  double* dvals = new double[cis_size];
  double* dummat= new double[cis_size*cis_size];
  double* corr  = new double[cis_size];
  double* muq   = new double[cis_size];
  double* exv   = new double[cis_size];
  double* ozstr = new double[cis_size];


  for(int i = 0; i < cis_size; i++){
   muq[i] = 0.;
  }

  
  datf.read((char *) ens, cis_size * sizeof(double));
  datf.read((char *) vecs, cis_size* cis_size * sizeof(double));

 
 int x = 0;
  for(int y = x+1; y < nros; y++){
   double val = 0.;
   for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
   muq[y] += pow(val,2);
  }


  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

 x = 0;
  for(int y = x+1; y < nros; y++){
   double val = 0.;
   for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
   muq[y] += pow(val,2);
  }


  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

 x = 0;
  for(int y = x+1; y < nros; y++){
   double val = 0.;
   for(int z = 0; z < cis_size; z++) val += dmat[z*cis_size+x]*dvals[z]*dmat[z*cis_size+y];
   muq[y] += pow(val,2);
  }


  datf.read((char *) dvals, cis_size * sizeof(double));
  datf.read((char *) dmat, cis_size* cis_size * sizeof(double));
  for(int x = 0; x < cis_size; x++) dummat[x*cis_size+x] = dvals[x];

  datf.read((char *) corr, cis_size * sizeof(double));

  for(int x = 0; x < cis_size; x++){
   ens[x] += corr[x];
   exv[x] = 219474.63*ens[x];
  }


  for(int x = 0; x < cis_size; x++){
   ozstr[x] = (2./3.)*ens[x]*muq[x];
  }
  datf.close();


  double h = 2*M_PI;
  double c0 = 137.036;
  double l_min = 0.;
  double l_max = 0.052917721*h*c0/ens[1] + 400;
  double incr = 0.002;
  int nrop = (int) ((l_max-l_min)/incr) + 1;
  double epsilon;
  double cl, cv;

  double dl = incr;
  double kappa = 4.3189984e-10;
  double nfac = sigma*sqrt(2*M_PI);
  double prefac = 0.1/(kappa*nfac);


  for(int i = 0; i < nrop; i++){
   epsilon = 0.;
   cl = l_min + i*dl;
   cv = 1./cl*10000000.;
   for(int s = 1; s < cis_size; s++){
    epsilon += ozstr[s]*prefac*exp(-0.5*pow((cv-exv[s])/sigma,2));
   }
   cout << cl << " " << epsilon << "\n"; 
  }
  

 }

