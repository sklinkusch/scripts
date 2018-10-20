#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f);
void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./sys2ecp <sysfile> <ip> <ecpfile>\n";
	exit(1);
    }
    ifstream inf(argv[1]);
    int nbasis, nstates;
    inf.read((char *) &nbasis, sizeof(int));
    inf.read((char *) &nstates, sizeof(int));
    double* configurations = new double[nbasis];
    double* gamma_csf = new double[nbasis];
    double* dx_csf = new double[nbasis*nbasis];
    double* dy_csf = new double[nbasis*nbasis];
    double* dz_csf = new double[nbasis*nbasis];
    double* rates = new double[nbasis*nbasis];
    double* dephasing = new double[nbasis*nbasis];
    double* energies = new double[nstates];
    double* cisvecs = new double[nbasis*nstates];
    double* tmpmat = new double[nbasis*nstates];
    inf.read((char *) configurations, sizeof(double)*nbasis);
    inf.read((char *) gamma_csf, sizeof(double)*nbasis);
    inf.read((char *) dx_csf, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dy_csf, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dz_csf, sizeof(double)*nbasis*nbasis);
    inf.read((char *) rates, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dephasing, sizeof(double)*nbasis*nbasis);
    inf.read((char *) energies, sizeof(double)*nstates);
    inf.read((char *) cisvecs, sizeof(double)*nbasis*nstates);
    ofstream outf;
    outf.open(argv[3]);
    outf.write((char *) &nstates, sizeof(int));
    outf.write((char *) energies, sizeof(double)*nstates);
    trans_mat(nstates, dx_csf, cisvecs, tmpmat, 1);
    trans_mat(nstates, dy_csf, cisvecs, tmpmat, 1);
    trans_mat(nstates, dz_csf, cisvecs, tmpmat, 1);
    double* dx = new double[nstates*nstates];
    double* dy = new double[nstates*nstates];
    double* dz = new double[nstates*nstates];
    for(int x = 0; x < nstates; x++){
	for(int y = 0; y< nstates; y++){
	dx[x*nstates+y] = dx_csf[x*nstates+y];
	dy[x*nstates+y] = dy_csf[x*nstates+y];
	dz[x*nstates+y] = dz_csf[x*nstates+y];
	}
    }
    outf.write((char *) dx, sizeof(double)*nstates*nstates);
    outf.write((char *) dy, sizeof(double)*nstates*nstates);
    outf.write((char *) dz, sizeof(double)*nstates*nstates);
    double ip = strtod(argv[2], NULL);
    double* ion = new double[nstates];
    double dumv;
    for(int x = 0; x < nstates; x++){
	if(energies[x] > ip){
	    dumv = 0.;
	    for(int y = 0; y < nstates; y++){
		dumv += pow(fabs(cisvecs[x*nstates+y]),2.)*gamma_csf[y];
	    }
	    ion[x] = dumv;
	}else{
	    ion[x] = 0.;
	}
    }
    outf.write((char *) ion, sizeof(double)*nstates);
    double* omega = new double[nstates*nstates];
    double* mutot = new double[nstates*nstates];
    for(int x = 0; x < nstates; x++){
	for(int y = 0; y < nstates; y++){
	    omega[x*nstates+y] = energies[x] - energies[y];
	    mutot[x*nstates+y] = sqrt(pow(dx_csf[x*nstates+y],2.)+pow(dy_csf[x*nstates+y],2.)+pow(dz_csf[x*nstates+y],2.));
	}
    }
    double* cpl = new double[nstates*nstates];
    for(int x = 0; x < nstates; x++){
	for(int y = 0; y < nstates; y++){
	    cpl[x*nstates+y] = (4.*pow(fabs(mutot[x*nstates+y]),2.)*pow(fabs(omega[x*nstates+y]),3.))/(3.*pow(137.036,3.));
	}
    }
    outf.write((char *) cpl, sizeof(double)*nstates*nstates);
    outf.flush();
}

void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir){

  if(dir == 0) //CALC mat*trans
    mat_mat( np,  mat,trans, tmpmat);
  else         //CALC trans*mat 
    mat_mat( np, trans, mat, tmpmat);
  
  
  if(dir == 0) //CALC trans^T*(mat*trans)
    mat_T_mat( np, trans, tmpmat, mat);
  else         //CALC (trans*mat)*trans^T
    mat_mat_T( np, tmpmat, trans, mat);
  
}

void mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[z*np+y];
      }
    }
  }
}


void mat_mat_T(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[x*np+z]*mat_i2[y*np+z];
      }
    }
  }
}


void mat_T_mat(int np, double* mat_i1, double* mat_i2, double* mat_f){
  for(int x = 0; x < np; x++){
    for(int y = 0; y < np; y++){
      mat_f[x*np+y] = 0.;
      for(int z = 0; z < np; z++){
        mat_f[x*np+y] += mat_i1[z*np+x]*mat_i2[z*np+y];
      }
    }
  }
}

