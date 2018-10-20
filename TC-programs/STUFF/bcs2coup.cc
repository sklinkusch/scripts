#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Usage: ./bcs2coup <nr of printout states> <eta> <omega_c> <bcs-file> <outfile>\n";
	exit(0);
    }

    int nros = atoi(argv[1]);
    double eta = strtod(argv[2], NULL);
    double omega_c = strtod(argv[3], NULL);
    char dumc[1024];
    int nroao, nroe, llim, ulim;
    ifstream datf(argv[4]);
    ofstream outf;
    outf.open(argv[5]);

    datf.read((char *) &nroao, sizeof(int));
    datf.read((char *) &nroe, sizeof(int));
    datf.read((char *) &llim, sizeof(int));
    datf.read((char *) &ulim, sizeof(int));

    if(llim < 0 || llim >= (nroe/2)){
	cerr << "Lower limit of orbitals out of range!\n";
	exit(1);
    }
    if(ulim < (nroe/2) || ulim >= nroao){
	cerr << "Upper limit of orbitals out of range!\n";
	exit(2);
    }
    if(nroao < (nroe/2)){
	cerr << "Wrong number of atomic orbitals!\n";
	exit(3);
    }

    int omo = nroe/2 - llim;
    int umo = ulim - nroe/2 + 1;
    long long int cis_size = omo * umo + 1;
    if(nros > cis_size) nros = cis_size;

    double val_x = 0.;
    double val_y = 0.;
    double val_z = 0.;
    double* ens  = new double[cis_size];
    double* corr = new double[cis_size];
    double* correns = new double[cis_size];
    double* omega = new double[cis_size*cis_size];
    double* vecs = new double[cis_size*cis_size];
    double* muvals_x = new double[cis_size];
    double* muvecs_x = new double[cis_size*cis_size];
    double* muvals_y = new double[cis_size];
    double* muvecs_y = new double[cis_size*cis_size];
    double* muvals_z = new double[cis_size];
    double* muvecs_z = new double[cis_size*cis_size];
    double* rates = new double[cis_size*cis_size];
    double* dephasing = new double[cis_size*cis_size];

    double* mumat_x = new double[cis_size*cis_size];
    double* mumat_y = new double[cis_size*cis_size];
    double* mumat_z = new double[cis_size*cis_size];

    datf.read((char *) ens, sizeof(double)*cis_size);
    datf.read((char *) vecs, sizeof(double)*cis_size*cis_size);
    datf.read((char *) muvals_x, sizeof(double)*cis_size);
    datf.read((char *) muvecs_x, sizeof(double)*cis_size*cis_size);
    datf.read((char *) muvals_y, sizeof(double)*cis_size);
    datf.read((char *) muvecs_y, sizeof(double)*cis_size*cis_size);
    datf.read((char *) muvals_z, sizeof(double)*cis_size);
    datf.read((char *) muvecs_z, sizeof(double)*cis_size*cis_size);
    datf.read((char *) corr, sizeof(double)*cis_size);
#pragma omp parallel for reduction(+:val_x,val_y,val_z)
    for(int x = 0; x < cis_size; x++){
	val_x = 0.;
	val_y = 0.;
	val_z = 0.;
	for(int y = 0; y < cis_size; y++){
	    val_x += muvecs_x[y*cis_size+x] * muvals_x[y] * muvecs_x[y*cis_size+x];
	    val_y += muvecs_y[y*cis_size+x] * muvals_y[y] * muvecs_y[y*cis_size+x];
	    val_z += muvecs_z[y*cis_size+x] * muvals_z[y] * muvecs_z[y*cis_size+x];
	}
	correns[x] = ens[x] + corr[x];
	mumat_x[x*cis_size+x] = val_x;
	mumat_y[x*cis_size+x] = val_y;
	mumat_z[x*cis_size+x] = val_z;
    }
#pragma omp parallel for reduction(+:val_x,val_y,val_z)
    for(int x = 0; x < cis_size-1; x++){
	for(int y = x+1; y < cis_size; y++){
	    val_x = 0.;
	    val_y = 0.;
	    val_z = 0.;
	    for(int z = 0; z < cis_size; z++){
		val_x += muvecs_x[z*cis_size+x] * muvals_x[z] * muvecs_x[z*cis_size+y];
		val_y += muvecs_y[z*cis_size+x] * muvals_y[z] * muvecs_y[z*cis_size+y];
		val_z += muvecs_z[z*cis_size+x] * muvals_z[z] * muvecs_z[z*cis_size+y];
	    }
	    mumat_x[x*cis_size+y] = val_x;
	    mumat_x[y*cis_size+x] = val_x;
	    mumat_y[x*cis_size+y] = val_y;
	    mumat_y[y*cis_size+x] = val_y;
	    mumat_z[x*cis_size+y] = val_z;
	    mumat_z[y*cis_size+x] = val_z;
	}
    }
#pragma omp parallel for
    for(int x = 0; x < cis_size; x++){
	for(int y = 0; y < cis_size; y++){
	    omega[x*cis_size+y] = fabs(correns[x] - correns[y]);
	}
    }

    double dep_pref = 0.00001 / (omega[1]*omega[1]);
#pragma omp parallel for
    for(int x = 0; x < cis_size; x++){
	for(int y = 0; y < cis_size; y++){
	    dephasing[x*cis_size+y] = dep_pref * pow(omega[x*cis_size+y],2);
	}
    }
#pragma omp parallel for
    for(int x = 0; x < cis_size; x++){
	for(int y = 0; y < cis_size; y++){
	    if(x != y){
		rates[x*cis_size+y] = eta*omega[x*cis_size+y]*exp(-1.*omega[x*cis_size+y]/omega_c);
	    }else{
		rates[x*cis_size+y] = 0.;
	    }
	}
    }

    for(int j = 0; j < nros; j++){
	for(int i = 0; i <= j; i++){
	    sprintf(dumc, "%4i  %4i   % 12.10f    % 12.10f     % 12.10f      % 22.20f    % 12.10f", i, j, mumat_x[i*cis_size+j], mumat_y[i*cis_size+j], mumat_z[i*cis_size+j], rates[i*cis_size+j], dephasing[i*cis_size+j]);
	    outf << dumc << "\n";
	    outf.flush();
	}
    }
    outf.close();
}

