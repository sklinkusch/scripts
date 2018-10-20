#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./weighted-dipole <ecp file> <pop file> <output prefix>\n";
	exit(1);
    }
    char dumc[512];
    double dumv;

    ifstream inf(argv[1]);
    int nros;
    inf.read((char *) &nros, sizeof(int));
    double* ens = new double[nros];
    double* dx  = new double[nros*nros];
    double* dy  = new double[nros*nros];
    double* dz  = new double[nros*nros];
    double* ion = new double[nros];
    inf.read((char *) ens, sizeof(double)*nros);
    inf.read((char *) dx, sizeof(double)*nros*nros);
    inf.read((char *) dy, sizeof(double)*nros*nros);
    inf.read((char *) dz, sizeof(double)*nros*nros);
    inf.read((char *) ion, sizeof(double)*nros);
    double t, popion;
    double* pop = new double[nros];
    double* wdip = new double[nros];

    ifstream popf(argv[2]);
    popf >> t;
    for(int x = 0; x < nros; x++){
	popf >> pop[x];
    }
    popf >> popion;

    double* dip = new double[nros*nros];
#pragma omp parallel for
    for(int x = 0; x < nros; x++){
	for(int y = x; y < nros; y++){
	    dip[x*nros+y] = sqrt(pow(dx[x*nros+y],2.)+pow(dy[x*nros+y],2.)+pow(dz[x*nros+y],2.));
	    dip[y*nros+x] = dip[x*nros+y];
	}
    }

#pragma omp parallel for reduction(+:dumv)
    for(int y = 0; y < nros; y++){
	dumv = 0.;
	if(ion[y] != 0.){
	    for(int x = 0; x < nros; x++){
		if(ens[x] < ens[y]){
		    dumv += pop[x]*dip[x*nros+y];
		}
	    }
	}
	wdip[y] = dumv;
    }

    ofstream outf;
    sprintf(dumc, "%s.wdip", argv[3]);
    outf.open(dumc);
    for(int y = 0; y < nros; y++){
	outf << y << " " << ens[y] << " " << wdip[y] << "\n";
    }
    outf.close();
}

