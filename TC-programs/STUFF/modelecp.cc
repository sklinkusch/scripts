# include <iostream>
# include <fstream>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 2){
	cerr << "Usage: ./modelecp <output ecp file>\n";
	exit(1);
    }
    int nros = 5;
    ofstream outf(argv[1]);
    outf.write((char *) &nros, sizeof(int));
    double* energies = new double[nros];
    double delta = 0.1;
    double eps = 0.01;
    energies[0] = 0.;
    energies[1] = delta;
    energies[2] = delta + eps;
    energies[3] = delta + (2.*eps);
    energies[4] = delta + (3.*eps);
    outf.write((char *) energies, sizeof(double)*nros);
    double* dx = new double[nros*nros];
    double* dy = new double[nros*nros];
    double* dz = new double[nros*nros];
    for(int x = 0; x < nros; x++){
	for(int y = 0; y < nros; y++){
	    dx[x*nros+y] = 0.;
	    dy[x*nros+y] = 0.;
	    dz[x*nros+y] = 0.;
	}
    }
    dx[3] = 1.;
    for(int x = 0; x < nros; x++){
	for(int y = x; y < nros; y++){
	    dx[y*nros+x] = dx[x*nros+y];
	}
    }
    outf.write((char *) dx, sizeof(double)*nros*nros);
    outf.write((char *) dy, sizeof(double)*nros*nros);
    outf.write((char *) dz, sizeof(double)*nros*nros);
    double* ion = new double[nros];
    for(int x = 0; x < nros; x++) ion[x] = 0.;
    ion[nros-1] = 0.1;
    outf.write((char *) ion, sizeof(double)*nros);
    double* omega = new double[nros*nros];
    for(int x = 0; x < nros; x++){
	for(int y = 0; y < nros; y++){
	    omega[x*nros+y] = energies[x] - energies[y];
	}
    }
    double* cpl = new double[nros*nros];
    double matsubara = 0.000237381;
    for(int x = 0; x < nros; x++){
	for(int y = 0; y < nros; y++){
	    cpl[x*nros+y] = 1.3711956e-7 / (pow(fabs(omega[x*nros+y]),2.) + pow(matsubara,2.));
	    if(x == 0) cpl[x*nros+y] = 0.;
	    if(y == 0) cpl[x*nros+y] = 0.;
	    if(x == y) cpl[x*nros+y] = 0.;
	}
    }
    outf.write((char *) cpl, sizeof(double)*nros*nros);
    double* kinens = new double[nros];
    for(int x = 0; x < nros; x++) kinens[x] = ion[x];
    outf.write((char *) kinens, sizeof(double)*nros);
    outf.flush();
    outf.close();
}


