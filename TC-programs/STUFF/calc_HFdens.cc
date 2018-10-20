#include <iostream>
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
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Need input-prefix, maximal input value, output-filename!\n";
	exit(0);
    }
    char filename[128];
    int nroa, nroa_rd, dim_rd;
    int imax = atoi(argv[2]);
    double* xgrid_rd = new double[3];
    double* xgrid = new double[3];
    int* nrop_rd = new int[3];
    int* nrop = new int[3];
    double* sgrid_rd = new double[9];
    double* sgrid = new double[9];
    sprintf(filename,"%s%i.bcb",argv[1],1);
    ifstream inf;
    inf.open(filename);
    inf.read((char *) &nroa_rd, sizeof(int));
    nroa = nroa_rd;
    inf.read((char *) xgrid_rd, sizeof(double)*3);
    inf.read((char *) nrop_rd, sizeof(int)*3);
    dim_rd = nrop_rd[0] * nrop_rd[1] * nrop_rd[2];
    inf.read((char *) sgrid_rd, sizeof(double)*9);
    int* charge_rd = new int[nroa];
    double* dumv_rd = new double[nroa];
    double* coord_rd = new double[3*nroa];
    long double* dens_rd = new long double[dim_rd];
    int* charge = new int[nroa];
    double* dumv = new double[nroa];
    double* coord = new double[3*nroa];
    long double* dens = new long double[dim_rd];
    for(int j = 0; j < 3; j++){
	xgrid[j] = xgrid_rd[j];
	nrop[j] = nrop_rd[j];
    }
    for(int j = 0; j < 9; j++) sgrid[j] = sgrid_rd[j];
    inf.read((char *) charge_rd, sizeof(int)*nroa);
    inf.read((char *) dumv_rd, sizeof(double)*nroa);
    for(int j = 0; j < nroa; j++){
	charge[j] = charge_rd[j];
	dumv[j] = dumv_rd[j];
    }
    inf.read((char *) coord_rd, sizeof(double)*3*nroa);
    for(int j = 0; j < (3*nroa); j++) coord[j] = coord_rd[j];
    inf.read((char *) dens_rd, sizeof(long double)*dim_rd);
    inf.close();
    for(int j = 0; j < dim_rd; j++) dens[j] = pow(dens_rd[j],2.);
    for(int i = 2; i <= imax; i++){
	sprintf(filename,"%s%i.bcb",argv[1],i);
	ifstream inf;
	inf.open(filename);
	inf.read((char *) &nroa_rd, sizeof(int));
	inf.read((char *) xgrid_rd, sizeof(double)*3);
	inf.read((char *) nrop_rd, sizeof(int)*3);
	dim_rd = nrop_rd[0] * nrop_rd[1] * nrop_rd[2];
	inf.read((char *) sgrid_rd, sizeof(double)*9);
	inf.read((char *) charge_rd, sizeof(int)*nroa);
	inf.read((char *) dumv_rd, sizeof(double)*nroa);
	inf.read((char *) coord_rd, sizeof(double)*3*nroa);
	if(nroa_rd != nroa){
	    cerr << "Nr of atoms is not equal\n";
	    exit(1);
	}
	if(xgrid_rd[0] != xgrid[0] || xgrid_rd[1] != xgrid[1] || xgrid_rd[2] != xgrid[2]){
	    cerr << "Size of (extra) grid not equal\n";
	    exit(2);
	}
	if(nrop_rd[0] != nrop[0] || nrop_rd[1] != nrop[1] || nrop_rd[2] != nrop[2]){
	    cerr << "Nr of grid points not equal\n";
	    exit(3);
	}
	for(int j = 0; j < 9; j++){
	    if(sgrid_rd[j] != sgrid[j]){
		cerr << "Stepwidth of grid points not equal\n";
		exit(4);
	    }
	}
	for(int j = 0; j < nroa; j++){
	    if(charge_rd[j] != charge[j]){
		cerr << "Nuclear charges not equal (" << charge_rd[j] << ", " << charge[j] << "\n";
		exit(5);
	    }
	    if(dumv_rd[j] != dumv[j]){
		cerr << "System data not equal\n";
		exit(6);
	    }
	}
	for(int j = 0; j < (3*nroa); j++){
	    if(coord_rd[j] != coord[j]){
		cerr << "Molecular coordinates not equal\n";
	    }
	}
	inf.read((char *) dens_rd, sizeof(long double)*dim_rd);
	inf.close();
	for(int j = 0; j < dim_rd; j++) dens[j] += pow(dens_rd[j],2.);
    }
    char outputfile[128];
    sprintf(outputfile,"%s",argv[3]);
    ofstream outf(outputfile);
    outf.write((char *) &nroa, sizeof(int));
    outf.write((char *) xgrid, sizeof(double)*3);
    outf.write((char *) nrop, sizeof(int)*3);
    outf.write((char *) sgrid, sizeof(double)*9);
    outf.write((char *) charge, sizeof(int)*nroa);
    outf.write((char *) dumv, sizeof(double)*nroa);
    outf.write((char *) coord, sizeof(double)*3*nroa);
    outf.write((char *) dens, sizeof(long double)*dim_rd);
    outf.flush();
}

