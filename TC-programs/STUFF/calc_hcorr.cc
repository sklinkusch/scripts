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

void check_sizes(int nroa_rd, int nroa, double* xgrid_rd, double* xgrid, int* nrop_rd, int* nrop, double* sgrid_rd, double* sgrid, int* charge_rd, int* charge, double* dumv_rd, double* dumv, double* coord_rd, double* coord);

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Need HFdensity file, input-prefix, output-prefix, llim, nroe!\n";
	exit(0);
    }
// read Hartree-Fock density
    char HFfilename[128];
    char input_pref[128];
    char output_pref[128];
    char filename[128];
    char outname[128];
    sprintf(HFfilename,"%s",argv[1]);
    sprintf(input_pref,"%s",argv[2]);
    sprintf(output_pref,"%s",argv[3]);
    int llim = atoi(argv[4]);
    int nroe = atoi(argv[5]);
    if(nroe%2 != 0){
	cerr << "Nr of electrons must be even!\n";
	exit(11);
    }
    int homo = nroe/2;
    int nroa, nroa_rd, dim_rd;
    double* xgrid_rd = new double[3];
    double* xgrid = new double[3];
    int* nrop_rd = new int[3];
    int* nrop = new int[3];
    double* sgrid_rd = new double[9];
    double* sgrid = new double[9];
    ifstream inf;
    inf.open(HFfilename);
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
    long double* dedens = new long double[dim_rd];
    long double* findens = new long double[dim_rd];
    int* charge = new int[nroa];
    double* dumv = new double[nroa];
    double* coord = new double[3*nroa];
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
    inf.close();
    for(int i = llim; i < homo; i++){
	sprintf(filename,"%s%i.bcb",input_pref,i);
	ifstream inf;
	inf.open(filename);
	inf.read((char *) &nroa_rd, sizeof(int));
	inf.read((char *) xgrid_rd, sizeof(double)*3);
	inf.read((char *) nrop_rd, sizeof(int)*3);
	inf.read((char *) sgrid_rd, sizeof(double)*9);
	inf.read((char *) charge_rd, sizeof(int)*nroa);
	inf.read((char *) dumv_rd, sizeof(double)*nroa);
	inf.read((char *) coord_rd, sizeof(double)*3*nroa);
	check_sizes(nroa_rd,nroa,xgrid_rd,xgrid,nrop_rd,nrop,sgrid_rd,sgrid,charge_rd,charge,dumv_rd,dumv,coord_rd,coord);
	inf.read((char *) dens_rd, sizeof(long double)*dim_rd);
	for(int j = 0; j < dim_rd; j++) dedens[j] = dens_rd[j];
	inf.close();
	for(int v = (i+1); v <= homo; v++){
	    sprintf(filename,"%s%i.bcb",input_pref,v);
	    ifstream inf;
	    inf.open(filename);
	    inf.read((char *) &nroa_rd, sizeof(int));
	    inf.read((char *) xgrid_rd, sizeof(double)*3);
	    inf.read((char *) nrop_rd, sizeof(int)*3);
	    inf.read((char *) sgrid_rd, sizeof(double)*9);
	    inf.read((char *) charge_rd, sizeof(int)*nroa);
	    inf.read((char *) dumv_rd, sizeof(double)*nroa);
	    inf.read((char *) coord_rd, sizeof(double)*3*nroa);
	    check_sizes(nroa_rd,nroa,xgrid_rd,xgrid,nrop_rd,nrop,sgrid_rd,sgrid,charge_rd,charge,dumv_rd,dumv,coord_rd,coord);
	    inf.read((char *) dens_rd, sizeof(long double)*dim_rd);
	    inf.close();
	    for(int j = 0; j < dim_rd; j++) findens[j] = dedens[j] * dens_rd[j];
	    sprintf(outname,"%s_%i_%i.bds",output_pref,i,v);
	    ofstream outf;
	    outf.open(outname);
	    outf.write((char *) &nroa, sizeof(int));
	    outf.write((char *) xgrid, sizeof(double)*3);
	    outf.write((char *) nrop, sizeof(int)*3);
	    outf.write((char *) sgrid, sizeof(double)*9);
	    outf.write((char *) charge, sizeof(int)*nroa);
	    outf.write((char *) dumv, sizeof(double)*nroa);
	    outf.write((char *) coord, sizeof(double)*3*nroa);
	    outf.write((char *) findens, sizeof(long double)*dim_rd);
	    outf.close();
	}
    }
}

void check_sizes(int nroa_rd, int nroa, double* xgrid_rd, double* xgrid, int* nrop_rd, int* nrop, double* sgrid_rd, double* sgrid, int* charge_rd, int* charge, double* dumv_rd, double* dumv, double* coord_rd, double* coord){
	if(nroa_rd != nroa){
	    cerr << "Nr of atoms is not equal\n";
	    exit(1);
	}
	if(xgrid_rd[0] != xgrid[0] || xgrid_rd[1] != xgrid[1] || xgrid_rd[2] != xgrid[2]){
	    cerr << "Origin of grid not equal\n";
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
}

