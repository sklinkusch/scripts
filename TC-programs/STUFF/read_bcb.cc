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
    if(argc != 2){
	cerr << "Need input-file\n";
	exit(0);
    }
    int nroa;
    double* xgrid = new double[3];
    int* nrop = new int[3];
    double* sgrid = new double[9];
    ifstream inf(argv[1]);
    inf.read((char *) &nroa, sizeof(int));
    int* charge = new int[nroa];
    double* dumv = new double[nroa];
    double* coord = new double[3*nroa];
    cout << "Number of atoms: " << nroa << "\n";
    inf.read((char *) xgrid, sizeof(double)*3);
    cout << "Extra grid size: " << xgrid[0] << " " << xgrid[1] << " " << xgrid[2] << "\n";
    inf.read((char *) nrop, sizeof(int)*3);
    int dim = nrop[0] * nrop[1] * nrop[2];
    long double* dens = new long double[dim];
    cout << "Nr of points: " << nrop[0] << " x " << nrop[1] << " x " << nrop[2] << " = " << dim << "\n";
    inf.read((char *) sgrid, sizeof(double)*9);
    cout << "Grid stepwidth: \n";
    for(int i = 0; i < 9; i++){
	cout << sgrid[i] << " ";
	if(i%3 == 2) cout << "\n";
    }
    inf.read((char *) charge, sizeof(int)*nroa);
    cout << "Charges: \n";
    for(int i = 0; i < nroa; i++){
	cout << charge[i] << " ";
	if(i%3 == 2) cout << "\n";
    }
    inf.read((char *) dumv, sizeof(double)*nroa);
    cout << "System data: \n";
    for(int i = 0; i < nroa; i++){
	cout << dumv[i] << " ";
	if(i%3 == 2) cout << "\n";
    }
    inf.read((char *) coord, sizeof(double)*3*nroa);
    cout << "Coordinates: \n";
    for(int i = 0; i < nroa; i++){
	cout << coord[i*nroa+0] << " " << coord[i*nroa+1] << " " << coord[i*nroa+2] << "\n";
    }
    inf.read((char *) dens, sizeof(long double)*dim);
    cout << "Density: \n";
    for(int i = 0; i < dim; i++){
	cout << dens[i] << " ";
	if(dim%6 == 5) cout << "\n";
    }
}

