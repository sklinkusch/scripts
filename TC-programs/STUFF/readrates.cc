#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 2){
	cerr << "Usage: ./readsys <sys file>\n";
	exit(1);
    }
    char infile[256];
    sprintf(infile, "%s", argv[1]);
    int nbasis, nstates;
    ifstream inf;
    inf.open(infile);
    inf.read((char *) &nbasis, sizeof(int));
    inf.read((char *) &nstates, sizeof(int));
    double* configurations = new double[nbasis];
    double* ionization = new double[nbasis];
    double* dx = new double[nbasis*nbasis];
    double* dy = new double[nbasis*nbasis];
    double* dz = new double[nbasis*nbasis];
    double* rates = new double[nbasis*nbasis];
    double* dephasing = new double[nbasis*nbasis];
    double* energy = new double[nstates];
    double* eigenvectors = new double[nbasis*nstates];
    inf.read((char *) configurations, sizeof(double)*nbasis);
    inf.read((char *) ionization, sizeof(double)*nbasis);
    inf.read((char *) dx, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dy, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dz, sizeof(double)*nbasis*nbasis);
    inf.read((char *) rates, sizeof(double)*nbasis*nbasis);
    inf.read((char *) dephasing, sizeof(double)*nbasis*nbasis);
    inf.read((char *) energy, sizeof(double)*nstates);
    inf.read((char *) eigenvectors, sizeof(double)*nbasis*nstates);
    inf.close();
    cout << "Rates \n";
    cout << "Normal (0->223): " << rates[223] << " " << rates[66231] << "\n";
    cout << "B3LYP (0->186): " << rates[186] << " " << rates[55242] << "\n";
}

