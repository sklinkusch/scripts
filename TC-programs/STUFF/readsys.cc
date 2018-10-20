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
    cout << "Nbasis: " << nbasis << "\n";
    cout << "Nstates: " << nstates << "\n";
    cout << "Configurations: \n";
    for(int i = 0; i < nbasis; i++){
	cout << i << "  " << configurations[i] << "\n";
    }
    cout << "Ionization: \n";
    for(int i = 0; i < nbasis; i++){
	cout << i << "  " << ionization[i] << "\n";
    }
    cout << "Dipole (x,y,z) \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    cout << i << "\t" << j << "\t" << dx[i*nbasis+j] << "\t" << dy[i*nbasis+j] << "\t" << dz[i*nbasis+j] << "\n";
	}
    }
    cout << "Rates \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    cout << i << " " << j << " " << rates[i*nbasis+j] << "\n";
	}
    }
    cout << "Dephasing \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    cout << i << " " << j << " " << dephasing[i*nbasis+j] << "\n";
	}
    }
    cout << "Energy: \n";
    for(int i = 0; i < nstates; i++){
	cout << i << "  " << energy[i] << "\n";
    }
    cout << "Eigenvectors: \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nstates; j++){
	    cout << i << " " << j << " " << eigenvectors[i*nstates+j] << "\n";
	}
    }
}

