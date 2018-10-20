#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void swap_f(double* d, int A, int B);
void swap_d(int* d, int A, int B);

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./weighted-dipole <ecp file> <pop file> <target state>\n";
	exit(1);
    }

    ifstream inf(argv[1]);
    int nros;
    inf.read((char *) &nros, sizeof(int));
    double* ens = new double[nros];
    double* dx  = new double[nros*nros];
    double* dy  = new double[nros*nros];
    double* dz  = new double[nros*nros];
    double* ion = new double[nros];
    int dumint;
    double dumdub;
    inf.read((char *) ens, sizeof(double)*nros);
    inf.read((char *) dx, sizeof(double)*nros*nros);
    inf.read((char *) dy, sizeof(double)*nros*nros);
    inf.read((char *) dz, sizeof(double)*nros*nros);
    inf.read((char *) ion, sizeof(double)*nros);
    double popion;
    double* pop = new double[nros];
    int x = atoi(argv[3]);

    ifstream popf(argv[2]);
    for(int x = 0; x < nros; x++){
	popf >> dumint >> dumdub >> pop[x];
    }
    popf >> popion;

    double* dip = new double[nros*nros];
    // calculate total dipole moment
    for(int j = 0; j < nros; j++){
	for(int y = j; y < nros; y++){
	    dip[j*nros+y] = sqrt(pow(dx[j*nros+y],2.)+pow(dy[j*nros+y],2.)+pow(dz[j*nros+y],2.));
	    dip[y*nros+j] = dip[j*nros+y];
	}
    }
    int k_max = x-2;

    double* dumf = new double[x];
    int* dumd = new int[x];
#pragma omp parallel for
    for(int y = 0; y < x; y++){
	dumf[y] = pop[y]*dip[y*nros+x];
	dumd[y] = y;
    }

    for(int k = 0; k < k_max; k++){
	if(k%2 == 0){
#pragma omp parallel for
	    for(int i = 0; i < (x/2); i++){
		if(dumf[2*i] < dumf[2*i+1]){
		    swap_f(dumf,2*i,2*i+1);
		    swap_d(dumd,2*i,2*i+1);
		}
	    }
	}else{
#pragma omp parallel for
	    for(int i = 0; i < (x/2-1); i++){
		if(dumf[2*i+1] < dumf[2*i+2]){
		    swap_f(dumf,2*i+1,2*i+2);
		    swap_d(dumd,2*i+1,2*i+2);
		}
	    }
	}
    }

    for(int i = 0; i < 5; i++) printf("%4d :   % 12.10f\n",dumd[i],dumf[i]);
}

void swap_f(double* d, int A, int B){
    double tmp;
    tmp = d[B];
    d[B] = d[A];
    d[A] = tmp;
}

void swap_d(int* d, int A, int B){
    int tmp;
    tmp = d[B];
    d[B] = d[A];
    d[A] = tmp;
}

