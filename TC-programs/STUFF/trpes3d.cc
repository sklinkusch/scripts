#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./trpes3d <inputfile> <nr of rows> <outputfile>\n";
	exit(1);
    }
    ifstream inf;
    ofstream outf;
    inf.open(argv[1]);
    int nror = atoi(argv[2]);
//    int nroc = 9;
    int nroc = 8;
    int nrox = nroc - 2;
    outf.open(argv[3]);
    double hartree;
    double evolt;
    double* delay = new double[nrox];
    delay[0] = 100.;
    delay[1] = 200.;
    delay[2] = 400.;
    delay[3] = 800.;
    delay[4] = 1600.;
    delay[5] = 2400.;
//    delay[5] = 3000.;
//    delay[6] = 6000.;
    double* signal = new double[nrox];
    for(int x = 0; x < nror; x++){
	inf >> hartree >> evolt;
	for(int y = 0; y < nrox; y++){
	    inf >> signal[y];
	    outf << hartree << " " << evolt << " " << delay[y] << " " << signal[y] << "\n";
	}
	outf << "\n";
    }
    inf.close();
    outf.close();
}

