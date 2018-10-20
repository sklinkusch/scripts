#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sstream>

using namespace std;

extern void read_bin(char* binfile, int nros, double& cisEgs, double* cisvals);

int main(int argc, char* argv[]){
    if(argc != 4){
     cerr << "Usage: ./readbin binfile nros outfile\n";
     exit(1);
    }

    char binfile[512];
    sprintf(binfile, "%s", argv[1]);
    int nros = atoi(argv[2]);
    double cisEgs = 0.;
    double* cisvals = new double[nros];
    read_bin(binfile, nros, cisEgs, cisvals);
    ofstream outf;
    outf.open(argv[3]);
    outf << nros << " States\n";
    outf << "GS Energy: " << cisEgs << "\n";
    outf << "Excitation Energies: " << "\n";
    for(int i = 0; i < nros; i++){
	outf << i << " " << cisvals[i] << "\n";
    }
    outf.close();
}

