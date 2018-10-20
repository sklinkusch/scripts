#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./energyscheme <nr of states> <input file> <output-prefix>\n";
	exit(1);
    }

    int nros = atoi(argv[1]);
    ifstream inf;
    int curr_st;
    char dumc[512];
    double* cisens = new double[nros];
//    double* ionrat = new double[nros];
//    double pop;
    inf.open(argv[2]);
    for(int i = 0; i < nros; i++){
	inf >> curr_st >> cisens[i]; //>> ionrat[i] >> pop;
    }
    inf.close();
    ofstream outf;
    sprintf(dumc, "%s.dat", argv[3]);
    outf.open(dumc);
    outf << "1.0  ";
    for(int i = 0; i < nros; i++){
	outf << cisens[i] << "   ";
    }
    outf << "\n";
    outf.flush();
    outf << "2.0  ";
    for(int i = 0; i < nros; i++){
	outf << cisens[i] << "   ";
    }
    outf << "\n";
    outf.flush();
    outf.close();

/*    ofstream irxf;
    sprintf(dumc, "%s.irx", argv[3]);
    irxf.open(dumc);
    for(int i = 0; i < nros; i++){
	irxf << i << "   " << cisens[i] << "   " << ionrat[i] << "\n";
	irxf.flush();
    }
    irxf.close();*/
}

