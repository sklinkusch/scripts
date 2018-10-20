# include <iostream>
# include <fstream>
# include <math.h>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){

    if(argc != 3){
	cerr << "Usage: ./writemat asciifile dimension!\n";
	exit(1);
    }
    ifstream inf;
    inf.open(argv[1]);
    int nros = atoi(argv[2]);
    char binfile[512];
    sprintf(binfile, "%s.bin", argv[1]);
    ofstream outf(binfile);
    double* matrix = new double[nros*nros];
    for(int i = 0; i < nros; i++){
	for(int j = 0; j < nros; j++){
	    inf >> matrix[i*nros+j];
	    if(i == 0 && j < 10){
	     cout << matrix[i*nros+j] << " ";
	     if(j == 9) cout << "\n";
	    }
	}
    }
    outf.write((char *) matrix, sizeof(double)*nros*nros);
    outf.close();
}

