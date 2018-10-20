#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 3){
     cerr << "Usage: ./calc-norm <size> <vecfile>\n";
     exit(1);
 }
int size = atoi(argv[1]);
char vecfile[256];
sprintf(vecfile, "%s", argv[2]);
ifstream vecf;
vecf.open(vecfile);
char outfile[256];
sprintf(outfile, "%s.trnorm", argv[2]);
ofstream datf;
datf.open(outfile);
double* k = new double[size*size];
double* norm = new double[size];
for(int i = 0; i < size; i++){
    norm[i] = 0;
    for(int j = 0; j < size; j++){
	vecf >> k[i*size+j];
	norm[i] += k[j*size+i]*k[j*size+i];
    }
}
vecf.close();
double pref;
for(int i = 0; i < size; i++){
    pref = sqrt(norm[i]);
    for(int j = 0; j < size; j++){
	k[j*size+i] *= pref;
	datf << k[j*size+i] << " ";
    }
    datf << "\n";
    datf.flush();
}
datf.close();
}

