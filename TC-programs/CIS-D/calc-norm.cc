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
sprintf(outfile, "%s.out", argv[2]);
ofstream datf;
datf.open(outfile);
double k;
double norm;
for(int i = 0; i < size; i++){
    norm = 0;
    for(int j = 0; j < size; j++){
	vecf >> k;
	norm += k*k;
    }
    datf << norm << "\n";
    datf.flush();
}
vecf.close();
datf.close();
}

