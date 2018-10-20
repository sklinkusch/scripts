#include<iostream>
#include<fstream>
#include<string.h>
#include<stdlib.h>

//Program to read binary ef-files 
//Stefan Klinkusch, 2008

using namespace std;

int main(int argc, char* argv[]){
 if(argc != 3 && argc != 4){
  cerr << "Need index input-file (output-prefix)!\n";
  exit(1);
  }
int index = atoi(argv[1]);   // if 0 output -> terminal, if 1 output -> file
 if(index == 1){
  if(argc == 3){
   cerr << "Need output-prefix!\n";
   exit(1);
  }
char dumc[2048];
sprintf(dumc, "%s.fld", argv[3]);
ofstream outf(dumc);
}

int nrots;
double dt;
int* pfac = new int[3];

ifstream ief(argv[2]);
ief.read((char *) &nrots, sizeof(int));
ief.read((char *) &dt, sizeof(double));
ief.read((char *) pfac, 3*sizeof(int));

double* field_h = new double[3*nrots];
double* shape = new double[nrots];

ief.read((char *) field_h, 3*nrots*sizeof(double));
ief.read((char *) shape, nrots*sizeof(double));

if(index == 1){
 if(argc == 3){
  cerr << "Need output-prefix! \n";
  exit(1);
 }
char dumc[2048];
sprintf(dumc, "%s.fld", argv[3]);
ofstream outf(dumc);
for(int t = 0; t < nrots; t++){
 outf << t*dt << " " << field_h[0*nrots+t] << " " << field_h[1*nrots+t] << " " << field_h[2*nrots+t] << "\n";
}
outf.close();
}else{
for(int t = 0; t < nrots; t++){
 cout << t*dt << " " << field_h[0*nrots+t] << " " << field_h[1*nrots+t] << " " << field_h[2*nrots+t] << "\n";
 }
}
}

