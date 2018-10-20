#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

// Extern Functions
extern "C" void dgesv_(const int *N, const int* NHRS, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int* info);

int main(int argc, char* argv[]){

if(argc != 3 && argc != 4){
 cerr << "Need input-index, nr of timesteps, (input-file)!\n";
 exit(1);
}
int index = atoi(argv[1]);
int nrots = atoi(argv[2]);
double* xvals = new double[nrots];
double* yvals = new double[nrots];
if(index == 1){
ifstream inf;

inf.open(argv[3]);

for(int x = 0; x < nrots; x++){
 inf >> xvals[x] >> yvals[x];
}
inf.close();
}else{
for(int x = 0; x < nrots; x++){
 cin >> xvals[x] >> yvals[x];
}


double sumx4 = 0.;
double sumx3 = 0.;
double sumx2 = 0.;
double sumx1 = 0.;
double sumx2y = 0.;
double sumxy = 0.;
double sumy2 = 0.;
double sumy1 = 0.;
double* quadparam = new double[3];
double* A = new double[9];
double* b = new double[3];

for(int x = 0; x < nrots; x++){
 sumx4 += pow(xvals[x],4);
 sumx3 += pow(xvals[x],3);
 sumx2 += pow(xvals[x],2);
 sumx1 += xvals[x];
 sumx2y += pow(xvals[x],2)*yvals[x];
 sumxy += xvals[x]*yvals[x];
 sumy2 += pow(yvals[x],2);
 sumy1 += yvals[x];
}

A[0] = sumx4;
A[1] = sumx3;
A[2] = sumx2;
A[3] = sumx3;
A[4] = sumx2;
A[5] = sumx1;
A[6] = sumx2;
A[7] = sumx1;
A[8] = nrots;

b[0] = sumx2y;
b[1] = sumxy;
b[2] = sumy1;

// Performing square regression

int N = 3;
int NHRS = 1;
int lda = 3;
int ipiv[3];
int ldb = 3;
int info;

dgesv_(&N, &NHRS, A, &lda, ipiv, b, &ldb, &info);

for(int x = 0; x < 3; x++){
 quadparam[x] = b[x];
}

// Computing coefficient of determination

 double xav = sumx1/nrots;
 double yav = sumy1/nrots;
 
 double r = (sumxy - ((double) nrots*xav*yav))/sqrt((sumx2 - ((double) nrots*pow(xav,2)))*(sumy2 - ((double) nrots*pow(yav,2))));
 double r2 = pow(r,2);

cout << "Coefficient of determination: " << r2 << "\n";

// Performing linear regression

double* C = new double[2];

double delta = nrots*sumx2 - pow(sumx1,2);
C[0] = (nrots*sumxy - (sumx1*sumy1))/delta;
C[1] = (sumx2*sumy1 - (sumx1*sumxy))/delta;

cout.precision(6);
cout << "Koeffizient  quadratische Regression  lineare Regression\n";
cout << "a:           " << b[0] << "               -----------\n";
cout << "b:           " << b[1] << "               " << C[0] << "\n";
cout << "c:           " << b[2] << "               " << C[1] << "\n";

}

}

