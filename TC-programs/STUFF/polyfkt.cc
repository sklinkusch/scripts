#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

// Extern Functions
extern "C" void dgesv_(const int *N, const int* NHRS, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int* info);

int main(int argc, char* argv[]){

if(argc != 4 && argc != 5){
 cerr << "Need input-index, max order of polynom, nr of timesteps, (input-file)!\n";
 exit(1);
}
int index = atoi(argv[1]);
int order = atoi(argv[2]);
int nrots = atoi(argv[3]);
int ord   = order + 1;
double* xvals = new double[nrots];
double* yvals = new double[nrots];
if(index == 1){
ifstream inf;

inf.open(argv[4]);

for(int x = 0; x < nrots; x++){
 inf >> xvals[x] >> yvals[x];
}
inf.close();
}else{
for(int x = 0; x < nrots; x++){
 cin >> xvals[x] >> yvals[x];
}


/*double sumx6 = 0.;
double sumx5 = 0.;
double sumx4 = 0.;
double sumx3 = 0.;
double sumx2 = 0.;
double sumx1 = 0.;
double sumx3y = 0.;
double sumx2y = 0.;
double sumxy = 0.;
double sumy2 = 0.;
double sumy1 = 0.;*/
double* param = new double[ord];
double* A = new double[ord*ord];
double* b = new double[ord];

/*for(int x = 0; x < nrots; x++){
 sumx6 += pow(xvals[x],6);
 sumx5 += pow(xvals[x],5);
 sumx4 += pow(xvals[x],4);
 sumx3 += pow(xvals[x],3);
 sumx2 += pow(xvals[x],2);
 sumx1 += xvals[x];
 sumx3y += pow(xvals[x],3)*yvals[x];
 sumx2y += pow(xvals[x],2)*yvals[x];
 sumxy += xvals[x]*yvals[x];
 sumy2 += pow(yvals[x],2);
 sumy1 += yvals[x];
}*/

for(int x = 0; x < nrots; x++){
 for(int i = 0; i < (ord*ord); i++){
  A[i] = 0.;
 }
 for(int i = 0; i < ord; i++){
  b[i] = 0.;
 }
 for(int i = 0; i < ord; i++){
  for(int j = 0; j < ord; j++){
   A[i*ord+j] += pow(xvals[x],(2*order-i-j));
  }
  b[i] += pow(xvals[x],(order-i))*yvals[x];
 }
}


/*A[0] = sumx6;
A[1] = sumx5;
A[2] = sumx4;
A[3] = sumx3;
A[4] = sumx5;
A[5] = sumx4;
A[6] = sumx3;
A[7] = sumx2;
A[8] = sumx4;
A[9] = sumx3;
A[10] = sumx2;
A[11] = sumx1;
A[12] = sumx3;
A[13] = sumx2;
A[14] = sumx1;
A[15] = nrots;

b[0] = sumx3y;
b[1] = sumx2y;
b[2] = sumxy;
b[3] = sumy1;*/

// Performing square regression

int N = ord;
int NHRS = 1;
int lda = ord;
int ipiv[ord];
int ldb = ord;
int info;

dgesv_(&N, &NHRS, A, &lda, ipiv, b, &ldb, &info);

for(int x = 0; x < ord; x++){
 param[x] = b[x];
}

// Computing coefficient of determination

/* double xav = sumx1/nrots;
 double yav = sumy1/nrots;
 
 double r = (sumxy - ((double) nrots*xav*yav))/sqrt((sumx2 - ((double) nrots*pow(xav,2)))*(sumy2 - ((double) nrots*pow(yav,2))));
 double r2 = pow(r,2);

cout << "Coefficient of determination: " << r2 << "\n";*/

// Performing linear regression

//double* C = new double[2];

//double delta = nrots*sumx2 - pow(sumx1,2);
//C[0] = (nrots*sumxy - (sumx1*sumy1))/delta;
//C[1] = (sumx2*sumy1 - (sumx1*sumxy))/delta;

cout.precision(6);
//cout << "Koeffizient  quadratische Regression  lineare Regression\n";
//cout << "a:           " << b[0] << "               -----------\n";
//cout << "b:           " << b[1] << "               " << C[0] << "\n";
//cout << "c:           " << b[2] << "               " << C[1] << "\n";

cout << "coefficient:     " << order << "\n";
for(int i = order; i >=0; i--){
cout <<  "a_" << i << ":    " <<  b[i] << "\n";
}

}

}

