#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

using namespace std;

double deltafkt(double x, double width);

int main(void){
    double increment = 1e-06;
    double matsubara = 2.37e-04;
    double max = 0.3;
    double min = 0.0;
//    double doub = 2./(M_PI*matsubara);
    double douba = 1./(pow(M_PI,2.)*pow(matsubara,2.));
    double x = min;
    double deltaa, deltab, diffen, prod;
    while(x <= max){
	deltaa = deltafkt(x,matsubara);
	diffen = fabs(x - matsubara);
	deltab = deltafkt(diffen,matsubara);
	prod = deltaa * deltab;
	printf("%8.6f   %12.6f   %12.6f   %20.6f   %20.6f\n",x,deltaa,deltab,prod,douba);
	x += increment;
    }
}
double deltafkt(double x, double width){
    double delta = 0.;
    delta = (1/M_PI)*((width/2.)/(pow(x,2.) + pow(width/2.,2.)));
    return(delta);
}
