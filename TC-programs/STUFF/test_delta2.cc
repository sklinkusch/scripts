#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

using namespace std;

double deltafkt(double x);

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Need resonance frequency and laser frequency!\n";
	exit(0);
    }
    double omega = strtod(argv[1],NULL);
    double laser = strtod(argv[2],NULL);
    double diff  = omega - laser;
    double delta = deltafkt(diff);
    printf("%8.6f   %12.6f\n",diff,delta);
}
double deltafkt(double x){
    double delta = 0.;
    delta = (1/M_PI)*(1/pow(x,2.));
    return(delta);
}
