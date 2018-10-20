#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Usage: ./genpipuls <frequency> <sigma> <amplitude> <nr of points> <t_final of propagation>\n";
	exit(1);
    }

    cout.precision(10);
    const long double attofs = 41.341373337;
    long double freq = strtod(argv[1],NULL);
    long double sigma = strtod(argv[2],NULL);
    long double amp = strtod(argv[3],NULL);
    int nrop = atoi(argv[4]);
    long double tf = strtod(argv[5],NULL);
    long double tfin = tf / attofs;
    long double dumx = 0.;
    long double dumy = 0.;
    long double incr = 2 * sigma / (nrop - 1);
    int z = 0;
    int nroptot = 0;

    if(tf > 2*sigma){
	nroptot = nrop + 10;
    }else if(tf == 2*sigma){
	nroptot = nrop;
    }else{
	cerr << "Pulse does not fit in propagation time\n";
	exit(2);
    }
    long double* xvals = new long double[nroptot];
    long double* yvals = new long double[nroptot];
    for(int x = 0; x < nrop; x++){
	dumx = x*incr;
	dumy = cos(freq*dumx)*pow(cos(M_PI*(dumx-sigma)/(2*sigma)),2.)*amp;
	xvals[x] = dumx;
	yvals[x] = dumy;
    }
    z += nrop;
    for(int x = nrop; x < nroptot; x++){
	dumx = x*incr;
	dumy = 0.000;
	if(dumx < tf){
	    xvals[x] = dumx;
	    yvals[x] = dumy;
	    z++;
	}
    }
    xvals[z] = tf;
    yvals[z] = 0.000;
    z++;
    cout << " " << z << "  " << tfin << "d0\n";
    for(int x = 0; x < z; x++){
	cout << " " << xvals[x] << "  " << yvals[x] << "\n";
    }
}

