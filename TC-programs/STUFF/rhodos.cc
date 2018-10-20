#include<iostream>
#include<fstream>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3 && argc != 4){
	cerr << "Usage ./rhodos <nr of frequencies> <infile> <outfile (optionally)>\n";
    }
    int modus = 0;
    if(argc == 4){
       modus = 1;
    }
    ofstream outf;
    int nrov = atoi(argv[1]);
    ifstream inf;
    inf.open(argv[2]);
    double* freq = new double[nrov];
    double* intens = new double[nrov];
    for(int v = 0; v < nrov; v++){
	inf >> freq[v] >> intens[v];
    }
    double width = freq[0];
    int nrost = nrov*1000;
    double extra_freq = 0.00023;
    double freq_max = freq[nrov-1] + extra_freq;
    double freq_step = freq_max / (nrost-1);
    double* freqs = new double[nrost];
    double* vals = new double[nrost];
    double sum1, sum2;

    for(int s = 0; s < nrost; s++){
	freqs[s] = s*freq_step;
	vals[s]  = 0.;
	for(int v = 0; v < nrov; v++){
	    sum1 = freqs[s]-freq[v];
	    sum2 = 0.5*width;
	    vals[s] += (1/M_PI)*((sum2)/(pow(sum1,2)+pow(sum2,2)));
	}
    }
    if(modus == 1){
	outf.open(argv[3]);
    }

    for(int s = 0; s < nrost; s++){
	if(modus == 0){
	    cout << freqs[s] << "    " << vals[s] << "\n";
	}else{
	    outf << freqs[s] << "    " << vals[s] << "\n";
	    outf.flush();
	}
    }
    if(modus == 1) outf.close();
}

