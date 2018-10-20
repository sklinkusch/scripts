# include <iostream>
# include <fstream>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <sstream>
# include <math.h>

using namespace std;

void swap_line(int &froma, int &toa, double &coeffa, int &fromb, int &tob, double &coeffb){
    int tempfrom, tempto;
    double tempcoeff;
    tempfrom = froma;
    froma = fromb;
    fromb = tempfrom;
    tempto = toa;
    toa = tob;
    tob = tempto;
    tempcoeff = coeffa;
    coeffa = coeffb;
    coeffb = tempcoeff;
}

int main(int argc, char* argv[]){
    if(argc != 2){
	cerr << "Usage: ./sort_excitations <infile>\n";
	exit(0);
    }
    ifstream inf;
    inf.open(argv[1]);
    int nros;
    double gsenergy;
    inf >> gsenergy >> nros;
    double* excenergies = new double[nros];
    for(int i = 0; i < nros; i++) inf >> excenergies[i];
    int limit, state; 
    double dndum;
    double* upcoeff = new double[500];
    int* frommo = new int[500];
    int* tomo = new int[500];
    for(int i = 0; i < nros; i++){
	inf >> limit;
	if(limit > 0){
	for(int j = 0; j < limit; j++){
	   inf >> state >> frommo[j] >> tomo[j] >> upcoeff[j] >> dndum;
	   cout << state << " " << frommo[j] << " " << tomo[j] << " " << upcoeff[j] << " " << dndum << "\n";
	}
	for(int k = 0; k < (limit-2); k++){
	    if(k%2 == 0){
#pragma omp parallel for
		for(int j = 0; j < (limit/2); j++){
		    if(fabs(upcoeff[2*j]) < fabs(upcoeff[2*j+1])){
			    swap_line(frommo[2*j], tomo[2*j], upcoeff[2*j], frommo[2*j+1], tomo[2*j+1], upcoeff[2*j+1]);
		    }
		}
	    }else{
#pragma omp parallel for
		for(int j = 0; j < (nros/2-1); j++){
		    if(fabs(upcoeff[2*j+1]) < fabs(upcoeff[2*j+2])){
			swap_line(frommo[2*j+1], tomo[2*j+1], upcoeff[2*j+1], frommo[2*j+2], tomo[2*j+2], upcoeff[2*j+2]);
		    }
		}
	    }
	}
	}
	cout << "State: " << i << "\n";
	if(limit >= 3){
	    for(int j = 0; j < 3; j++) cout << frommo[j] << "  " << tomo[j] << "  " << upcoeff[j] << "\n";
	}else if(limit > 0){
	    for(int j = 0; j < limit; j++) cout << frommo[j] << "  " << tomo[j] << "  " << upcoeff[j] << "\n";
	}else{
	    cout << "0  0  1.00000\n";
	}
    }
	    
}

