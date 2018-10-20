# include <iostream>
# include <fstream>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <sstream>
# include <math.h>

using namespace std;

void swap_line(int &statea, double &energya, double &popa, int &stateb, double &energyb, double &popb){
    int tempstate;
    double tempenergy, temppop;
    tempstate = statea;
    statea = stateb;
    stateb = tempstate;
    tempenergy = energya;
    energya = energyb;
    energyb = tempenergy;
    temppop = popa;
    popa = popb;
    popb = temppop;
}

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage: ./sort_populations <#states> <infile> <outfile>\n";
	exit(0);
    }

    int nros = atoi(argv[1]);
    int* state = new int[nros];
    double* energy   = new double[nros];
    double* pop = new double[nros];
    ifstream inf(argv[2]);
    ofstream outf;
    outf.open(argv[3]);
    for(int s = 0; s < nros; s++){
	inf >> state[s] >> energy[s] >> pop[s];
    }
    for(int k = 0; k < (nros-2); k++){
	if(k%2 == 0){
	    #pragma omp parallel for
	    for(int x = 0; x < (nros/2); x++){
		if(pop[2*x+1] > pop[2*x]){
		    swap_line(state[2*x], energy[2*x], pop[2*x], state[2*x+1], energy[2*x+1], pop[2*x+1]);
		}
	    }
	}else{
	    #pragma omp parallel for
	    for(int x = 0; x < (nros/2-1); x++){
		if(pop[2*x+2] > pop[2*x+1]){
		    swap_line(state[2*x+1], energy[2*x+1], pop[2*x+1], state[2*x+2], energy[2*x+2], pop[2*x+2]);
		}
	    }
	}
    }
    char dumc[1024];
    for(int s = 0; s < nros; s++){
	sprintf(dumc, "%4d   %10.6f    %11.10f", state[s], energy[s], pop[s]);
	outf << dumc << "\n";
	outf.flush();
    }
    outf.close();
}

