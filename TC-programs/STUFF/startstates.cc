#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Usage: ./startstates <nr of timesteps> <ecp file> <pop file> <target state> <outfile>\n";
	exit(0);
    }

    int nrots = atoi(argv[1]);                 // number of time steps
    int nros;                                  // number of states
    int nrobs;                                 // nr of bound states
    int x, t, s;
    char ecpfile[256];                         // name of ecp file
    sprintf(ecpfile, "%s",argv[2]);
    ifstream ecpf;
    ecpf.open(ecpfile);
    ecpf.read((char *) &nros, sizeof(int));
    cout << "Number of total states: " << nros << "\n";
    double* ens = new double[nros];            // eigenenergies
    ecpf.read((char *) ens, sizeof(double)*nros);
    double* dx = new double[nros*nros];        // dipole moment, x component
    double* dy = new double[nros*nros];        // dipole moment, y component
    double* dz = new double[nros*nros];        // dipole moment, z component
    ecpf.read((char *) dx, sizeof(double)*nros*nros);
    ecpf.read((char *) dy, sizeof(double)*nros*nros);
    ecpf.read((char *) dz, sizeof(double)*nros*nros);
    double* ion = new double[nros];            // ionization rates
    ecpf.read((char *) ion, sizeof(double)*nros);
    ecpf.close();
    int count = 0;
    for(x = 0; x < nros; x++){
	if(ion[x] == 0.){
	    count++;
	}else{
	    break;
	}
    }
    cout << "Number of bound states: " << count << "\n";
    nrobs = count;
    double* musq = new double[nrobs];      // square absolute of total dipole moment
    int target = atoi(argv[4]);
    for(x = 0; x < nrobs; x++){
	musq[x] = pow(dx[x*nros+target],2.) + pow(dy[x*nros+target],2.) + pow(dz[x*nros+target],2.);
    }
    char popfile[256];
    sprintf(popfile, "%s", argv[3]);           // name of population file
    ifstream popf;
    popf.open(popfile);
    char outfile[256];
    sprintf(outfile, "%s",argv[5]);
    ofstream outf;
    outf.open(outfile);
    double time;
    double pop;
    double prob;
    char dumc[256];
    for(t = 0; t < nrots; t++){
	popf >> time;
	sprintf(dumc, "%12.5f  ",time);
	outf << dumc;
	for(s = 0; s <= nros; s++){
	    popf >> pop;
	    if(s < target && s < nrobs){
		prob = pop * musq[s];
		sprintf(dumc, "%10.6f  ",prob);
		outf << dumc;
	    }else if(s < target){
		continue;
	    }else if(s < nrobs){
		continue;
	    }else{
		continue;
	    }
	}
	outf << "\n";
	outf.flush();
    }
    outf.close();
    popf.close();
}

