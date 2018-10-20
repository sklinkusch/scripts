#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex>
#include<sstream>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Usage: ./read_ecp <input ecp file> <output log>\n";
	exit(1);
    }
    ifstream inf(argv[1]);
    ofstream outf(argv[2]);
    int nros;
    inf.read((char *) &nros, sizeof(int));
    outf << "Number of states: " << nros << "\n";
    outf.flush();
    double* ens = new double[nros];
    inf.read((char *) ens, sizeof(double)*nros);
    for(int i = 0; i < nros; i++){
	outf << ens[i] << " ";
	if(i%20 == 19) outf << "\n";
	outf.flush();
    }
    outf << "\n";
    outf << "Energies ranging from " << ens[0] << " to " << ens[nros-1] << " Hartree\n";
    outf.flush();
    delete [] ens;
    double* dx = new double[nros*nros];
    inf.read((char *) dx, sizeof(double)*nros*nros);
    outf << "Dipole moments (x) from " << dx[0] << " to " << dx[nros*nros-1] << " ea_0\n";
    outf.flush();
    delete [] dx;
    double* dy = new double[nros*nros];
    inf.read((char *) dy, sizeof(double)*nros*nros);
    outf << "Dipole moments (y) from " << dy[0] << " to " << dy[nros*nros-1] << " ea_0\n";
    outf.flush();
    delete [] dy;
    double* dz = new double[nros*nros];
    inf.read((char *) dz, sizeof(double)*nros*nros);
    outf << "Dipole moments (z) from " << dz[0] << " to " << dz[nros*nros-1] << " ea_0\n";
    outf.flush();
    delete [] dz;
    double* ion = new double[nros];
    inf.read((char *) ion, sizeof(double)*nros);
    outf << "Ionization rates from " << ion[0] << " to " << ion[nros-1] << " hartree/hbar\n";
    outf.flush();
    delete [] ion;
    double* cpl = new double[nros*nros];
    inf.read((char *) cpl, sizeof(double)*nros*nros);
    for(int i = 0; i < nros; i++){
	for(int j = 0; j < nros; j++){
	    outf << i << " " << j << " " << cpl[i*nros+j] << "\n";
	}
    }
    outf.flush();
    delete [] cpl;
    double* kin = new double[nros];
    inf.read((char *) kin, sizeof(double)*nros);
    outf << "Kinetic energies from " << kin[0] << " to " << kin[nros-1] << " hartree\n";
    outf.flush();
    for(int i = 0; i < nros; i++){
	outf << i << " " << kin[i] << "\n";
	outf.flush();
    }
    outf.close();
}
