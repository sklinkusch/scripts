# include <iostream>
# include <fstream>
# include <math.h>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 5){
	cerr << "Usage: ./process_dftinfo <infile (ascii)> <#occ. mos> <#virt. mos> <outfile (binary)>\n";
	exit(1);
    }
    double gsenergy;
    int nros;
    int nroc;
    int omo = atoi(argv[2]);
    int umo = atoi(argv[3]);
//    int llim = atoi(argv[4]);
    int nroao = omo+umo;
    int nroctot = omo*umo+1;
    int state, idum, jdum;
    double updum, downdum;
    ifstream inf;
    inf.open(argv[1]);
    ofstream outf;
    outf.open(argv[4]);
    inf >> gsenergy;
    inf >> nros;
    double* exenergy = new double[nros];
    for(int x = 0; x < nros; x++) inf >> exenergy[x];
    double* upcoeffs = new double[nros*nroctot];
    double* dncoeffs = new double[nros*nroctot];
    for(int z = 0; z < (nros*nroctot); z++){
	upcoeffs[z] = 0.;
	dncoeffs[z] = 0.;
    }
    for(int x = 0; x < nros; x++){
	inf >> nroc;
	for(int y = 0; y < nroc; y++){
	    inf >> state >> idum >> jdum >> updum >> downdum;
	    if(state != x){
		cerr << "Wrong configuration in file " << argv[1] << "(" << state << "," << x << "," << idum << "," << jdum << "," << ")\n";
		exit(2);
	    }
	    upcoeffs[x*nroctot+(idum-1)*umo+(jdum-omo)] = updum;
	    dncoeffs[x*nroctot+(idum-1)*umo+(jdum-omo)] = downdum;
	}
    }
    upcoeffs[0] = 1.;
    dncoeffs[0] = 1.;
//    outf.write((char *) &gsenergy, sizeof(double));
//    outf.write((char *) &omo, sizeof(int));
//    outf.write((char *) &umo, sizeof(int));
//    outf.write((char *) &llim, sizeof(int));
    outf.write((char *) &nroao, sizeof(int));
    outf.write((char *) &nros, sizeof(int));
    outf.write((char *) &nroctot, sizeof(int));
    outf.write((char *) exenergy, sizeof(double)*nros);
    outf.write((char *) upcoeffs, sizeof(double)*nros*nroctot);
//    outf.write((char *) dncoeffs, sizeof(double)*nros*nroctot);
    outf.flush();
    outf.close();
    inf.close();
}

