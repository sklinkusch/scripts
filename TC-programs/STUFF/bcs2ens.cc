#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 4){
	cerr << "Usage ./bcs2ens <input-file> <offset in hartrees> <output-prefix>\n";
	exit(1);
    }
    char fin[128];
    sprintf(fin, "%s", argv[1]);
    ifstream inf;
    inf.open(fin);
    int nroao, nroe, llim, ulim;
    double offset = strtod(argv[2], NULL);
    char fout[128];
    sprintf(fout, "%s.dat", argv[3]);
    inf.read((char *) &nroao, sizeof(int));
    inf.read((char *) &nroe, sizeof(int));
    inf.read((char *) &llim, sizeof(int));
    inf.read((char *) &ulim, sizeof(int));
    int omo = nroe/2 - llim;
    int umo = ulim - llim - omo + 1;
    int cis_size = umo * omo + 1;
    double* cisens = new double[cis_size];
    double* vecs = new double[cis_size*cis_size];
    double* mus = new double[cis_size];
    double* corr = new double[cis_size];
    double* cisdens = new double[cis_size];
    inf.read((char *) cisens, sizeof(double)*cis_size);
    inf.read((char *) vecs, sizeof(double)*cis_size*cis_size);
    inf.read((char *) mus, sizeof(double)*cis_size);
    inf.read((char *) vecs, sizeof(double)*cis_size*cis_size);
    inf.read((char *) mus, sizeof(double)*cis_size);
    inf.read((char *) vecs, sizeof(double)*cis_size*cis_size);
    inf.read((char *) mus, sizeof(double)*cis_size);
    inf.read((char *) vecs, sizeof(double)*cis_size*cis_size);
    inf.read((char *) corr, sizeof(double)*cis_size);
    inf.close();
    delete [] mus;
    delete [] vecs;
    ofstream outf;
    outf.open(fout);
    for(int i = 0; i < cis_size; i++){
	cisdens[i] = cisens[i] + corr[i];
	cisdens[i] += offset;
	outf << i << "       " << cisdens[i] << "\n";
	outf.flush();
    }
    outf.close();
}

