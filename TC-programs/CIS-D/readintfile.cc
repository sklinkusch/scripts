# include <iostream>
# include <fstream>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdint.h>
# include <time.h>
# include <unistd.h>
# include <ctype.h>
# include <sys/time.h>
# include <sys/resource.h>
# include <sstream>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Usage: Need input-file output-prefix\n";
	exit(1);
    }
    char inpfile[512];
    char outfile[512];
    sprintf(inpfile, "%s", argv[1]);
    sprintf(outfile, "%s.out", argv[2]);
    uint32_t nonzero;
    uint32_t x_max;
    uint32_t y_max;

    ifstream inf(inpfile);
    ofstream outf;
    outf.open(outfile);
    inf.read((char *) &nonzero, sizeof(uint32_t));
    outf << "Number of nonzero elements: " << nonzero << "\n";
    outf.flush();
    inf.read((char *) &x_max, sizeof(uint32_t));
    inf.read((char *) &y_max, sizeof(uint32_t));
    outf << "Maximum values for x and y: " << x_max << " " << y_max << "\n";
    outf.flush();
    double* prec_vals = new double[nonzero];
    uint32_t* prec_cols = new uint32_t[nonzero];
    uint32_t* prec_rows = new uint32_t[x_max+1];
    inf.read((char *) prec_vals, sizeof(double)*nonzero);
    inf.read((char *) prec_cols, sizeof(uint32_t)*nonzero);
    inf.read((char *) prec_rows, sizeof(uint32_t)*(x_max+1));
    for(uint32_t i = 0; i < nonzero; i++){
       outf << prec_vals[i] << "\n";
       outf.flush();
    }
    outf << "...............................................\n";
    outf.flush();
    for(uint32_t i = 0; i < nonzero; i++){
       outf << prec_cols[i] << "\n";
       outf.flush();
    }
    outf << ":::::::::::::::::::::::::::::::::::::::::::::::\n";
    outf.flush();
    for(uint32_t i = 0; i <= x_max; i++){
	outf << prec_rows[i] << "\n";
	outf.flush();
    }
    outf.close();
}


