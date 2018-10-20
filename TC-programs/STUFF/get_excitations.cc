# include <iostream>
# include <fstream>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <sstream>
# include <math.h>

using namespace std;

void read_line(ifstream* inf, char* line_buf, int max_count, long long int &linecount){
    int count = 0;
    char c;
    int end = 0;
    while(end == 0 && inf->get(c)){
	line_buf[count++] = c;
	if(c == '\n') end = 1;
	if(count == max_count){
	    cerr << "Buffer overflow while parsing\n";
	    exit(1);
	}
    }
    line_buf[count-1] = 0;
    linecount++;
}

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
    if(argc != 4){
	cerr << "Usage: ./get_excitations <#states> <infile> <outfile>\n";
	exit(0);
    }

    int nros = atoi(argv[1]);
    
    int* frommo = new int[nros];
    int* tomo   = new int[nros];
    double* mocoeff = new double[nros];
    char line_buf[64000];
    char dumdata[1024];
    ifstream inf(argv[2]);
    ofstream outf;
    outf.open(argv[3]);
    long long int linecount = 0;

    for(int s = 1; s <= nros; s++){
	while(strcmp("          FROM MO     TO MO", line_buf)!=0){
	    read_line(&inf, line_buf, 64000, linecount);
	    if(linecount%10000 == 0) cout << "Line " << linecount << " read\n";
	}
	read_line(&inf, line_buf, 64000, linecount);
	if(linecount%10000 == 0) cout << "Line " << linecount << " read\n";
	for(int i = 0; i < nros; i++){
	    read_line(&inf, line_buf, 64000, linecount);
	    if(linecount%10000 == 0) cout << "Line " << linecount << " read\n";
	    istringstream ist(line_buf);
	    ist >> frommo[i] >> tomo[i] >> mocoeff[i];
	}
	/*
	for(int k = 0; k < (nros-2); k++){
	    if(k%2 == 0){
		#pragma omp parallel for
		for(int x = 0; x < (nros/2); x++){
		    if(fabs(mocoeff[2*x+1]) > fabs(mocoeff[2*x])){
    			swap_line(frommo[2*x], tomo[2*x], mocoeff[2*x], frommo[2*x+1], tomo[2*x+1], mocoeff[2*x+1]);
		    }
		}
	    }else{
		#pragma omp parallel for
		for(int x = 0; x < (nros/2-1); x++){
		    if(fabs(mocoeff[2*x+2]) > fabs(mocoeff[2*x+1])){
			swap_line(frommo[2*x+1], tomo[2*x+1], mocoeff[2*x+1], frommo[2*x+2], tomo[2*x+2], mocoeff[2*x+2]);
		    }
		}
	    }
	}
	*/
	for(int i = 0; i < nros; i++){
	    sprintf(dumdata, "%4d     %2d   %3d    % 11.8f", s, frommo[i], tomo[i], mocoeff[i]);
	    outf << dumdata << "\n";
	    outf.flush();
	}
    }
    outf.close();
}

