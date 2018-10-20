#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

void read_line(ifstream* inf, char* line_buf, int max_count){                   // routine to read strings from a file (a line is saved into line_buf)
  int count = 0;
  char c;
  int end = 0;
  while(end == 0 && inf->get(c)){
    line_buf[count++] = c;
    if(c == '\n') end = 1;
    if(count == max_count){
      cerr << "Buffer overflow while parsing\n"; exit(1);
    }
  }
  line_buf[count-1] = 0;
}

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Need input-file output-prefix\n";
	exit(0);
    }
    char line_buf[64000];                                                        // string where a single line from the GAMESS punch file is saved
    char commentA[64000];
    char commentB[64000];
    int nroa;
    double* xgrid = new double[3];
    int* nrop = new int[3];
    double* sgrid = new double[9];
    ifstream inf(argv[1]);
    char numchar[14];                                                            // character string to store the respective MO vector elements
    numchar[13] = 0;                                                             // is first set to zero
    for(int i = 0; i < 2 ; i++){
	read_line(&inf, line_buf,64000);
	if(i == 0){
	    sprintf(commentA,"%s",line_buf);
	}else{
	    sprintf(commentB,"%s",line_buf);
	}
    }
    read_line(&inf, line_buf,64000);
    istringstream ist(line_buf);
    ist >> nroa;
    for(int i = 0; i < 3 ; i++) ist >> xgrid[i];
    for(int i = 0; i < 3; i++){
	read_line(&inf,line_buf,64000);
        istringstream ist(line_buf);
	ist >> nrop[i];
	for(int j = 0; j < 3; j++) ist >> sgrid[i*3+j];
    }
    int* charge = new int[nroa];
    double* dumv = new double[nroa];
    double* coord = new double[3*nroa];
    int dim = nrop[0] * nrop[1] * nrop[2];
    long double* dens = new long double[dim];
    for(int i = 0; i < nroa; i++){
	read_line(&inf,line_buf,64000);
        istringstream ist(line_buf);
	ist >> charge[i] >> dumv[i];
	for(int j = 0; j < 3; j++) ist >> coord[3*i+j];
    }
    int count = 0;
    while(!inf.eof()){                                       // read lines until " $END   " is reached into line_buf
	read_line(&inf, line_buf,64000);
	int line_count = 0;
	while(line_count < 6){
	    for(int a = 0; a < 13; a++) numchar[a] = line_buf[a+13*line_count];                            // read 5 values with a width of 15 characters into numchar
	    line_count++;
	    if(count < dim){
		dens[count] = strtod(numchar, NULL);                   // MO vectors are transformed char -> double and written to MOs
		count++;
	    }
	    if(count%nrop[2] == 0){
		line_count = 6; 
	    }
	}
    }
    ofstream outf;
    char dumc[128];
    sprintf(dumc,"%s.bcb",argv[2]);
    outf.open(dumc);
    outf.write((char *) &nroa, sizeof(int));
//    cout << "Number of atoms: " << nroa << "\n";
    outf.write((char *) xgrid, sizeof(double)*3);
//    cout << "Grid origin: " << xgrid[0] << " " << xgrid[1] << " " << xgrid[2] << "\n";
    outf.write((char *) nrop, sizeof(int)*3);
//    cout << "Nr of points: " << nrop[0] << " x " << nrop[1] << " x " << nrop[2] << " = " << dim << "\n";
    outf.write((char *) sgrid, sizeof(double)*9);
/*    cout << "Grid stepwidth: \n";
    for(int i = 0; i < 9; i++){
	cout << sgrid[i] << " ";
	if(i%3 == 2) cout << "\n";
    }*/
    outf.write((char *) charge, sizeof(int)*nroa);
/*    cout << "Charges: \n";
    for(int i = 0; i < nroa; i++){
	cout << charge[i] << " ";
    if(i%3 == 2) cout << "\n";
    }*/
    outf.write((char *) dumv, sizeof(double)*nroa);
    cout << "System data: \n";
    for(int i = 0; i < nroa; i++){
	cout << dumv[i] << " ";
	if(i%3 == 2) cout << "\n";
    }
    outf.write((char *) coord, sizeof(double)*3*nroa);
/*    cout << "Coordinates: \n";
    for(int i = 0; i < nroa; i++){
	cout << coord[i*nroa+0] << " " << coord[i*nroa+1] << " " << coord[i*nroa+2] << "\n";
    }*/
    outf.write((char *) dens, sizeof(long double)*dim);
/*    cout << "Density: \n";
    for(int i = 0; i < dim; i++){
	cout << dens[i] << " ";
	if(dim%6 == 5) cout << "\n";
    }*/
    outf.flush();
}

