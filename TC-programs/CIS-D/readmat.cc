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

/******************************************************************************
 * readMAT                                                                    *
 * developped by Stefan Klinkusch, Free University of Berlin (2014)           *
 * read CI eigenenergies and the CI matrix from a binary GAMESS matrix file   *
 ******************************************************************************/

// Extern Functions
extern void read_matrix(char* matfile, int nros, double* cisvecs, int type);

int main(int argc, char* argv[]){
 if(argc != 4 && argc != 5){                                                                     // if #arguments is not 4 and not 5 -> emergency exit (wrong #parameters)
  cerr << "Usage: ./readmat #states matfile logfile (optional: type)\n";
  cerr << "type=1 for Fortran binary files, type=0 for C++\n";
  cerr << "Default: type = 0 (if none or a wrong one given)\n";
  exit(1);
 }
 int type;                                                                                       // type of the binfile, 0: c++, 1: fortran
 if(argc == 4){                                                                                  // if no type is explicitly given, then c++ is assumed
  type = 0;
 }else{
  type = atoi(argv[4]);                                                                          // otherwise type is read from the shell
 }
 if(type != 0 && type != 1) type = 0;                                                            // wrong types are set to zero
 char matfile[256];                                                                              // name of the binary file containing the matrix
 char outrow[1024];                                                                              // string for formatted output
 sprintf(matfile, "%s", argv[2]);                                                                // store name of matfile in matfile
 int nros = atoi(argv[1]);                                                                       // number of states
 int rows = nros / 10;                                                                           // number of complete rows (needed for formatted output)
 int rest = nros%10;                                                                             // number of extra elements (needed for formatted output)
 double* cisvecs = new double[nros*nros];                                                        // array for the CI matrix elements
 read_matrix(matfile, nros, cisvecs, type);                                                      // read the matrix elements using an external routine
 ofstream outf;                                                                                  // output stream for the logfile
 outf.open(argv[3]);                                                                             // output filename is read from the shell, output file is opened
 outf.precision(10);                                                                             // output precision is set to ten numbers
 outf << "Number of States: " << nros << "\n";                                                   // first the number of states is printed
 outf << "----------------------------------------------------\n";
 outf << "CI matrix: \n";                                                                        // now the elements of the matrix are printed using formatted output
 for(int x = 0; x < (nros-1); x++){
  for(int y = 0; y < rows; y++){
  sprintf(outrow, "% 11.10f % 11.10f % 11.10f % 11.10f % 11.10f % 11.10f % 11.10f % 11.10f % 11.10f % 11.10f", cisvecs[x*(nros-1)+10*y+0],  cisvecs[x*(nros-1)+10*y+1],  cisvecs[x*(nros-1)+10*y+2],  cisvecs[x*(nros-1)+10*y+3],  
	  cisvecs[x*(nros-1)+10*y+4],  cisvecs[x*(nros-1)+10*y+5],  cisvecs[x*(nros-1)+10*y+6],  cisvecs[x*(nros-1)+10*y+7],  cisvecs[x*(nros-1)+10*y+8], cisvecs[x*(nros-1)+10*y+9]);
  outf << outrow << "\n";
  }
  for(int z = 0; z < rest; z++){
      sprintf(outrow, "% 11.10f ", cisvecs[x*(nros-1)+10*rows+z]);
      outf << outrow;
  }
  outf << "\n";
 }

 outf.close();
}

