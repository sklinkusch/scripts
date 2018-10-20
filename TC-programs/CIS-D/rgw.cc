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
#include <iomanip>

using namespace std;

/****************************************************************************************
 * rgw (Read GAMESS Wavefunction)                                                       *
 * developped by Tillmann Klamroth, University of Potsdam (2004)                        *
 * modified by Stefan Klinkusch, Free University of Berlin (2014)                       *
 * reads MO energies and MO vectors from a modified GAMESS punch file and writes it to  *
 * a binary .dat.hwf file                                                               *
 ****************************************************************************************/

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
    cerr << "Need nroao dat-file\n";
    exit(0);
  }
  int nroao = atoi(argv[1]);                                                   // number of AOs
  double* MOens = new double[nroao];                                           // array for MO energies
  double* MOs   = new double[nroao*nroao];                                     // array for MO vectors
  for(int x = 0; x < nroao; x++) MOens[x] = 0.;                                // MO energies are set to zero
  for(int x = 0; x < nroao*nroao; x++) MOs[x] = 0.;                            // MO vectors are set to zero
  double* tmpmat = new double[3*nroao];                                        // temporary array
  int maxcount = nroao / 6;                                                    // number of lines with MO energies
 
  ifstream inf;                                                                // input stream from file
  inf.open(argv[2]);
  ifstream iof;
  char numchar[16];                                                            // character string to store the respective MO vector elements
  numchar[15] = 0;                                                             // is first set to zero
  
  clog << "Reading from " << argv[2] << " Nroao is " << nroao << "\n";
  char line_buf[64000];                                                        // string where a single line from the GAMESS punch file is saved
  line_buf[1] = 0;                                                             // is first set to zero
  while(strcmp(" $EIGENVALUES", line_buf)!=0) read_line(&inf, line_buf,64000); // read until line with " $EIGENVALUES" is reached
    for(int i = 0; i <= maxcount; i++){                                        // loop over all lines with MO energies
     read_line(&inf, line_buf,64000);                                          // read each single line and write it to line_buf
     istringstream ist(line_buf);                                              // input stream from the buffer line_buf, able to read doubles
     for(int z = 0; z < 6; z++){                                               // loop over all elements in a line
      if((6*i+z) < nroao){                                                     // condition: only possible elements are read, no more than # AOs
       ist >> tmpmat[6*i+z];                                                   // read the MO energies into the temporary vector
      }
     }
    }
  for(int x = 0; x < nroao; x++){                                              // loop over all elements
   MOens[x] = tmpmat[x];						       // write MO energies into MOens
  }
  inf.close();
  iof.open(argv[2]);
  while(strcmp(" $VEC   ",line_buf)!=0)   read_line(&iof, line_buf,64000);     // read lines until " $VEC   " is reached
  int MOcount = 0;
  int AOcount = 0;
  while(strcmp(" $END   ",line_buf)!=0){                                       // read lines until " $END   " is reached into line_buf
    read_line(&iof, line_buf,64000);
    if(strcmp(" $END   ", line_buf)!=0){
      int line_count = 0;
      while(line_count < 5){
	for(int x = 0; x < 15; x++)
	  numchar[x] = line_buf[5+x+15*line_count];                            // read 5 values with a width of 15 characters into numchar
	line_count++;
	MOs[MOcount*nroao+AOcount] =  strtod(numchar, NULL);                   // MO vectors are transformed char -> double and written to MOs
	AOcount++;
	if(AOcount == nroao){
	  line_count = 5; 
	  AOcount = 0;
	  MOcount++;
	}
      }
    }
  }
  sprintf(line_buf,"%s.hwf",argv[2]);                                         // write the data to a .dat.hwf file (n binary format)
  ofstream outf(line_buf);
  outf.write((char *) &nroao, sizeof(int));
  outf.write((char *) MOens,  sizeof(double)*nroao);
  outf.write((char *) MOs,    sizeof(double)*nroao*nroao);

  outf.close();
  iof.close();

}
