#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <sstream>

using namespace std;

/******************************************************************************
 * readLOG                                                                    *
 * developped by Stefan Klinkusch, Free University of Berlin (2014)           *
 * reads # states, RHF energy, and CIS excitation energies from a GAMESS      *
 * .log file and writes the data to a binary .bin.dat file                    *
 ******************************************************************************/

void read_line(ifstream* inf, char* line_buf, int max_count){                   // routine to read strings from a file (a line is saved into line_buf)
 int count = 0;                                                                 // #characters is set to zero
 char c;                                                                        // the read character
 int end = 0;                                                                   // dummy integer to see the end of line
 while(end == 0 && inf->get(c)){                                                // loop as long as the end of line is not yet reached, reads each character c
  line_buf[count++] = c;                                                        // saves the current character to the line_buf array
  if(c == '\n') end = 1;                                                        // sets end to one if the end of line is reached
  if(count == max_count){                                                       // emergency exit if the line is too long
    cerr << "Buffer overflow while parsing\n"; exit(1);
  }
 }
 line_buf[count-1] = 0;
}


int main(int argc, char* argv[]){

 if(argc != 4){
     cerr << "Usage: ./readlog <GAMESS logfile> <# states> <output file>\n";
     exit(1);
 }
     char line_buf[64000];                                                       // line_buffer to save each read line as long as it is used
     char ev[128];                                                               // dummy string to save the excitation energy in eV
     char kcal[128];                                                             // dummy string to save the excitation energy in kcal/mol
     char cm[128];                                                               // dummy string to save the excitation wavenumber in cm^-1
     char nm[128];                                                               // dummy string to save the excitation wavelength in nm
     char binfile[512];                                                          // filename for binary outfile
     sprintf(binfile, "%s.bin.dat", argv[1]);                                    // appends .bin.dat to the input logfile name
     ofstream datf(binfile);                                                     // defines datf as the output stream to the new binfile
     double RHFen = 0.;                                                          // Final RHF Energy
     char dumc[1024];                                                            // dummy string for all nasty things to read from the input file
     char sym[128];                                                              // dummy string to save the symmetry of an excited state
     ifstream inf(argv[1]);                                                      // defines inf as the input stream from the GAMESS logfile
     int nros = atoi(argv[2]);                                                   // reads #states from the shell
     datf.write((char *) &nros, sizeof(int));                                    // writes #states to the binfile (binary format)
     double* cisvals = new double[nros];                                         // defines an array for the CIS excitation energies
     ofstream outf;                                                              // defines an output stream for a readLOG logfile
     cisvals[0] = 0.;                                                            // sets the CI excitation energy for the ground state to zero
     outf.open(argv[3]);                                                         // opens the readLOG logfile (filename from the shell)
     while(strcmp("          DENSITY CONVERGED", line_buf)!=0) read_line(&inf, line_buf,64000);  // reads each line until "DENSITY CONVERGED" is reached
     for(int i = 0; i < 5; i++) read_line(&inf, line_buf,64000);                 // reads four lines of junk data and one important line (the last line of the five ones)
     istringstream ist(line_buf);                                                // defines ist as a streamstring from the line_buf
     for(int i = 0; i < 4; i++) ist >> dumc;                                     // reads four words of junk and saves it in dumc
     ist >> RHFen;                                                               // reads the RHF energy
     outf << RHFen << "\n";                                                      // writes the RHF energy to the readLOG logfile
     datf.write((char *) &RHFen, sizeof(double));                                // writes the RHF energy to the binfile (binary format)
     while(strcmp("                    CI-SINGLES EXCITATION ENERGIES",line_buf)!=0) read_line(&inf,line_buf,64000); // reads each line until "CI-SINGLES EXCITATION ENERGIES" is reached
     for(int i = 0; i < 2; i++){                                                 // reads two lines of junk
	 read_line(&inf, line_buf,64000);
     }
     for(int i = 1; i < (nros); i++){                                            // read now (nros-1) lines with excitation energies
	 read_line(&inf, line_buf,64000);
         istringstream ist(line_buf);                                            // defines a stringstream reading from line_buf
	 ist >> sym >> cisvals[i] >> ev >>  kcal >> cm >> nm;                    // reads the symmetry, excitation energies (E_h, eV, kcal/mol, cm^-1, nm), only E_h is permanently saved
     }
     for(int i = 0; i < nros; i++){
	 outf << i << " " << cisvals[i] << "\n";                                 // writes CI excitation energies to the readLOG logfile
     }
     datf.write((char *) cisvals, sizeof(double)*nros);                          // writes CI excitation energies to the binfile (binary format)
     datf.close();                                                               // closes the binfile
     outf.close();                                                               // closes the readLOG logfile
}

