#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sstream>
#include <stdint.h>

using namespace std;

/********************************************************************************************************
 * CIS(D) Program Suite                                                                                 *
 * developped by Tillmann Klamroth, University of Potsdam (2004)                                        *
 * modified by Stefan Klinkusch, Free University of Berlin (2014)                                       *
 * reads GAMESS one- and two-electron integrals, GAMESS CIS matrices and eigenenergies, MO energies and *
 * vectors and calculates a doubles correction to the eigenenergies                                     *
 ********************************************************************************************************/

//Variables

int           nroao;                        // Number of atomic orbitals/basis functions
int           nroa;                         // Number of atoms
int           nroe;                         // Number of electrons
int           nros;                         // Number of electronic states
int           llim;                         // lower limit of MOs used for correlation
int           ulim;                         // upper limit of MOs used for correlation
int           umo;                          // Number of unoccupied (virtual) MOs
int           omo;                          // Number of occupied MOs
int           cis_size;                     // size of CIS matrix
double        dftEgs;                       // DFT ground state energy

double*             cistmpmat;              // temporary matrix for CI
char                intfile[256];           // MO Integrals file
char                dftwavfile[256];        // file containing MOs and MO energies (DFT)

//Extern Functions

extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);                                                   // diagonalizes a matrix
extern int calc_mo2el_ind_x(int i, int j, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index x
extern int calc_mo2el_ind_y(int k, int l, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index y
extern double get_precalc_ints_sd(int i, int j, int k, int l,
	                           int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, 
				   uint32_t* prec_rows);                                                                    // extract precalculated integrals safely
extern int* init_prex( int omo, int umo, int llim);
extern double get_precalc_ints_ovov(int i, int j, int k, int l,
	                             int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols,
				    uint32_t* prec_rows);                                                                   // extract precalculated integrals fast ovov type
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                                              // reads MO data from the wavfile
extern void read_moint_sizes(uint32_t& nonzero, uint32_t& xmax, uint32_t& ymax, ifstream* inf);
extern void read_moints(uint32_t nonzero, uint32_t xmax, uint32_t ymax, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, ifstream* inf);
extern int rem_com(char* filename, char* streamstring, int string_length);                                                  // routine to read the input file
extern void status(ofstream* outf);                                                                                         // write status (date, time, and host) to the logfile

//Functions

int main(int argc, char* argv[]){
 if(argc != 3 && argc != 4){
  cerr << "Usage: ./grimme-20 input-file output-prefix or ./grimme-20 input-file output-prefix 1 for degenerate states\n";                                                                                // error message if the number of input parameters is wrong
  exit(1);                                                                                                                  // emergency exit
 }
    
 const double evfactor = 27.211383;                                                                                         // conversion hartree <-> eV
 const double ca = 0.317;                                                                                                   // prefactor c_1 
 const double cb = 0.033;                                                                                                   // prefactor c_2
 const double cc = 12700000;                                                                                                // prefactor c_3
 const double threshold = 1.e-07;
 const double prefactor = -0.025;
 int typus;
 if(argc == 4){
   typus = atoi(argv[3]);
 }else{
   typus = 0;
 }
 if(typus != 1) typus = 0;

 char dumc[2048];                                                                                                           // summy string for the logfile name
 sprintf(dumc, "%s.log", argv[2]);                                                                                          // appends log to output-prefix as name for the logfile
 ofstream outf(dumc);                                                                                                       // opens logfile
status(&outf);                                                                                                              // write status (see above)
outf << ".................................................................................\n";
outf << "Reading input from: " << argv[1] << "\n";

int buff_length = 65536;
char* file_buff = new char[buff_length];                                                                                    // file buffer to read input file and to exclude comments marked using '#'
rem_com(argv[1], file_buff, buff_length);                                                                                   // read input file to write everything from the input file (excluding comments) to the file buffer
istringstream ist(file_buff);                                                                                               // stringstream to read from the file buffer

ist >> intfile >> dftwavfile;                                                                                               // read the name of the F08 file and the wavefunction files (HF and DFT)
ist >> dftEgs >> nroao >> nroa;                                                                                             // read the RHF ground state energy, #AOs, #atoms, bytes per integral, integrals per record, #records, #integrals
ist >> nroe >> nros >> llim >> ulim;                                                                                        // read the #electrons, #electronic states, the lower and the upper limit of orbitals used for correlation 

outf << ".................................................................................\n";                              // prints the read data to log file (in order to show that everything is working right)
outf << "Integral file : " << intfile << "\n";
//outf << "Wavefunction file (RHF): " << hfwavfile << "\n";
outf << "Wavefunction file (DFT): " << dftwavfile << "\n";
outf << "DFT ground state energy: " << dftEgs << " E_h (" << dftEgs*evfactor << " eV)\n";
outf << "Nr of basis functions: " << nroao << "\n";
outf << "Nr of atoms: " << nroa << "\n";
outf << "................................................................................\n";
outf.flush();
// read data for perturbation part
outf << "CIS-DFT correction routine is starting\n";
if(nroe < 0)  nroe *= -1;

if(nroe%2 != 0){                                                                                                            // emergency exit if an unrestricted system is treated
 cerr << "Number of electrons not even!\n";
 exit(7);
}
outf.flush();

double*             KSMOs = new double[nroao*nroao];
double*             KSMOens = new double[nroao];
read_wav_HF(dftwavfile, nroao, KSMOens, KSMOs);                                                                                // read Kohn-Sham MOs and MO energies from wavefunction file
outf << "Number of electrons: " << nroe << "\n";
outf << "Lower limit of MOs used for correlation: " << llim << "\n";
outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
outf << "Number of electronic states: " << nros << "\n";
outf.flush();
double*             dftmat  = new double[nros*nros];
double*             dftvecs = new double[nros*nros];
double*             dftens  = new double[nros];
cis_size = nros;

#pragma omp parallel for
for(int x = 0; x < nros; x++){
 for(int y = 0; y < nros; y++){
     dftmat[x*cis_size+y] = 0.;
 }
}

  outf << "Calculating upper triangle of CIS-matrix and CIS-DFT matrix\n";
  outf.flush();
long long int Cis_size = (long long int) cis_size;
omo = (nroe / 2) - llim;
umo = ulim + 1 - omo - llim;
double tempen;
double* tempdiffao = new double[nroao];
for(int i = 0; i < (nroao-1); i++){
    if((KSMOens[i+1] - KSMOens[i]) < threshold){
	tempen = 0.5*(KSMOens[i] + KSMOens[i+1]);
	KSMOens[i] = tempen;
	KSMOens[i+1] = tempen;
	for(int j = 0; j < nroao; j++) tempdiffao[j] = 1./(sqrt(2.)) * (KSMOs[i*nroao+j] - KSMOs[(i+1)*nroao+j]);
	for(int j = 0; j < nroao; j++){
	   KSMOs[i*nroao+j] += KSMOs[(i+1)*nroao+j];
	   KSMOs[i*nroao+j] *= 1./(sqrt(2.));
	   KSMOs[(i+1)*nroao+j] = tempdiffao[j];
	}
    }
}
delete [] tempdiffao;

ifstream intf(intfile);
uint32_t nonzero, x_max, y_max;
read_moint_sizes(nonzero, x_max, y_max, &intf);
double* prec_vals = new double[nonzero];
uint32_t* prec_cols = new uint32_t[nonzero];
uint32_t* prec_rows = new uint32_t[x_max+1];
read_moints(nonzero, x_max, y_max, prec_vals, prec_cols, prec_rows, &intf);
double KSMOi, KSMOf, Jar, Kas, Kar;
    //// CODE FOR NORMAL LOOPS
#pragma omp parallel for
  for(int x = 1 ; x < cis_size; x++){
   int i1 = (x-1)/umo+llim;                         // min(i1) = llim, max(i1) = HOMO, loop over occupied orbitals
   int f1 = (x-1)%umo+omo+llim;                     // min(f1) = llim + omo = LUMO, max(f1) = HUMO, loop over virtual orbitals
   for(int y = x ; y < cis_size; y++){
    int i2 = (y-1)/umo+llim;                        // min(i2) = i1 =(min)=> llim, max(i2) = HOMO, loop over (higher) occupied orbitals
    int f2 = (y-1)%umo+omo+llim;                    // min(f2) = i2 =(min)=> LUMO, max(f2) = HUMO, loop over (higher) virtual orbitals
    //## PRECALCULATED INTEGRAL --> biggest memory requierd
//    cismat[x*Cis_size+y] +=   -get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (oovv)
//    cismat[x*Cis_size+y] += 2.*get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (ovov)
//    dftmat[x*Cis_size+y] +=   -ca*get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
//    dftmat[x*Cis_size+y] += 2.*get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
    //## SEMI-DIRECT TRANSFORMATION
    //cismat[x*cis_size+y] +=   -mo2int_op(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //cismat[x*cis_size+y] += 2.*mo2int_op(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //## BRUT FORCE CALCULATIONS  --> lowest memory requierd
    //cismat[x*cis_size+y] +=   -calc_mo2int_bf(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    //cismat[x*cis_size+y] += 2.*calc_mo2int_bf(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    if(i1==i2 && f1==f2){
      if(typus == 1){
      for(int p = -1; p < 2; p+=2){
	  for(int q = -1; q < 2; q+=2){
//	    if((i1+p) < llim || f1+q > ulim) continue;
	    if(fabs(KSMOens[i1+p]-KSMOens[i1]) < threshold && fabs(KSMOens[f1+q]-KSMOens[f1]) < threshold){
	      //KSMOi = 0.5*(KSMOens[i1]+KSMOens[i1+p]);
	      KSMOi = KSMOens[i1];
	      //KSMOf = 0.5*(KSMOens[f1]+KSMOens[f1+q]);
	      KSMOf = KSMOens[f1];
	      //Jar = 0.25*(get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, i2+p, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows)
	      //+ get_precalc_ints_sd(i1, i2, f1+q, f2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, i2+p, f1+q, f2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kas = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kar = 0.25*(get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, f1, i2+p, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows)
	      + get_precalc_ints_sd(i1, f1+q, i2, f2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, f1+q, i2+p, i2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      dftmat[x*cis_size + y] += -ca * Jar;
	      dftmat[x*cis_size + y] += 2. * Kas;
	      dftmat[x*cis_size + y] += KSMOf - KSMOi;
	      dftmat[x*cis_size + y] += prefactor*KSMOi + cb * exp(-cc*pow(Kar,4.));
	      goto newloop;
	    }else if(fabs(KSMOens[i1+p]-KSMOens[i1]) < threshold && fabs(KSMOens[f1+q]-KSMOens[f1]) >= threshold){
	      if(q == 1){
	      //KSMOi = 0.5*(KSMOens[i1]+KSMOens[i1+p]);
	      KSMOi = KSMOens[i1];
	      KSMOf = KSMOens[f1];
	      //Jar = 0.5*(get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, i2+p, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kas = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kar = 0.5*(get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1+p, f1, i2+p, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      dftmat[x*cis_size + y] += -ca * Jar;
	      dftmat[x*cis_size + y] += 2. * Kas;
	      dftmat[x*cis_size + y] += KSMOf - KSMOi;
	      dftmat[x*cis_size + y] += prefactor*KSMOi + cb * exp(-cc*pow(Kar,4.));
	      goto newloop;
	      }else{
	      continue;
	      }
	    }else if(fabs(KSMOens[i1+p]-KSMOens[i1]) >= threshold && fabs(KSMOens[f1+q] - KSMOens[f1]) < threshold){
	      if(p == 1){
	      KSMOi = KSMOens[i1];
	      //KSMOf = 0.5*(KSMOens[f1]+KSMOens[f1+q]);
	      KSMOf = KSMOens[f1];
	      //Jar = 0.5*(get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1, i2, f1+q, f2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kas = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kar = 0.5*(get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows) + get_precalc_ints_sd(i1, f1+q, i2, f2+q, omo, umo, llim, prec_vals, prec_cols, prec_rows));
	      dftmat[x*cis_size + y] += -ca * Jar;
	      dftmat[x*cis_size + y] += 2. * Kas;
	      dftmat[x*cis_size + y] += KSMOf - KSMOi;
	      dftmat[x*cis_size + y] += prefactor*KSMOi + cb * exp(-cc*pow(Kar,4.));
	      goto newloop;
	      }else{
	      continue;
	      }
	    }else if(p == -1 && q == -1){
	      continue;
	    }else if(p == -1 && q == +1){
	      goto newp;
	    }else if(p == +1 && q == -1){
	      continue;
	    }else{
	      KSMOi = KSMOens[i1];
	      KSMOf = KSMOens[f1];
	      Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kar = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	      Kas = Kar;
	      dftmat[x*cis_size + y] += -ca * Jar;
	      dftmat[x*cis_size + y] += 2. * Kas;
	      dftmat[x*cis_size + y] += KSMOf - KSMOi;
	      dftmat[x*cis_size + y] += prefactor*KSMOi + cb * exp(-cc*pow(Kar,4.));
	      goto newloop;
	    }
	}
	newp:
	;
      }
      }else{
	KSMOi = KSMOens[i1];
	KSMOf = KSMOens[f1];
	Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	Kar = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	Kas = Kar;
	dftmat[x*cis_size + y] += -ca*Jar;
	dftmat[x*cis_size + y] += 2.*Kas;
	dftmat[x*cis_size + y] += KSMOf - KSMOi;
	dftmat[x*cis_size + y] += prefactor*KSMOi + cb * exp(-cc*pow(Kar,4.));
	goto newloop;
      }
      }else{
	Jar = get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	Kar = get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows);
	dftmat[x*cis_size+y] += -ca * Jar;
	dftmat[x*cis_size+y] += 2. * Kar;
	goto newloop;
   }
 newloop: 
 ;
  }
  }
#pragma omp parallel for
 for(int x = 0; x < cis_size; x++){
  for(int y = x; y < cis_size; y++){
   dftmat[y*cis_size+x] = dftmat[x*cis_size+y];
  }
 }
outf << "\n";
outf.flush();

diag_mat(cis_size, dftmat, dftens, dftvecs);
delete [] dftmat;
status(&outf);
char dumline[2048];
 outf << "State        Eexc(DFT) [E_h]         Eexc(DFT) [eV]   \n";
 outf << "------------------------------------------------------\n";
for(int x = 1; x < cis_size; x++){                                                                                          // prints out CIS(D) energies for all states
 sprintf(dumline, "%5d        % 9.6f             % 9.4f", x, dftens[x], dftens[x]*evfactor);
 outf << dumline << "\n";
 outf.flush();
}
status(&outf);
outf.flush();
outf.close();
}

