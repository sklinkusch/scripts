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
int           nrods;                        // Number of electronic states desired for CIS(D) correction
int           llim;                         // lower limit of MOs used for correlation
int           ulim;                         // upper limit of MOs used for correlation
int           umo;                          // Number of unoccupied (virtual) MOs
int           omo;                          // Number of occupied MOs
int           bpi;                          // Bytes per integral (from GAMESS logfile)
int           ipr;                          // Integrals per record (from GAMESS logfile)
int           nrorec;                       // Number of records (from GAMESS logfile)
int           cis_size;                     // size of CIS matrix
long long int nrofint;                      // Number of non-zero two-electron integrals
long long int sortcount[4];                 // boundaries for different permutation patterns
double        cisEgs;                       // RHF ground state energy

double*             intval;                 // two-electron integrals
unsigned short int* intnums;                // two-electron indices
double*             cistmpmat;              // temporary matrix for CI
char                f08file[256];           // GAMESS F08 file
char                wavfile[256];           // file containing MOs and MO energies
//char                binfile[256];           // file containing RHF energy and CI excitation energies
char                bcsfile[256];           // file to store binary data for time-dependent calculations

//Extern Functions

extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);                                                   // diagonalizes a matrix
//extern void diag_mat_mt(int argc, char** argv, int nroao, double* mat, double* vals, double* vecs);                        // diagonalizes a matrix
extern double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);                                               // calculates molecular orbital energy difference
extern inline double RO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);                                    // calculates prefactors for CSF Slater determinants
extern int pos_i_a(int i, int a, int omo, int umo, int llim);                                                               // calculates position of excitation i->a in CIS matrix
extern double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens);    
                                                                                                                            // calculates triples correction for i->a
extern double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim,
	                      double* cis_vec, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);                // calculates doubles correction for ij->ab
extern double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en,
	                  double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens);                      // calculates first sum in equation 14, CPL 219, 21-29 (1994)
extern double RO_calc_MP2(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, 
	                                       uint32_t* prec_rows, double* MOens);                                         // calculates MP2 energy correction to the ground state
extern double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
	                      double* cis_vec, double cis_en, double* t1, double* t2);                                      // calculates second sum in equation 14, CPL 219, 21-29 (1994)
extern void SL_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
	                    double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2,                     // calculates CIS(D) correction in one single loop for all states
			    int cis_size, ofstream* outf);
extern int calc_mo2el_ind_x(int i, int j, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index x
extern int calc_mo2el_ind_y(int k, int l, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index y
extern double get_precalc_ints_sd(int i, int j, int k, int l,
	                           int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, 
				   uint32_t* prec_rows);                                                                    // extract precalculated integrals safely
extern int* init_prex( int omo, int umo, int llim);
extern double get_precalc_ints_ovov(int i, int j, int k, int l,
	                             int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols,
				    uint32_t* prec_rows);                                                                   // extract precalculated integrals fast ovov type
extern "C" void dafrd_(double* V, long long int* LEN, long long int* RECN, char* FNAME);                                    // reads data from F10 files
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                                              // reads MO data from the wavfile
extern double* precalc_ints_sd(int omo, int umo, int llim, int nroao, double* MOs, long long int nrofint, 
	                       long long int* sortcount,  double* intval, unsigned short* intnums, ofstream* outf);         // routine to get the precalculated two-electron integrals
extern int rem_com(char* filename, char* streamstring, int string_length);                                                  // routine to read the input file
//extern void read_bin(char* binname, int nros, double& cisEgs, double* cisvals);                                             // routine to read by binfile
extern void status(ofstream* outf);                                                                                         // write status (date, time, and host) to the logfile
extern uint32_t  find_nonzeros(uint32_t dim, double* matrix);
extern void sparse_matrix(uint32_t nx, uint32_t ny, double* matrix, uint32_t nz, double* nonzero, 
	                  uint32_t* col_ind, uint32_t* row_ptr);
extern void SK_calc_cis_dc(int x, int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens, double* cisvec, double* cisen, int cis_size, double* t1);
extern void SK_calc_cis_tc(int x, int omo, int umo, int llim, double* cisvec, int cis_size, double *t2, double* twoel, double* MOdelta);
extern long long int calc_mapvar(int i, int j, int a, int b, int omo, int umo, int llim);
extern void calc_twoel_MOdelta(int omo, int umo, int llim, double* twoel, double* MOdelta, double* MOens, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);

//Functions

void read_integrals_12(void);                                                                                               // reads two-electron integrals from F08 file with 12 bytes per integral
void read_integrals_16(void);                                                                                               // reads two-electron integrals from F08 file with 16 bytes per integral
void resort_integrals(ofstream *outf);                                                                                      // routine to perform permutations of the two-electron integrals

int main(int argc, char* argv[]){
 if(argc != 3){
  cerr << "Need input-file output-prefix\n";                                                                                // error message if the number of input parameters is wrong
  exit(1);                                                                                                                  // emergency exit
 }
    
 const double evfactor = 27.211383;                                                                                         // conversion hartree <-> eV

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

ist >> f08file >> wavfile;                                                                                                  // read the name of the F08 file and the wavefunction file
ist >> cisEgs >> nroao >> nroa >> bpi >> ipr >> nrorec >> nrofint;                                                          // read the RHF ground state energy, #AOs, #atoms, bytes per integral, integrals per record, #records, #integrals
ist >> nroe >> nros >> nrods >> llim >> ulim;                                                                               // read the #electrons, #electronic states, #desired electronic states the lower and the upper limit of orbitals used for correlation 

outf << ".................................................................................\n";                              // prints the read data to log file (in order to show that everything is working right)
outf << "F08 file : " << f08file << "\n";
outf << "Wavefunction file: " << wavfile << "\n";
outf << "Hartree-Fock ground state energy: " << cisEgs << " E_h (" << cisEgs*evfactor << " eV)\n";
outf << "Nr of basis functions: " << nroao << "\n";
outf << "Nr of atoms: " << nroa << "\n";
outf.flush();
outf << "................................................................................\n";
outf << "Two-electron data: \n";
outf << "Bytes per integral: " << bpi << "\n";
outf << "Integrals per record: " << ipr << "\n";
outf << "Nr of records: " << nrorec << "\n";
outf << "Nr of two-electron integrals: " << nrofint << "\n";
outf << "................................................................................\n";
outf.flush();

if(bpi != 12 && bpi != 16){                                                                                                 // error message if neither 12 nor 16 bytes per integral, emergency exit
 outf << bpi << " bytes per integral not implemented yet\n";
 outf.flush();
 exit(2);
}

outf << "Allocating memory\n";

outf << "Need " << nrofint*8 << " bytes for two electron indices, and the same for two electron values.\n";

intval = new double[nrofint];
intnums = new unsigned short[nrofint*4];

outf << "Done \n";
outf << ".................................................................................\n";


outf.flush();
outf << "Reading 2el integrals\n";                                                                                          // read two-electron integrals from F08 file
if(bpi == 12) read_integrals_12();
if(bpi == 16) read_integrals_16();
outf << "Done \n";
outf.flush();
outf << ".................................................................................\n";
resort_integrals(&outf);                                                                                                    // final sorting of integrals
// read data for perturbation part
outf << "CIS(D) correction routine is starting\n";
if(nroe < 0)  nroe *= -1;

if(nroe%2 != 0){                                                                                                            // emergency exit if an unrestricted system is treated
 cerr << "Number of electrons (" << nroe << ") not even!\n";
 exit(7);
}
outf.flush();

double*             MOs = new double[nroao*nroao];                                                                          // molecular orbitals
double*             MOens = new double[nroao];                                                                              // MO energies
read_wav_HF(wavfile, nroao, MOens, MOs);                                                                                    // read MOs and MO energies from wavefunction file
outf << "Number of electrons: " << nroe << "\n";
outf << "Lower limit of MOs used for correlation: " << llim << "\n";
outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
outf << "Number of electronic states: " << nros << "\n";
outf << "Number of desired electronic states: " << nrods << "\n";
outf.flush();
double*             cismat  = new double[nros*nros];                                                                        // CI Matrix
double*             cisvecs = new double[nros*nros];                                                                        // CI Eigenvectors
double*             cisens  = new double[nros];
cis_size = nros;

#pragma omp parallel for
for(int x = 0; x < nros; x++){
 for(int y = 0; y < nros; y++){
     cismat[x*cis_size+y] = 0.;                                                                                             // sets all CI matrix elements to zero
 }
}

  outf << "Calculating upper triangle of CIS-matrix\n";
  outf.flush();
long long int Cis_size = (long long int) cis_size;
omo = (nroe / 2) - llim;
umo = ulim + 1 - omo - llim;
double* prec_ints = precalc_ints_sd(omo, umo, llim, nroao, MOs, nrofint,  sortcount, intval, intnums, &outf);               // writes the precalculated integrals into the array prec_ints
delete [] MOs;
delete [] intval;
delete [] intnums;
uint32_t x_max = calc_mo2el_ind_x(omo-1, omo+umo, omo, umo, 0); 
uint32_t y_max = calc_mo2el_ind_y(omo+umo,omo+umo, omo,  umo, 0); 
uint32_t nonzero = find_nonzeros(x_max*y_max, prec_ints);
double* prec_vals = new double[nonzero];
uint32_t* prec_cols = new uint32_t[nonzero];
uint32_t* prec_rows = new uint32_t[x_max+1];
sparse_matrix(x_max, y_max, prec_ints, nonzero, prec_vals, prec_cols, prec_rows);
delete [] prec_ints;
char dumchar[512];
sprintf(dumchar, "%s.moint", argv[2]);
ofstream datf(dumchar);
datf.write((char *) &nonzero, sizeof(uint32_t));
datf.write((char *) &x_max, sizeof(uint32_t));
datf.write((char *) &y_max, sizeof(uint32_t));
datf.write((char *) prec_vals, sizeof(double)*nonzero);
datf.write((char *) prec_cols, sizeof(uint32_t)*nonzero);
datf.write((char *) prec_rows, sizeof(uint32_t)*(x_max+1));
datf.flush();


    //// CODE FOR NORMAL LOOPS
#pragma omp parallel for
  for(int x = 1 ; x < cis_size; x++){
   int i1 = (x-1)/umo+llim;                         // min(i1) = llim, max(i1) = HOMO, loop over occupied orbitals
   int f1 = (x-1)%umo+omo+llim;                     // min(f1) = llim + omo = LUMO, max(f1) = HUMO, loop over virtual orbitals
   for(int y = x ; y < cis_size; y++){
    int i2 = (y-1)/umo+llim;                        // min(i2) = i1 =(min)=> llim, max(i2) = HOMO, loop over (higher) occupied orbitals
    int f2 = (y-1)%umo+omo+llim;                    // min(f2) = i2 =(min)=> LUMO, max(f2) = HUMO, loop over (higher) virtual orbitals
    //## PRECALCULATED INTEGRAL --> biggest memory requierd
    cismat[x*Cis_size+y] +=   -get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (oovv)
    cismat[x*Cis_size+y] += 2.*get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (ovov)
    //## SEMI-DIRECT TRANSFORMATION
    //cismat[x*cis_size+y] +=   -mo2int_op(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //cismat[x*cis_size+y] += 2.*mo2int_op(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //## BRUT FORCE CALCULATIONS  --> lowest memory requierd
    //cismat[x*cis_size+y] +=   -calc_mo2int_bf(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    //cismat[x*cis_size+y] += 2.*calc_mo2int_bf(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    if(i1==i2 && f1==f2)  cismat[x*cis_size+y] += MOens[f1]-MOens[i1];
   }
//   outf  << x << "\t";
//   if(x%10==0) outf << "\n";
//   outf.flush();
  }
#pragma omp parallel for
 for(int x = 0; x < cis_size; x++){
  for(int y = x; y < cis_size; y++){
   cismat[y*cis_size+x] = cismat[x*cis_size+y];                                                                             // copys the upper right triangle to the lower left triangle (symmetric matrix)
  }
 }
outf << "\n";
outf.flush();
// Perturbation theory part -> only implemented for precalculated integrals so far
outf << "Calculating doubles corrections\n";
double* cisd_corr = new double[nros];                                                                                       // CIS(D) correction
for(int x = 0; x < cis_size; x++) cisd_corr[x] = 0.;                                                                        // sets initially all correction values to zero
double E2 = RO_calc_MP2(omo, umo, llim, prec_vals, prec_cols, prec_rows, MOens);                                                                    // calculates ground state MP2 energy
outf << "Ground state MP2 energy (RO-OPS): T2 " << E2 << " Hartree (" << E2*evfactor << " eV)  TOT " << E2+cisEgs << " Hartree (" << (E2+cisEgs)*evfactor << " eV)\n";
outf << "...............................................................................\n";
status(&outf);
outf.flush();
double* C_t1 = new double[nros];                                                                                            // first term of the CIS(D) correction (double corrections)
double* C_t2 = new double[nros];                                                                                            // second term of the CIS(D) correction (triple corrections)
for(int i = 0; i < nros; i++){
    C_t1[i] = 0.;
    C_t2[i] = 0.;
}
//SINGLE LOOP CODE

// Output of non-corrected CIS energies
diag_mat(cis_size, cismat, cisens, cisvecs);
delete [] cismat;
double t1, t2;
int count_doub = 0;
int count_det = 0;
//SL_calc_cis_dc(omo, umo, llim, prec_vals, prec_cols, prec_rows, MOens, cisvecs, cisens, cistmpmat, C_t1, C_t2, cis_size, &outf);                  // calculates terms for CIS(D) correction for each state
#pragma omp parallel for
for(int x = 1; x <= nrods; x++){
    SK_calc_cis_dc(x, omo, umo, llim, prec_vals, prec_cols, prec_rows, MOens, cisvecs, cisens, cis_size, &t1);
    C_t1[x] = t1;
    outf << x << "\t";
    count_doub++;
    if(count_doub%10 == 0){
	outf << count_doub << "\n";
	outf.flush();
    }
}
outf << "\n";
outf << "............................................................................\n";
status(&outf);
outf << "Calculating triple corrections\n";
outf.flush();
long long int dimens = calc_mapvar(llim+omo-1, llim+omo-1, llim+omo+umo-1, llim+omo+umo, omo, umo, llim);
double* twoel = new double[dimens];
double* MOdelta = new double[dimens];
calc_twoel_MOdelta(omo, umo, llim, twoel, MOdelta, MOens, prec_vals, prec_cols, prec_rows); 
delete [] prec_vals;
delete [] prec_cols; 
delete [] prec_rows;
delete [] MOens;
for(int x = 1; x <= nrods; x++){
    SK_calc_cis_tc(x, omo, umo, llim, cisvecs, cis_size, &t2, twoel, MOdelta);
    C_t2[x] = t2;
    outf << x << "\t";
    count_det++;
    if(count_det%10 == 0){
    outf << count_det << "\n";
	outf.flush();
    }
}
outf << "\n";
status(&outf);
outf.flush();
#pragma omp parallel for
for(int x = 0; x < cis_size; x++) cisd_corr[x] = C_t1[x] + C_t2[x];                                                         // calculates total CIS(D) correction for each state
for(int x = 1; x < cis_size; x++){                                                                                          // prints out CIS(D) energies for all states
 outf << "State " << x << " : \n";
 outf << "doubles correction is " << C_t1[x] << " E_h (" << C_t1[x]*evfactor << " eV)\n";
 outf << "triples correction is " << C_t2[x] << " E_h (" << C_t2[x]*evfactor << " eV)\n";
 outf << "Total corr " << cisd_corr[x] << " E_h (" << cisd_corr[x] << " eV), excit.en. " << cisens[x]+cisd_corr[x] << " E_h (" << (cisens[x]+cisd_corr[x])*evfactor << "eV) from " << cisens[x] << " E_h (" << cisens[x]*evfactor << "eV)\n";
 outf << ".............................................................................\n";
}
status(&outf);
outf.flush();
outf.close();
}

// internal functions

void read_integrals_12(void){
 char dumchar[256];
 unsigned char num_buf[15000*4];   //(!!!)
 double val_buf[15000];

 long long int count = 0;

 ifstream inf(f08file);                                                                                                     // defines input stream and opens F08 file
 for(int x = 0; x < nrorec; x++){
  inf.read(dumchar, 12);                                                                                                    // reads 12 bytes of junk data at the beginning
  inf.read((char *) num_buf,ipr*4);                                                                                         // read characters
  inf.read((char *) val_buf,ipr*sizeof(double));                                                                            // read floats of double precision
  int buffcount = 0;
  int offset = ipr;
  if(count + ipr >  nrofint) offset = nrofint - count;
   for(long long int y = count; y < count+offset; y++)                                                                      // loop reading the two-electron integrals
    intval[y] = val_buf[buffcount++];
  buffcount = 0;
  //Bloody high low word kacke !!! (comment by Tillmann Klamroth, left here for historical reasons)
  for(long long int y = count; y < count+offset; y+=2){                                                                     // loop performing a permutation of the two-electron indices (always 8 indices, taking the second quartet before the first quartet)
   intnums[y*4+4] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+5] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+6] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+7] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+0] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+1] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+2] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+3] = (unsigned short ) num_buf[buffcount++]-1;
  }
  count += offset;
  inf.read( dumchar , 4);                                                                                                   // reads 4 bytes of junk data at the end
 }
}

void read_integrals_16(void){
 char dumchar[256];
 unsigned short num_buf[15000*4];   //(!!!)
 double val_buf[15000];

 long long int count = 0;

 ifstream inf(f08file);
 for(int x = 0; x < nrorec; x++){
  inf.read(dumchar, 12);                                                                                                    // reads 12 bytes of junk data at the beginning 
  inf.read((char *) num_buf,ipr*8);                                                                                         // read charakters
  inf.read((char *) val_buf,ipr*sizeof(double));                                                                            // read floats of double precision
  int buffcount = 0;
  int offset = ipr;
  if(count + ipr >  nrofint) offset = nrofint - count;
  for(long long int y = count; y < count+offset; y++)                                                                       // loop reading the two-electron integrals
   intval[y] = val_buf[buffcount++];
  buffcount = 0;
  //Bloody high low word kacke !!! (comment by Tillmann Klamroth, left here for historical reasons)
  for(long long int y = count; y < count+offset; y++){                                                                      // loop performing a permutation of the two-electron indices (always 4 indices, taking the second pair before the first pair)
   intnums[y*4+2] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+3] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+0] = (unsigned short ) num_buf[buffcount++]-1;
   intnums[y*4+1] = (unsigned short ) num_buf[buffcount++]-1;
  }
  count += offset;
  inf.read( dumchar , 4);                                                                                                   // reads 4 bytes of junk data at the end 
 }
}

inline void swap_ints(long long int a, long long int b, double* intvals, unsigned short* intnums){
 //temp = a
 unsigned short ta = intnums[a*4+0];
 unsigned short tb = intnums[a*4+1];
 unsigned short tc = intnums[a*4+2];
 unsigned short td = intnums[a*4+3];
 double tval = intvals[a];
 
 //a = b
 intnums[a*4+0] =   intnums[b*4+0];
 intnums[a*4+1] =   intnums[b*4+1];
 intnums[a*4+2] =   intnums[b*4+2];
 intnums[a*4+3] =   intnums[b*4+3];
 intvals[a]     =   intvals[b]    ;
 
 //b = temp
 intnums[b*4+0] =    ta;
 intnums[b*4+1] =    tb;
 intnums[b*4+2] =    tc;
 intnums[b*4+3] =    td;
 intvals[b]     =    tval;
}

void resort_integrals(ofstream *outf){
 // RESORT INTEGRALS

 // STEP 1: BRING to basis types
 for(long long int x = 0; x < nrofint; x++){
  unsigned short a,b,c,d;
  //INPUT AUCH  CHEMIKER ?? (comment by Tillmann Klamroth, left here for historical reasons)
  a = intnums[x*4+0];  
  b = intnums[x*4+1];  
  c = intnums[x*4+2];  
  d = intnums[x*4+3];  

  //TYP IIb: reorder to IIa  <ab|cd> = <aa|ca> -> <aa|ac>
  if(a==b && b==d && c!=d){
   intnums[x*4+0] = a;
   intnums[x*4+1] = a;
   intnums[x*4+2] = a; 
   intnums[x*4+3] = c;
  }

  //TYP IIc: reorder to IIa <ab|cd> = <ab|aa> -> <aa|ab>
  if(a!=b && a==c && c==d){
   intnums[x*4+0] = a;  //a
   intnums[x*4+1] = a;  //b
   intnums[x*4+2] = a;  //c
   intnums[x*4+3] = b;  //d
  }

  //TYP IId: reorder to IIa <ab|cd> = <ab|bb> -> <bb|ba>
  if(a!=b && b==c && c==d){
   intnums[x*4+0] = b;  //a
   intnums[x*4+1] = b;  //b
   intnums[x*4+2] = b;  //c
   intnums[x*4+3] = a;  //d
  }

  //TYP IIIc: reorder to IIIb <ab|cd> (= <ab|ba>) -> <ab|dc>
  if(a==d && b==c && a!=b){
   intnums[x*4+0] = a;  //a
   intnums[x*4+1] = b;  //b
   intnums[x*4+2] = d;  //c
   intnums[x*4+3] = c;  //d      
  }

  //TYP IVc:  reorder IVb <ab|cd> (= <ab|bd>) -> <ba|cd>
  if(b==c && b!=a && b!=d && a!=d){
   intnums[x*4+0] = b;  //a
   intnums[x*4+1] = a;  //b
   intnums[x*4+2] = c;  //c
   intnums[x*4+3] = d;  //d      
  }

  //TYP IVd:  reorder to IVa <ab|cd> = <ab|cc> -> <cc|ab>
  if(c==d && c!=a && c!=b && a!=b){
   intnums[x*4+0] = c;  //a
   intnums[x*4+1] = c;  //b
   intnums[x*4+2] = a;  //c
   intnums[x*4+3] = b;  //d      
  }

  //TYP IVe:  reorder IVb <ab|cd> (= <ab|ca>) -> <ab|dc>
  if(a==d && a!=b && a!=c && b!=c){
   intnums[x*4+0] = a;  //a
   intnums[x*4+1] = b;  //b
   intnums[x*4+2] = d;  //c
   intnums[x*4+3] = c;  //d      
  }

  //TYP IVf: perm_all,  reorder IVb <ab|cd> (= <ab|cb>) -> <ba|dc>
  if(b==d && b!=a && b!=c && a!=c){
   intnums[x*4+0] = b;  //a
   intnums[x*4+1] = a;  //b
   intnums[x*4+2] = d;  //c
   intnums[x*4+3] = c;  //d      
  }
 }

 //STEP2: BRING PERM1 to front
 long long int sorted = 0;
 for(long long int x = 0; x < nrofint; x++){
  unsigned short a,b,c,d;
 
  //INPUT AUCH  CHEMIKER ??
  a = intnums[x*4+0];  //a
  b = intnums[x*4+1];  //b
  c = intnums[x*4+2];  //c
  d = intnums[x*4+3];  //d
  if(a==b && b==c && c==d){
   swap_ints(x, sorted,  intval,  intnums);
   sorted++;
  }
 }
 *outf << sorted << " integrals after search for perm_1\n";
 sortcount[0] = sorted;

 //STEP3: BRING PERM1_5 thereafter
 for(long long int x = sorted; x < nrofint; x++){
  unsigned short a,b,c,d;

  //INPUT AUCH  CHEMIKER ??
  a = intnums[x*4+0];  //a
  b = intnums[x*4+1];  //b
  c = intnums[x*4+2];  //c
  d = intnums[x*4+3];  //d
  if(a==b && c==d){
   swap_ints(x, sorted,  intval,  intnums);
   sorted++;
  }
 }
 *outf << sorted << " integrals after search for perm_15\n";
 sortcount[1] = sorted;

 //STEP4: BRING PERM1234 thereafter
 for(int x = sorted; x < nrofint; x++){
  int a,b,c,d;
 
  //INPUT AUCH  CHEMIKER ??
  a = intnums[x*4+0];  //a
  b = intnums[x*4+1];  //b
  c = intnums[x*4+2];  //c
  d = intnums[x*4+3];  //d
  if(a==c && b==d){
   swap_ints(x, sorted,  intval,  intnums);
   sorted++;
  }
 }
 *outf << sorted << " integrals after search for perm_1234\n";
 sortcount[2] = sorted;

 //STEP5: BRING PERM1256 thereafter
 for(int x = sorted; x < nrofint; x++){
  int a,b; // c,d;
  
  //INPUT AUCH  CHEMIKER ??
  a = intnums[x*4+0];  //a
  b = intnums[x*4+1];  //b
//  c = intnums[x*4+2];  //c
//  d = intnums[x*4+3];  //d
  if(a==b){
   swap_ints(x, sorted,  intval,  intnums);
   sorted++;
  }
 }
 *outf << sorted << " integrals after search for perm_1256\n";
 sortcount[3] = sorted;
 
 outf->flush();
 *outf << "-------------------------------------------------------------------------------\n";
 *outf << "Undoing GAMESS - scaling\n";
 
 //PERM_1
 for(long long int x = 0; x < sortcount[0]; x++)
  intval[x] *= 8.;
 
 //PERM_15
 for(long long int x = sortcount[0]; x < sortcount[1]; x++)
  intval[x] *= 4.;
 
 //PERM_1234
 for(long long int x = sortcount[1]; x < sortcount[2]; x++)
  intval[x] *= 2.;
 
 //PERM_1256
 for(long long int x = sortcount[2]; x < sortcount[3]; x++)
  intval[x] *= 2.;
 
 *outf << "-------------------------------------------------------------------------------\n";
 outf->flush();
}
                                   
