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
int           nroe;                         // Number of electrons
int           nros;                         // Number of electronic states
int           nrods;                        // Number of electronic states desired for CIS(D) correction
int           llim;                         // lower limit of MOs used for correlation
int           ulim;                         // upper limit of MOs used for correlation
int           umo;                          // Number of unoccupied (virtual) MOs
int           omo;                          // Number of occupied MOs
int           cis_size;                     // size of CIS matrix
uint32_t      nonzero;                      // number of nonzero two-electron integrals
uint32_t      x_max;
uint32_t      y_max;
double        cisEgs;                       // RHF ground state energy

double*             cistmpmat;              // temporary matrix for CI
char                intfile[256];           // GAMESS F08 file
char                wavfile[256];           // file containing MOs and MO energies

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
extern double RO_calc_MP2(int omo, int umo, int llim, double* ints_ovov, double* MOdelta);                                         // calculates MP2 energy correction to the ground state
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
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                                              // reads MO data from the wavfile
extern int rem_com(char* filename, char* streamstring, int string_length);                                                  // routine to read the input file
//extern void read_bin(char* binname, int nros, double& cisEgs, double* cisvals);                                             // routine to read by binfile
extern void status(ofstream* outf);                                                                                         // write status (date, time, and host) to the logfile
extern void SK_calc_cis_dc(int x, int omo, int umo, int llim, double* ints_ooov, double* ints_ovvv, double* MOdelta, double* cisvec, double* cisen, int cis_size, double* t1);
extern void SK_calc_cis_tc(int x, int omo, int umo, int llim, double* cisvec, int cis_size, double *t2, double* twoel, double* MOdelta);
extern void read_moint_sizes(uint32_t& nonzero, uint32_t& xmax, uint32_t& ymax, ifstream* inf);
extern void read_moints(uint32_t nonzero, uint32_t xmax, uint32_t ymax, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, ifstream* inf);
extern void calc_MOdelta(int omo, int umo, int llim, double* MOdelta, double* MOens);
extern long long int map_oovv(int i, int j, int a, int b, int omo, int umo, int llim);
extern long long int map_ovov(int i, int a, int j, int b, int omo, int umo, int llim);
extern long long int map_ooov(int i, int j, int k, int a, int omo, int umo, int llim);
extern void map_ints_ooov(int omo, int umo, int llim, double* ints_ooov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
extern long long int map_ovvv(int i, int a, int b, int c, int omo, int umo, int llim);
extern void map_ints_ovvv(int omo, int umo, int llim, double* ints_ovvv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
extern void map_ints_ovov(int omo, int umo, int llim, double* ints_ovov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
extern void map_ints_oovv(int omo, int umo, int llim, double* ints_oovv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);

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

// Read input file
int buff_length = 65536;
char* file_buff = new char[buff_length];                                                                                    // file buffer to read input file and to exclude comments marked using '#'
rem_com(argv[1], file_buff, buff_length);                                                                                   // read input file to write everything from the input file (excluding comments) to the file buffer
istringstream ist(file_buff);                                                                                               // stringstream to read from the file buffer

ist >> intfile >> wavfile;                                                                                                  // read the name of the F08 file and the wavefunction file
ist >> cisEgs >> nroao;                                                          // read the RHF ground state energy, #AOs, #atoms, bytes per integral, integrals per record, #records, #integrals
ist >> nroe >> nros >> nrods >> llim >> ulim;                                                                               // read the #electrons, #electronic states, #desired electronic states the lower and the upper limit of orbitals used for correlation 

outf << ".................................................................................\n";                              // prints the read data to log file (in order to show that everything is working right)
outf << "Integrals file : " << intfile << "\n";
outf << "Wavefunction file: " << wavfile << "\n";
outf << "Hartree-Fock ground state energy: " << cisEgs << " E_h (" << cisEgs*evfactor << " eV)\n";
outf << "Nr of basis functions: " << nroao << "\n";
outf << "Number of electrons: " << nroe << "\n";
outf << "Lower limit of MOs used for correlation: " << llim << "\n";
outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
outf << "Number of electronic states: " << nros << "\n";
outf << "Number of desired electronic states: " << nrods << "\n";
outf.flush();

// make number of electrons positive
if(nroe < 0)  nroe *= -1;
// make sure that the system is restricted
if(nroe%2 != 0){                                                                                                            // emergency exit if an unrestricted system is treated
 cerr << "Number of electrons not even!\n";
 exit(7);
}

// read WF file
double*             MOs = new double[nroao*nroao];                                                                          // molecular orbitals
double*             MOens = new double[nroao];                                                                              // MO energies
read_wav_HF(wavfile, nroao, MOens, MOs);                                                                                    // read MOs and MO energies from wavefunction file
double*             cismat  = new double[nros*nros];                                                                        // CI Matrix
cis_size = nros;
delete [] MOs;

#pragma omp parallel for
for(int x = 0; x < nros; x++){
 for(int y = 0; y < nros; y++){
     cismat[x*cis_size+y] = 0.;                                                                                             // sets all CI matrix elements to zero
 }
}
long long int Cis_size = (long long int) cis_size;
omo = (nroe / 2) - llim;
umo = ulim + 1 - omo - llim;

// read data from MO-integral file
ifstream intf(intfile);
read_moint_sizes(nonzero, x_max, y_max, &intf);
double* prec_vals = new double[nonzero];
uint32_t* prec_cols = new uint32_t[nonzero];
uint32_t* prec_rows = new uint32_t[x_max+1];
read_moints(nonzero, x_max, y_max, prec_vals, prec_cols, prec_rows, &intf);

// Write (oo|vv) and (ov|ov) integrals into separate arrays with their own mappings
outf << "Storing (oo|vv) and (ov|ov) MO integrals in special arrays\n";
outf.flush();
long long int dim_oovv = map_oovv(llim+omo-1, llim+omo-1, llim+omo+umo-1, llim+omo+umo, omo, umo, llim);
long long int dim_ovov = map_ovov(llim+omo-1, llim+omo+umo-1, llim+omo-1, llim+omo+umo, omo, umo, llim);
double* ints_oovv = new double[dim_oovv];
double* ints_ovov = new double[dim_ovov];
map_ints_oovv(omo, umo, llim, ints_oovv, prec_vals, prec_cols, prec_rows);
map_ints_ovov(omo, umo, llim, ints_ovov, prec_vals, prec_cols, prec_rows);
long long int mapvar;
outf << "Building CI matrix\n";
outf.flush();
    //// CODE FOR NORMAL LOOPS
#pragma omp parallel for
  for(int x = 1 ; x < cis_size; x++){
   int i1 = (x-1)/umo+llim;                         // min(i1) = llim, max(i1) = HOMO, loop over occupied orbitals
   int f1 = (x-1)%umo+omo+llim;                     // min(f1) = llim + omo = LUMO, max(f1) = HUMO, loop over virtual orbitals
   for(int y = x ; y < cis_size; y++){
    int i2 = (y-1)/umo+llim;                        // min(i2) = i1 =(min)=> llim, max(i2) = HOMO, loop over (higher) occupied orbitals
    int f2 = (y-1)%umo+omo+llim;                    // min(f2) = i2 =(min)=> LUMO, max(f2) = HUMO, loop over (higher) virtual orbitals
    //## PRECALCULATED INTEGRAL --> biggest memory requierd
    mapvar = map_oovv(i1, i2, f1, f2, omo, umo, llim);
    cismat[x*Cis_size+y] += -1.*ints_oovv[mapvar]; //  -get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (oovv)
    mapvar = map_ovov(i1, f1, i2, f2, omo, umo, llim);
    cismat[x*Cis_size+y] += 2.*ints_ovov[mapvar]; // get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_vals, prec_cols, prec_rows); // get prec_ints for (ovov)
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
double* MOdelta = new double[dim_oovv];
calc_MOdelta(omo, umo, llim, MOdelta, MOens);
delete [] MOens;
#pragma omp parallel for
 for(int x = 0; x < cis_size; x++){
  for(int y = x; y < cis_size; y++){
   cismat[y*cis_size+x] = cismat[x*cis_size+y];                                                                             // copys the upper right triangle to the lower left triangle (symmetric matrix)
  }
 }
outf << "Done\n";
outf << "Diagonalizing CI matrix\n";
outf.flush();
delete [] ints_oovv;
double*             cisvecs = new double[nros*nros];                                                                        // CI Eigenvectors
double*             cisens  = new double[nros];
diag_mat(cis_size, cismat, cisens, cisvecs);
delete [] cismat;
outf << "Calculating ground state correction\n";
double* cisd_corr = new double[nros];                                                                                       // CIS(D) correction
for(int x = 0; x < cis_size; x++) cisd_corr[x] = 0.;                                                                        // sets initially all correction values to zero
double E2 = RO_calc_MP2(omo, umo, llim, ints_ovov, MOdelta);                                                                    // calculates ground state MP2 energy
double* C_t1 = new double[nros];                                                                                            // first term of the CIS(D) correction (double corrections)
double* C_t2 = new double[nros];                                                                                            // second term of the CIS(D) correction (triple corrections)
for(int i = 0; i < nros; i++){
    C_t1[i] = 0.;
    C_t2[i] = 0.;
}
double t2;
int count_det = 0;
outf << "............................................................................\n";
status(&outf);
outf << "Calculating triple corrections\n";
outf.flush();
for(int x = 1; x <= nrods; x++){
    SK_calc_cis_tc(x, omo, umo, llim, cisvecs, cis_size, &t2, ints_ovov, MOdelta);
    C_t2[x] = t2;
    outf << x << "\t";
    count_det++;
    if(count_det%10 == 0) outf << count_det << "\n";
    outf.flush();
}
outf << "\n";
status(&outf);
outf.flush();
delete [] ints_ovov;
outf << "Starting mapping for (oo|ov) and (ov|vv) integrals: \n";
outf.flush();
double* ints_ooov = new double[omo*omo*omo*umo];
double* ints_ovvv = new double[omo*umo*umo*umo];
map_ints_ooov(omo, umo, llim, ints_ooov, prec_vals, prec_cols, prec_rows);
map_ints_ovvv(omo, umo, llim, ints_ovvv, prec_vals, prec_cols, prec_rows);
delete [] prec_vals;
delete [] prec_cols;
delete [] prec_rows;
double t1;
outf << "Calculating doubles corrections\n";
outf.flush();
int count_doub = 0;
for(int x = 1; x <= nrods; x++){
    SK_calc_cis_dc(x, omo, umo, llim, ints_ooov, ints_ovvv, MOdelta, cisvecs, cisens, cis_size, &t1);
    C_t1[x] = t1;
    outf << x << "\t";
    count_doub++;
    if(count_doub%10 == 0) outf << count_doub << "\n";
    outf.flush();
}
delete [] ints_ovvv;
delete [] ints_ooov;
delete [] MOdelta;
// Perturbation theory part -> only implemented for precalculated integrals so far
//SINGLE LOOP CODE
//long long int dimens = calc_mapvar(llim+omo-1, llim+omo-1, llim+omo+umo-1, llim+omo+umo, omo, umo, llim);
//SL_calc_cis_dc(omo, umo, llim, prec_vals, prec_cols, prec_rows, MOens, cisvecs, cisens, cistmpmat, C_t1, C_t2, cis_size, &outf);                  // calculates terms for CIS(D) correction for each state
outf << "\n";
outf << "Ground state MP2 energy (RO-OPS): T2 " << E2 << " Hartree (" << E2*evfactor << " eV)  TOT " << E2+cisEgs << " Hartree (" << (E2+cisEgs)*evfactor << " eV)\n";
outf << "...............................................................................\n";
status(&outf);
outf.flush();
outf << "Calculating CIS(D) corrections\n";
outf.flush();
#pragma omp parallel for
for(int x = 0; x < cis_size; x++) cisd_corr[x] = C_t1[x] + C_t2[x];                                                         // calculates total CIS(D) correction for each state
for(int x = 1; x < cis_size; x++){                                                                                          // prints out CIS(D) energies for all states
outf << "Calculated CIS(D) corrections\n";
outf.flush();
 outf << "State " << x << " : \n";
 outf << "doubles correction is " << C_t1[x] << " E_h (" << C_t1[x]*evfactor << " eV)\n";
 outf << "triples correction is " << C_t2[x] << " E_h (" << C_t2[x]*evfactor << " eV)\n";
 outf << "Total corr " << cisd_corr[x] << " E_h (" << cisd_corr[x]*evfactor << " eV), excit.en. " << cisens[x]+cisd_corr[x] << " E_h (" << (cisens[x]+cisd_corr[x])*evfactor << "eV) from " << cisens[x] << " E_h (" << cisens[x]*evfactor << "eV)\n";
 outf << ".............................................................................\n";
 outf.flush();
}
status(&outf);
outf.flush();
outf.close();
}
                                   
