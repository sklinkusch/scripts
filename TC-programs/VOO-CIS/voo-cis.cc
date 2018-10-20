
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
int           bpi;                          // Bytes per integral (from GAMESS logfile)
int           ipr;                          // Integrals per record (from GAMESS logfile)
int           nrorec;                       // Number of records (from GAMESS logfile)
int           cis_size;                     // size of CIS matrix
int           dst;                          // Number of desired states (states that are corrected)
int           type;                         // Type of correction to CIS energies
long long int nrofint;                      // Number of non-zero two-electron integrals
long long int sortcount[4];                 // boundaries for different permutation patterns
double        cisEgs;                       // CIS Ground State Energy

double*             intval;                 // two-electron integrals
unsigned short int* intnums;                // two-electron indices
double*             cistmpmat;              // temporary matrix for CI
double*             hmat;                   // one electron hamiltonian
double*             smat;                   // overlap matrix
char                f08file[256];           // GAMESS F08 file (contains two electron integrals)
char                f10file[256];           // GAMESS F10 file (contains one electron integrals)
char                wavfile[256];           // file containing MOs and MO energies
//char                binfile[256];           // file containing RHF energy and CI excitation energies

//Extern Functions

extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);                                                   // diagonalizes a matrix
//extern double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);                                               // calculates molecular orbital energy difference
//extern inline double RO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);                                    // calculates prefactors for CSF Slater determinants
//extern int pos_i_a(int i, int a, int omo, int umo, int llim);                                                               // calculates position of excitation i->a in CIS matrix
//extern double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens);    // calculates triples correction for i->a
/*extern double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim,
	                      double* cis_vec, double* prec_ints);                                                          // calculates doubles correction for ij->ab*/
/*extern double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en,
	                  double* prec_ints, double* MOens);                                                                // calculates first sum in equation 14, CPL 219, 21-29 (1994)*/
//extern double RO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens);                                    // calculates MP2 energy correction to the ground state
/*extern double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
	                      double* cis_vec, double cis_en, double* t1, double* t2);                                      // calculates second sum in equation 14, CPL 219, 21-29 (1994)*/
/*extern void SL_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
	                    double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2,                     // calculates CIS(D) correction in one single loop for all states
			                        int cis_size, ofstream* outf);*/
extern int calc_mo2el_ind_x(int i, int j, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index x
extern int calc_mo2el_ind_y(int k, int l, int omo, int umo, int llim);                                                      // routine to map ijkl for needed integrals to 2D file index y
extern double get_precalc_ints_sd(int i, int j, int k, int l,
	                           int omo, int umo, int llim, double* prec_ints);                                          // extract precalculated integrals safely
extern int* init_prex( int omo, int umo, int llim);
extern double get_precalc_ints_ovov(int i, int j, int k, int l,
	                             int omo, int umo, int llim, double* prec_ints);                                        // extract precalculated integrals fast ovov type
extern "C" void dafrd_(double* V, long long int* LEN, long long int* RECN, char* FNAME);                                    // reads data from F10 files
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);                                              // reads MO data from the wavfile
extern double* precalc_ints_sd(int omo, int umo, int llim, int nroao, double* MOs, long long int nrofint, 
	                       long long int* sortcount,  double* intval, unsigned short* intnums, ofstream* outf);         // routine to get the precalculated two-electron integrals
extern int rem_com(char* filename, char* streamstring, int string_length);                                                  // routine to read the input file
//extern void read_bin(char* binname, int nros, double& cisEgs, double* cisvals);                                             // routine to read by binfile
extern void status(ofstream* outf);                                                                                         // write status (date, time, and host) to the logfile
extern double calc_cisvec(int state, int initial, int final, int omo, int umo, int llim, int cis_size, 
	double* cisvecs);
extern double calc_Y(int I, int J, int a, int i, int omo, int umo, int llim, int cis_size, double* cis_vecs, 
	double* prec_ints);
extern double calc_theta(int I, int J, int a, int i, int omo, int umo, int llim, int cis_size, double* MOens, 
	double* cisens, double* cis_vecs, double* prec_ints);
extern void transformationCIS_SE(double* wav, double* SEwav, double* cisvecs, int cis_size);
extern void SEmax(double* SEwav, double* SEwavn, int cis_size, int nroao, int i, int a);
extern void renorm(double* SEwav, int cis_size);
extern int calc_diffs_dd(int i, int j, int k, int l, int a, int b, int c, int d, int nroao, int nroe);
extern void show_diff_dd(int i, int j, int k, int l, int a, int b, int c, int d, int p, int q, int nroao, int nroe);
extern void show_diffs_dd(int i, int j, int k, int l, int a, int b, int c, int d, int p, int q, int r, int s, int nroao, int nroe);
extern void produce_pops(int* pops, int i, int j, int a , int b, int nroao, int nroe);

//Functions

void read_integrals_12(void);                                                                                               // reads two-electron integrals from F08 file with 12 bytes per integral
void read_integrals_16(void);                                                                                               // reads two-electron integrals from F08 file with 16 bytes per integral
void resort_integrals(ofstream *outf);                                                                                      // routine to perform permutations of the two-electron integrals
void read_gamess(double* data, int recnr, int nromo, int type, char* filename);                                             // reads and resorts data from GAMESS F10 file

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

ist >> f08file >> f10file >> wavfile;                                                                                       // read the name of the F08 file and the wavefunction file
ist >> cisEgs >> nroao >> nroa >> bpi >> ipr >> nrorec >> nrofint;                                                          // read the name of the binfile, #AOs, #atoms, bytes per integral, integrals per record, #records, #integrals
ist >> nroe >> nros >> dst >> type >> llim >> ulim;                                                                         // read the #electrons, #electronic states, the lower and the upper limit of orbitals used for correlation 

outf << ".................................................................................\n";                              // prints the read data to log file (in order to show that everything is working right)
outf << "F08 file : " << f08file << "\n";
outf << "F10 file : " << f10file << "\n";
outf << "Wavefunction file: " << wavfile << "\n";
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

int aointl = nroao*nroao*6+nroa*5;
double* dumd = new double[aointl+1024];
int incre = 0;
hmat = &(dumd[incre]); incre += nroao*nroao;
smat = &(dumd[incre]); incre += nroao*nroao;

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
outf << "Reading one electron matrix elements\n";
read_gamess(hmat,11,nroao,1,f10file);
read_gamess(smat,12,nroao,1,f10file);
outf.flush();
outf << ".................................................................................\n";
resort_integrals(&outf);                                                                                                    // final sorting of integrals
// read data for perturbation part
outf << "VOO-CIS correction routine is starting\n";
if(nroe < 0)  nroe *= -1;

if(nroe%2 != 0){                                                                                                            // emergency exit if an unrestricted system is treated
 cerr << "Number of electrons not even!\n";
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
outf << "Number of electronic states to be corrected: " << dst << "\n";
switch(type){
    case 1: 
	outf << "CIS energies are corrected in a 1st order perturbation theory algorithm.\n";
	break;
    case 2:
	outf << "CIS energies are corrected through intrastate and interstate orbital relaxation.\n";
	break;
    default:
	outf << "CIS energies are not corrected.\n";
	type = 0;
}
//outf << "Starting to read binfile\n";
//outf.flush();
//double*             cisvals = new double[nros];                                                                             // CIS excitation energies
//read_bin(binfile, nros, cisEgs, cisvals);
//outf << "Read binfile\n";
//outf.flush();
outf << "CIS Ground State Energy: " << cisEgs << " E_h (" << cisEgs*evfactor << " eV)\n";                                   // write RHF ground state energy
outf.flush();
// prepatory part for perturbation
// read binfile and CI matrix file
double*             cismat  = new double[nros*nros];                                                                        // CI Matrix
double*             cisvecs = new double[nros*nros];                                                                        // CI Eigenvectors
double*             cisens  = new double[nros];
cis_size = nros;

for(int x = 0; x < nros; x++){
 for(int y = 0; y < nros; y++){
     cismat[x*cis_size+y] = 0.;                                                                                             // sets all CI matrix elements to zero
 }
}

  outf << "Calculating upper triangle of CIS-matrix\n";
long long int Cis_size = (long long int) cis_size;
omo = (nroe / 2) - llim;
umo = ulim + 1 - omo - llim;
double* prec_ints = precalc_ints_sd(omo, umo, llim, nroao, MOs, nrofint,  sortcount, intval, intnums, &outf);               // writes the precalculated integrals into the array prec_ints

                                                                                                                            // calculates the upper right triangle of the symmetric CIS matrix
    //// CODE FOR NORMAL LOOPS
  for(int x = 1 ; x < cis_size; x++){
   int i1 = (x-1)/umo+llim;
   int f1 = (x-1)%umo+omo+llim;
   for(int y = x ; y < cis_size; y++){
    int i2 = (y-1)/umo+llim;
    int f2 = (y-1)%umo+omo+llim;
    //## PRECALCULATED INTEGRAL --> biggest memory requierd
    cismat[x*Cis_size+y] +=   -get_precalc_ints_sd(i1, i2, f1, f2, omo, umo, llim, prec_ints);
    cismat[x*Cis_size+y] += 2.*get_precalc_ints_sd(i1, f1, i2, f2, omo, umo, llim, prec_ints);
    //## SEMI-DIRECT TRANSFORMATION
    //cismat[x*cis_size+y] +=   -mo2int_op(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //cismat[x*cis_size+y] += 2.*mo2int_op(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
    //## BRUT FORCE CALCULATIONS  --> lowest memory requierd
    //cismat[x*cis_size+y] +=   -calc_mo2int_bf(i1, i2, f1, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    //cismat[x*cis_size+y] += 2.*calc_mo2int_bf(i1, f1, i2, f2, nroao, MOs, nrofint, sortcount, intval, intnums);
    if(i1==i2 && f1==f2)  cismat[x*cis_size+y] += MOens[f1]-MOens[i1];
   }
   outf  << x << "\t";
   if(x%10==0) outf << "\n";
   outf.flush();
  }
 for(int x = 0; x < cis_size; x++){
  for(int y = x; y < cis_size; y++){
   cismat[y*cis_size+x] = cismat[x*cis_size+y];                                                                             // copys the upper right triangle to the lower left triangle (symmetric matrix)
  }
 }
outf << "\n";
outf.flush();
outf << "Starting to diagonalize CI matrix\n";
outf.flush();
diag_mat(cis_size, cismat, cisens, cisvecs);  
outf << "CI matrix diagonalized\n";
outf.flush();
// diagonalizes the CIS matrix to obtains the excitation energies (eigenvalues) ans the CIS vectors (eigenvectors)
// set up Hamiltonian for VOO correction
long long int dim;
long long int sqdim;
switch(type){
    case 0:
	dim = dst;
	outf << "Dimension of VOO Hamiltonian matrix is " << dim << "x1\n";
	outf.flush();
	break;
    case 1:
	dim = (1 + 2*dst)*(1 + 2*dst);
	sqdim = (1 + 2*dst);
	outf << "Dimension of VOO Hamiltonian matrix is " << sqdim << "x" << sqdim << "\n";
	outf.flush();
	break;
    case 2:
	dim = 1 + dst + dst*dst*dst + dst*dst;
	break;
}
outf << "Setting up VOO Hamiltonian matrix\n";
outf.flush();
double* voomat = new double[dim];       // VOO Hamiltonian matrix
//double* voosmat = new double[dim];      // VOO overlap matrix
#pragma omp parallel for
for(long long int i = 0; i < dim; i++){
    voomat[i] = 0.;                      // set all elements to zero
    //voosmat[i] = 0.;
}
if(type == 1){
  outf << "Setting up upper left block (ground state + single excitations)\n";
  outf.flush();
// upper left part (ground state and CIS states)
for(int i = 0; i < dst; i++){
#pragma omp parallel for
    for(int j = 0; j < dst; j++){
	voomat[i*sqdim+j] = cisens[i];      // all diagonal elements of the block contain the CIS excitation energies, others remain zero 
	//voosmat[i*sqdim+j] = 1.;            // all diagonal elements of the block are one, others remain zero
    }
}
outf << "Set up upper left block\n";
outf.flush();
// upper right part (couplings between single and double excitations)
outf << "Setting up upper right block (couplings between single and double excitations)\n";
outf.flush();
#pragma omp parallel for
for(int I = 0; I < dst; I++){
    for(int J = 0; J < dst; J++){
	for(int a = llim+omo; a < (llim+omo+umo); a++){
	    for(int i = llim; i < (llim+omo); i++){
		voomat[I*sqdim+dst+1+J] += calc_theta(J,J,i,a,omo,umo,llim,cis_size,MOens,cisens,cisvecs,prec_ints) * calc_Y(I,J,i,a,omo,umo,llim,cis_size,cisvecs,prec_ints);
	    }
        }
    }
}
outf << "Set up upper right block\n";
int diffs;
double integral;
double* wav_I = new double[cis_size];
double* wav_J = new double[cis_size];
double* SEwav_I = new double[cis_size];
double* SEwav_J = new double[cis_size];
int za, zb, zc, zd;
int Ii, Ji, If, Jf;
int* pops_I = new int[nroao];
int* pops_J = new int[nroao];
outf << "Setting up lower right block (double excitations)\n";
outf.flush();
#pragma omp parallel for
for(int x = dst; x < cis_size; x++){
    wav_I[x] = 0.;
    wav_J[x] = 0.;
}
outf << "Wavefunctions in higher states set to zero\n";
outf.flush();
for(int I = 1; I < dst; I++){
    for(int x = 0; x < dst; x++){
	if(x == I){
	    wav_I[x] = 1.;
	}else{
	    wav_I[x] = 0.;
	}
    }
    transformationCIS_SE(wav_I, SEwav_I, cisvecs, cis_size);
    for(int J = 1; J < dst; J++){
	for(int x = 0; x < dst; x++){
	    if(x == J){
		wav_J[x] = 1.;
	    }else{
		wav_J[x] = 0.;
	    }
	}
	transformationCIS_SE(wav_J, SEwav_J, cisvecs, cis_size); 
	for(int i = llim; i < (llim+omo); i++){
	    for(int a = llim+omo; a < (llim+omo+umo); a++){
		for(int j = llim; j < (llim+omo); j++){
		    for(int b = llim+omo; b < (llim+omo+umo); b++){
			integral = 0.;
			for(int x = 1; x < cis_size; x++){
			    for(int y = 1; y < cis_size; y++){
			      int Ii = (x-1)/umo+llim;
			      int If = (x-1)%umo+omo+llim;
			      int Ji = (y-1)/umo+llim;
			      int Jf = (y-1)%umo+omo+llim;
			      if(Ii == i) continue;
			      if(If == a) continue;
			      if(Ji == j) continue;
			      if(Jf == b) continue;
			      produce_pops(pops_I, Ii, i, If, a, nroao, nroe);
			      produce_pops(pops_J, Ji, j, Jf, b, nroao, nroe);
			      diffs = calc_diffs_dd(Ii, Ji, i, j, If, Jf, a, b, nroao, nroe);
			      if(diffs == 0){
				for(int A = 0; A < nroao; A++){
				  if(pops_I[A] == 1){
				  integral += SEwav_I[x] * SEwav_J[y] * hmat[A*nroao+A];
				  for(int B = 0; B < nroao; B++){
				    if(pops_J[B] == 1 && A != B && A < (llim+omo)){
     				    outf << "B = " << B << "\n";
				    outf.flush();
				    integral += SEwav_I[x] * SEwav_J[y] * 0.5 * get_precalc_ints_sd(A,B,A,B,omo,umo,llim,prec_ints);
				    }
				  }
				  }
				}
			      }else if(diffs == 1){
				show_diff_dd(Ii,Ji,i,j,If,Jf,a,b,za,zb, nroao, nroe);
				integral += SEwav_I[x] * SEwav_J[y] * hmat[za*nroao+zb];
				for(int B = 0; B < nroao; B++){
				  integral += SEwav_I[x] * SEwav_J[y] * get_precalc_ints_sd(za,B,zb,B,omo,umo,llim,prec_ints);
				}
			      }else if(diffs == 2){
				show_diffs_dd(Ii,Ji,i,j,If,Jf,a,b,za,zb,zc,zd,nroao,nroe);
				integral += SEwav_I[x] * SEwav_J[y] * get_precalc_ints_sd(za,zb,zc,zd,omo,umo,llim,prec_ints);
			      }else{
				integral += 0.;
			      }
			    }
			}
			voomat[(I+dst+1)*sqdim+dst+1+J] += calc_theta(I,I,j,b,omo,umo,llim,cis_size,MOens,cisens,cisvecs,prec_ints)
			* calc_theta(J,J,i,a,omo,umo,llim,cis_size,MOens,cisens,cisvecs,prec_ints) * integral;
		    }
		}
	    }
	}
    }
}
outf << "Set up lower right block\n";
// construct lower left part of the hamiltonian matrix from the upper right part
outf << "Setting up lower left block (transpose of upper right block)\n";
outf.flush();
#pragma omp parallel for
for(int I = 0; I < (1+dst); I++){
  for(int J = (1+dst); J < sqdim; J++){
    voomat[J*sqdim+I] = voomat[I*sqdim+J];
  }
}
outf << "Set up lower left block\n";
outf.flush();
double* voovecs = new double[dim];
double* voovals = new double[dst];
outf << "Diagonalizing VOO Hamiltonian matrix\n";
outf.flush();
diag_mat(sqdim,voomat,voovals, voovecs);
outf << "VOO Hamiltonian matrix diagonalized\n";
outf << "Printing corrected energes...\n";
outf.flush();
for(int i = 0; i < dst; i++){
  outf << "State " << i << " has a corrected energy of " << voovals[i] << " Hartree (" << voovals[i]*evfactor << " eV) from " << cisens[i] << " Hartree (" << cisens[i]*evfactor << " eV)\n";
  outf.flush();
}
status(&outf);
outf.flush();
outf.close();
}
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
                                  
void read_gamess(double* data, int recnr, int nromo, int type, char* filename){
 long long int len;
 //SET CORRECT LENGTH
 if(type == 0) len = nromo;
 if(type == 1) len = (nromo*nromo+nromo)/2;
 if(type == 2) len = nromo*nromo;

 //len = 2;

 char dumn[513];
 int fnamelen = strlen(filename);
 sprintf(dumn,"%s",filename);
 for(int x = fnamelen; x < 512; x++) dumn[x] = ' ';
 dumn[512] = 0;

 long long int Recnr = recnr;

 dafrd_(data, &len, &Recnr, dumn);


 //RESORT MATRIX
 if(type == 1){
  int dcount = len-1;
  for(int x = nromo-1; x >= 0; x--){
   for(int y = x; y >= 0; y--){
    data[x*nromo+y] = data[dcount--];
   }
  }
  for(int x = 0; x < nromo; x++){
   for(int y = x; y < nromo; y++) data[x*nromo+y] = data[y*nromo+x];
  }
 }
}

