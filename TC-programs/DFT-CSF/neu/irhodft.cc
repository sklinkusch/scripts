/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: istochdft.cc                                                           *
 *                                                                              *
 * RHF program                                                                  *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 *                                                    Stefan Klinkusch 2015     *
 ********************************************************************************/ 
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

//Functions
void cp_mus(int cis_size, int nstates, int nroao, int omo, int umo, int llim,  int nroe, double* cismat_x, 
	    double* cismat_y, double* cismat_z, double* cisvecs, double* Dx, double* Dy, double* Dz, 
	    double* mu_core, double* MOs, double* Pmat, double* tmpvecs, double* dipxcsf, double* dipycsf,
	    double* dipzcsf, ofstream* outf, long int &nonzero, int mode);


//Extern Functions
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void status(ofstream* outf);
extern void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
		     double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		     double *Dz); //, long long int* sortcount, double* intval, 
//		     unsigned short* intnums);
extern void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
			 double* mu_core);
extern double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old, 
			      double* MOs, double damp);
extern double calc_op_1el(int nroao, double* opmat, double* Pmat);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
                     double* cismat_x, double* cismat_y, double* cismat_z, double* mumat_x, 
		     double* mumat_y, double* mumat_z, double* mucore, 
                     double* MOs, double* Pmat, double* tmpvec_ao, ofstream* outf);
extern void calc_mu_mat_cis_crs(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
			    double *cismat, double* mumat, double mucore, 
			    double* MOs, double* Pmat, double* tmpvec_ao, double* dip_val, ofstream* outf);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mt(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf);
extern void trans_mt_mult(int nr_mat, int nc_mat, int nr_trans, int nc_trans, double* matx, double* maty, double* matz, double* trans, double* tmpmatx, double* tmpmaty,
	double* tmpmatz, int dir, ofstream* outf);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mt(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* mat, double* trans, double* tmpmat, int dir, ofstream* outf);
extern void trans_mt_mult(long long int nr_mat, long long int nc_mat, long long int nr_trans, long long int nc_trans, double* matx, double* maty, double* matz, double* trans, 
	double* tmpmatx, double* tmpmaty, double* tmpmatz, int dir, ofstream* outf);
extern long int  find_nonzeros(int dim, int umo, int omo, int llim);
extern void rows_cols_crs(int cis_size, int umo, int omo, int llim, long int* col_ind, long int* row_ptr);
extern void sparse_matrix(int nx, int ny, double* matrix, long int nz, double* nonzero, long int* col_ind, long int* row_ptr);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  int   nbasis;           //Nr of configurations
  int  nstates;           //Nr of electronic states
  int     nroe;           //Nr of electrons  (if negative, read in center of mass)
  int     llim;           //first MO used for correlation
  int     ulim;           //last  MO used for correlation
  char sysfile[256];      //binary system file
  char wavfile[256];      //HF-Wavefunction file
  char vecfile[256];      //Eigenvalue and Eigenvector file
  int calc_COM = 0;       // 0 -> calc COM otherwise read in
  double damp = 1.;

  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  sprintf(dumc,"%s.sys",argv[2]);
  ofstream syf(dumc);

  outf << "DFT-CSF [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input from " << argv[1] << "\n\n";
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  
  ist >> sysfile >> nbasis >> nstates >> nroe;

  double center_of_mass[3]; 

  if(nroe < 0){
    nroe *= -1;
    calc_COM = 1;
    for(int x = 0; x <3 ; x++) ist  >> center_of_mass[x];
    outf << "Center of mass read in!!!\n";
  }
  
  if(nroe%2 != 0){
    cerr << "Nr of electrons not even !!!\n";
    exit(2);
  }
  
  ist >> llim >> ulim >> wavfile >> vecfile;
  outf << "System data \nNr of electrons: " << nroe << "\n";
  outf << "Number of configurations: " << nbasis << "\n";
  outf << "Number of electronic states: " << nstates << "\n";
  outf << "CIS data\nLower limit of MOs used for correlation: " << llim << "\n";
  outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
  outf << "Reading HF-wave function from  " << wavfile << "\n";
  outf << "Reading eigenvalues and vectors from " << vecfile << "\n";
  ifstream vecf(vecfile);
  outf << "-------------------------------------------------------------------------------\n";

  //SCF VARIABLES 

  //system size
  int    nroao;                 //Nr of basis functions
  int    nroa;                  //Nr of atoms 
  long long  int nrofint;      //Nr of two electron Integrals
  
  get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

  outf << "System sizes read from " << sysfile << "\n";
  outf << "Nr of basis functions: " << nroao   << "\n";
  outf << "Nr of atoms: " << nroa << "\n";
  outf << "Nr of non zero 2el integrals " << nrofint << "\n";
  outf << "Allocating Memory\n";
  int omo = nroe/2 - llim;
  int umo = ulim + 1 - llim - omo;
  int cis_size = omo * umo + 1;
  if(cis_size != nbasis){
      cerr << "Wrong limits for the orbitals\n";
      exit(24);
  }
  long long int Cis_size = (long long int) cis_size;

  //atoms
  double*        coord;         //atomic coordinats               3*nroa
  double*        charges;       //atomic charges                    nroa
  double*        mass;          //atomic masses                     nroa

  //one electron mat&vecs
  double*        Hmat;          //Hamiltonian matrix                nroao*nroao
  double*        Tmat;          //Kinetic energy matrix             nroao*nroao
  double*        Smat;          //Overlap matrix                    nroao*nroao

  double*        Pmat;          //density matrix                    nroao*nroao
  double*        Pmat_old;      //Old density matrix                nroao*nroao

  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
//  double*        MOens_old;     //Old MO energies;                  nroao


  //Temporary memory spaces
  double*  tmpmat1;                //                               nroao*nroao
 // double*  tmpmat2;                //                               nroao*nroao
//  double*  tmpvecs;                //                               nroao*nroao

//  double*  tmpvals;                //                               nroao
 
  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+14*nroao*nroao+3*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]); inc += 3*nroa; charges = &(dumd[inc]); inc += nroa; mass  = &(dumd[inc]); inc += nroa;
  Hmat = &(dumd[inc]); inc += nroao*nroao; Tmat = &(dumd[inc]); inc += nroao*nroao; Smat = &(dumd[inc]); inc += nroao*nroao;
  Pmat = &(dumd[inc]); inc += nroao*nroao;Pmat_old = &(dumd[inc]); inc += nroao*nroao;
  MOs = &(dumd[inc]); inc += nroao*nroao; Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao; 
  
  MOens = &(dumd[inc]); inc += nroao; 
 // MOens_old = &(dumd[inc]); inc += nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; 
//  tmpmat2 = &(dumd[inc]); inc += nroao*nroao; 
//  tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
//  tmpvals = &(dumd[inc]); inc += nroao;

//  double*         intval        = new double[nrofint];                     //two electron integrals
//  unsigned short* intnums       = new unsigned short[nrofint*4];           //two electron indices
//  long long int   sortcount[4];                                            //num of two electron Integrals in each perm. type

  
  //Read system data 
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz); //, sortcount, intval, intnums);

  //System Properties for HF
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Ground state Properties:\n";
  outf << "...............................................................................\n";
  outf << "Ion Cores: coordinates mass charge\n";
  int count = 0;
  for(int x = 0; x < nroa; x++){
    outf << x << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%.5f", mass[x]);           outf  << dumc << "\t";
    sprintf(dumc,"%.0f", charges[x]);        outf  << dumc << "\n";     
  }
  
  outf.precision(10);

  double mu_core[3];
  if(calc_COM == 0)
    calc_center_of_mass( nroa, coord,  mass,  center_of_mass);
  calc_mu_core( nroa, coord, charges, center_of_mass, mu_core);
  outf << "Center of mass (x,y,z): " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << "\n";
  outf << "Core dipole moment (at C.o.M.): " << mu_core[0] << " " << mu_core[1] << " " << mu_core[2] << "\n";
  outf << "...............................................................................\n";
  outf << "Init HF operators:\n";
  read_wav_HF(wavfile, nroao, MOens, MOs);
  outf << "HF-wave function  read from " << wavfile << "\n";
  outf << "Building density and Fock matrix\n";
  build_Pmat_dscf( nroao, nroe,  Pmat,  Pmat_old,  MOs, 0.);
  outf.flush();
  outf << "MO-energies:\n";
  for(int x = 0; x < nroao; x++){
    sprintf(dumc,"%+.5f ",MOens[x]);
    outf << dumc << "\t";
    if((x+1)%5==0) outf << "\n";
  }
  
  outf << "\n\nDipole moment (AU): " 
       <<  -calc_op_1el( nroao, Dx, Pmat) +mu_core[0] << " " 
       <<  -calc_op_1el( nroao, Dy, Pmat) +mu_core[1] << " " 
       <<  -calc_op_1el( nroao, Dz, Pmat) +mu_core[2] << "\n";
  outf << "Dipole moment (Debye): " 
       <<  -(calc_op_1el( nroao, Dx, Pmat) -mu_core[0])*2.5417462 << " " 
       <<  -(calc_op_1el( nroao, Dy, Pmat) -mu_core[1])*2.5417462 << " " 
       <<  -(calc_op_1el( nroao, Dz, Pmat) -mu_core[2])*2.5417462 << "\n";
  
  outf << "-------------------------------------------------------------------------------\n";
  sprintf(dumc,"%s.bcs",argv[2]);
  outf << "Writing binary data to " << dumc << "\n";
  outf.flush();
  ofstream datf(dumc);
  datf.write((char *) &nroao, sizeof(int));
  datf.write((char *) &nroe, sizeof(int));
  datf.write((char *) &llim, sizeof(int));
  datf.write((char *) &ulim, sizeof(int));
  datf.write((char *) &nstates, sizeof(int));
  outf << "System sizes written to bcs-file\n";
  outf.flush();

  //Write data to data-file

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Allocating Memory for CIS:\n";
  outf.flush();
 
//  double*  cismat    = new double[Cis_size*Cis_size];     //Matrix for Operator representations in CIS-Basis
  double*  corr      = new double[Cis_size*Cis_size];     //CIS(D) correction dummy (not needed here, filled with zeroes)
//  double*  cistmpmat = new double[Cis_size*Cis_size];     //Temp space
  double*  cisvals   = new double[cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++) corr[i*cis_size+j] = 0.;
  }
  outf << "Corrections set to zero\n";
  outf.flush();

  int real_nroao, nros, nroctot;
  vecf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao > nroao){
      cerr << "Wrong nr of AOs in vec file (" << real_nroao << " vs. " << nroao << ")!\n";
      exit(111);
  }
  vecf.read((char *) &nros, sizeof(int));
  if(nros != nstates){
      cerr << "Wrong nr of states in vec file!\n";
      exit(112);
  }
  vecf.read((char *) &nroctot, sizeof(int));
  if(nroctot != nbasis){
      cerr << "Wrong nr of configurations in vec file!\n";
      exit(113);
  }
  vecf.read((char *) cisvals, sizeof(double)*nstates);
  vecf.read((char *) cisvecs, sizeof(double)*nstates*Cis_size);
  outf << "Eigenvalues and vectors read\n";
  outf.flush();
  datf.write((char *) cisvals, nstates*(long long int) sizeof(double));
  datf.write((char *) cisvecs, nstates*Cis_size*(long long int) sizeof(double));
  datf.flush();
  outf << "Eigenvalues and vectors written\n";
  outf.flush();

  double* cistmpvals = new double[Cis_size*Cis_size];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  if(x == y && x < nstates){
	      cistmpvals[x*cis_size+y] = cisvals[x];
	  }else{
	      cistmpvals[x*cis_size+y] = 0.;
	  }
      }
  }
// trans_mt(nstates, nbasis, nbasis, cistmpvals, cisvecs, cismat, 1);
// outf << "CI matrix built\n";


  outf << "-------------------------------------------------------------------------------\n";
  outf.flush();

   //##CALCULATION OF IONIZATION RATES / WRITE OUT TO IRX-FILE / WITHOUT DOUBLES CORRECTIONS
  double* gamma_csf = new double[cis_size];
  double* gamma = new double[nstates];
  double dumm = 0.;
  int MO_r;
  double ip = MOens[nroe/2-1];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      gamma_csf[x] = 0.;
      MO_r = (x-1)%umo+omo+llim;
      if(MOens[MO_r] > 0.) gamma_csf[x] = damp*sqrt(MOens[MO_r]);
  }
outf << "Ionization rates for CSFs calculated\n";
outf.flush();
  for(int x = 0; x < nstates; x++){
      if(cisvals[x] > ip){
	  dumm = 0.;
	  for(int y = 0; y < nbasis; y++) dumm += pow(cisvecs[x*nbasis+y],2.)*gamma_csf[y];
	  gamma[x] = dumm;
      }else{
	  gamma[x] = 0.;
      }
  }
  sprintf(dumc, "%s.irx", argv[2]);
  ofstream irxf;
  irxf.open(dumc);
  for(int x = 0; x < nstates; x++) irxf << x << " " << gamma[x] << "\n";
  irxf.close();
  outf << "Ionization rates for eigenstates calculated and written\n";
  outf.flush();

  double* energies = new double[cis_size];
  for(int x = 0; x < cis_size; x++) energies[x] = cisvals[x];

    outf << "-------------------------------------------------------------------------------\n";
  long int nonzero;
  outf << "Building dipole transition matrices!\n";
  outf << "X\n";
  outf.flush();
  double* cismat_x = new double[cis_size*cis_size];
  double* cismat_y = new double[cis_size*cis_size];
  double* cismat_z = new double[cis_size*cis_size];
  double* dipxcsf  = new double[cis_size*cis_size];
  double* dipycsf  = new double[cis_size*cis_size];
  double* dipzcsf  = new double[cis_size*cis_size];
  cp_mus(cis_size, nstates, nroao,  omo,  umo, llim, nroe,  cismat_x, cismat_y, cismat_z,  cisvecs, Dx, Dy, Dz,
	 mu_core, MOs, Pmat, tmpmat1, dipxcsf, dipycsf, dipzcsf, &outf, nonzero, 2);
  double* dipx = new double[nstates*nstates];
  double* dipy = new double[nstates*nstates];
  double* dipz = new double[nstates*nstates];
  double* dipxvals = new double[nstates];
  double* dipyvals = new double[nstates];
  double* dipzvals = new double[nstates];
  double* dipxvecs = new double[nstates*nstates];
  double* dipyvecs = new double[nstates*nstates];
  double* dipzvecs = new double[nstates*nstates];
  for(int x = 0; x < nstates; x++){
      for(int y = 0; y < nstates; y++){
	 dipx[x*nstates+y] = cismat_x[x*nstates+y];
	 dipy[x*nstates+y] = cismat_y[x*nstates+y];
	 dipz[x*nstates+y] = cismat_z[x*nstates+y];
      }
  }
  diag_mat(nstates, dipx, dipxvals, dipxvecs);
  diag_mat(nstates, dipy, dipyvals, dipyvecs);
  diag_mat(nstates, dipz, dipzvals, dipzvecs);
  datf.write((char *) dipxvals, nstates*sizeof(double));
  datf.write((char *) dipxvecs, nstates*nstates*(long long int) sizeof(double));
  datf.write((char *) dipyvals, nstates*sizeof(double));
  datf.write((char *) dipyvecs, nstates*nstates*(long long int) sizeof(double));
  datf.write((char *) dipzvals, nstates*sizeof(double));
  datf.write((char *) dipzvecs, nstates*nstates*(long long int) sizeof(double));
  datf.flush();

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Execution started on/at:\n";
  status(&outf);
  outf.flush();
}

void cp_mus(int cis_size, int nstates, int nroao, int omo, int umo, int llim, int nroe, double* cismat_x, double* cismat_y, double* cismat_z, double* cisvecs, double* Dx,
	    double* Dy, double* Dz, double* mu_core, double* MOs, double* Pmat, double* tmpvecs, double* dipxcsf, double* dipycsf, double* dipzcsf, ofstream* outf,
	    long int &nonzero, int mode){
  long long int Cis_size = cis_size;
  long long int Nstates = (long long int) nstates;
/*  if(mode == 0){
      int cpo = cis_size+1;
      long int* dip_row = new long int[cpo];
      nonzero = find_nonzeros(cis_size, umo, omo, llim);
      cout << "Nonzeroes: " << nonzero << "\n";
      long int* dip_col = new long int[nonzero];
      rows_cols_crs(cis_size, umo, omo, llim, dip_col, dip_row);
  }
  if(mode == 1){
  double* dip_val = new double[nonzero];
  calc_mu_mat_cis_crs(cis_size, nroao, omo, umo, llim, nroe, cismat, D, mu_core, MOs, Pmat, tmpvecs, dip_val, outf);
  if(cis_size < 10000)
    trans_mt(cis_size, cis_size, nstates, cis_size, cismat, cisvecs, cistmpmat, 1, outf); //Forward transfornation
  else
    trans_mt(Cis_size, Cis_size, Nstates, Cis_size, cismat, cisvecs, cistmpmat, 1, outf); //Forward transfornation
  }*/
  if(mode == 2){
      calc_mu_mat_cis(cis_size, nroao, omo, umo, llim, nroe, cismat_x, cismat_y, cismat_z, Dx, Dy, Dz, mu_core, MOs, Pmat, tmpvecs, outf);
      *outf << "Dipole matrix in CSF space calculated\n";
      outf->flush();
      if(cis_size < 10000){
	  trans_mt_mult(cis_size, cis_size, nstates, cis_size, cismat_x, cismat_y, cismat_z, cisvecs, dipxcsf, dipycsf, dipzcsf, 1, outf);
//	  trans_mt(cis_size, cis_size, nstates, cis_size, cismat, cisvecs, cistmpmat, 1, outf);
      }else{
	  trans_mt_mult(Cis_size, Cis_size, Nstates, Cis_size, cismat_x, cismat_y, cismat_z, cisvecs, dipxcsf, dipycsf, dipzcsf, 1, outf);
//	  trans_mt(Cis_size, Cis_size, Nstates, Cis_size, cismat, cisvecs, cistmpmat, 1, outf);
      }
      *outf << "Dipole matrix transformed to eigenspace\n";
      outf->flush();
  }
}

