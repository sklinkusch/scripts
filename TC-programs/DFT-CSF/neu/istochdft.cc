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
void cp_mus(int cis_size, int nroao, int omo, int umo, int llim,  int nroe, double* cismat, double* cisvecs, double* D,
	    double mu_core, double* MOs, double* Pmat, double* tmpvecs, double* cistmpmat,
	    ofstream *syf, long int &nonzero, int mode);


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
                     double *cismat, double* mumat, double mucore, 
                     double* MOs, double* Pmat, double* tmpvec_ao);
extern void calc_mu_mat_cis_crs(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
			    double *cismat, double* mumat, double mucore, 
			    double* MOs, double* Pmat, double* tmpvec_ao, double* dip_val);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
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
  
  syf.write((char *) &nbasis, sizeof(int));   // writing nbasis to sys file for stochastic TDCI
  syf.write((char *) &nstates, sizeof(int));  // writing nstates to sys file for stochastic TDCI
  double* configurations = new double[nbasis];
  for(int x = 0; x < nbasis; x++) configurations[x] = 0.;
#pragma omp parallel for
  for(int x = 1; x < nbasis; x++){
      int i = (x-1)/umo+llim;
      int f = (x-1)%umo+omo+llim;
      configurations[x] = MOens[f] - MOens[i];
  }
  syf.write((char *) configurations, sizeof(double)*cis_size);

  //Write data to data-file

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Allocating Memory for CIS:\n";
 
  double*  cismat    = new double[Cis_size*Cis_size];     //Matrix for Operator representations in CIS-Basis
  double*  corr      = new double[Cis_size*Cis_size];     //CIS(D) correction dummy (not needed here, filled with zeroes)
  double*  cistmpmat = new double[Cis_size*Cis_size];     //Temp space
  double*  cisvals   = new double[cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors
  for(int i = 0; i < cis_size; i++){
      for(int j = 0; j < cis_size; j++) corr[i*cis_size+j] = 0.;
  }

  int real_nroao, nros, nroctot;
  vecf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
      cerr << "Wrong nr of AOs in vec file!\n";
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
  vecf.read((char *) cisvals, sizeof(double)*cis_size);
  vecf.read((char *) cisvecs, sizeof(double)*Cis_size*Cis_size);
  datf.write((char *) cisvals, Cis_size*(long long int) sizeof(double));
  datf.write((char *) cisvecs, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  double* cistmpvals = new double[Cis_size*Cis_size];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  if(x == y){
	      cistmpvals[x*cis_size+y] = cisvals[x];
	  }else{
	      cistmpvals[x*cis_size+y] = 0.;
	  }
      }
  }
trans_mat(cis_size, cistmpvals, cisvecs, cismat, 1);


  outf << "-------------------------------------------------------------------------------\n";

   //##CALCULATION OF IONIZATION RATES / WRITE OUT TO IRX-FILE / WITHOUT DOUBLES CORRECTIONS
  double* gamma_csf = new double[cis_size];
  int MO_r;
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      gamma_csf[x] = 0.;
      MO_r = (x-1)%umo+omo+llim;
      if(MOens[MO_r] > 0.) gamma_csf[x] = damp*sqrt(MOens[MO_r]);
  }
  syf.write((char *) gamma_csf, sizeof(double)*cis_size);

  double* energies = new double[cis_size];
  for(int x = 0; x < cis_size; x++) energies[x] = cisvals[x];

    outf << "-------------------------------------------------------------------------------\n";
 
  outf << "Building dipole transition matrices!\n";

  outf.flush();
  long int nonzero;
  cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dx,
	 mu_core[0], MOs, Pmat, tmpmat1, cistmpmat, &syf, nonzero, 2);
  double* dx = new double[cis_size*cis_size];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  dx[x*cis_size+y] = cismat[x*cis_size+y];
      }
  }
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dy,
	 mu_core[1], MOs, Pmat, tmpmat1, cistmpmat, &syf, nonzero, 2);
  double* dy = new double[cis_size*cis_size];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  dy[x*cis_size+y] = cismat[x*cis_size+y];
      }
  }
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dz,
	 mu_core[2], MOs, Pmat, tmpmat1, cistmpmat, &syf, nonzero, 2);
  double* dz = new double[cis_size*cis_size];
#pragma omp parallel for
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  dz[x*cis_size+y] = cismat[x*cis_size+y];
      }
  }
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  double* rates = new double[cis_size*cis_size];
  double* dephasing = new double[cis_size*cis_size];
  double* mu_tot = new double[cis_size*cis_size];
  double* omega = new double[cis_size*cis_size];
#pragma omp parallel for
  for(int i = 0; i < cis_size; i++){
        for(int j = 0; j < cis_size; j++){
            mu_tot[i*cis_size+j] = sqrt(pow(dx[i*cis_size+j],2.)+pow(dy[i*cis_size+j],2.)+pow(dz[i*cis_size+j],2.));
            omega[i*cis_size+j] = configurations[i] - configurations[j]; 
        }
  }
#pragma omp parallel for
  for(int i = 0; i < cis_size; i++){
        for(int j = 0; j < cis_size; j++){
            rates[i*cis_size+j] = (4.*pow(fabs(mu_tot[i*cis_size+j]),2.)*pow(fabs(omega[i*cis_size+j]),3.))/(3.*pow(137.036,3.));
        }
  }
  double sum = 0.;
  double dum = (10./1.7776791e+09)*(1/pow(0.2254,2.));
  for(int i = 0; i << cis_size; i++){
        for(int j = 0; j < cis_size; j++){
            sum = 0.;
            for(int k = 0; k < cis_size; k++){
                sum += (rates[i*cis_size+k]+rates[j*cis_size+k]);
            }
            dephasing[i*cis_size+j] = sum/2. + dum*omega[i*cis_size+j];
        }
  }
  syf.write((char *) rates, sizeof(double)*cis_size*cis_size);
  syf.write((char *) dephasing, sizeof(double)*cis_size*cis_size);
  syf.write((char *) energies, sizeof(double)*cis_size);
  syf.write((char *) cisvecs, sizeof(double)*cis_size*cis_size);
  syf.flush();

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Execution started on/at:\n";
  status(&outf);
  outf.flush();
}

void cp_mus(int cis_size, int nroao, int omo, int umo, int llim, int nroe, double* cismat, double* cisvecs, double* D,
	    double mu_core, double* MOs, double* Pmat, double* tmpvecs, double* cistmpmat,
	    ofstream* syf, long int &nonzero, int mode){
  long long int Cis_size = cis_size;
  if(mode == 0){
      int cpo = cis_size+1;
      long int* dip_row = new long int[cpo];
      nonzero = find_nonzeros(cis_size, umo, omo, llim);
      cout << "Nonzeroes: " << nonzero << "\n";
      long int* dip_col = new long int[nonzero];
      rows_cols_crs(cis_size, umo, omo, llim, dip_col, dip_row);
      syf->write((char *) dip_row, sizeof(long int)*cpo);
      syf->write((char *) dip_col, sizeof(long int)*nonzero);
  }
  if(mode == 1){
  double* dip_val = new double[nonzero];
  calc_mu_mat_cis_crs(cis_size, nroao, omo, umo, llim, nroe, cismat, D, mu_core, MOs, Pmat, tmpvecs, dip_val);
  syf->write((char *) dip_val, sizeof(double)*nonzero);
//  syf->write((char *) cismat, sizeof(double)*cis_size*cis_size);
  syf->flush();

  if(cis_size < 10000)
    trans_mat(cis_size, cismat, cisvecs, cistmpmat, 1); //Forward transfornation
  else
    trans_mat(Cis_size, cismat, cisvecs, cistmpmat, 1); //Forward transfornation
  }
  if(mode == 2){
      calc_mu_mat_cis(cis_size, nroao, omo, umo, llim, nroe, cismat, D, mu_core, MOs, Pmat, tmpvecs);
      syf->write((char *) cismat, sizeof(double)*cis_size*cis_size);
      syf->flush();
      if(cis_size < 10000){
	  trans_mat(cis_size, cismat, cisvecs, cistmpmat, 1);
      }else{
	  trans_mat(cis_size, cismat, cisvecs, cistmpmat, 1);
      }
  }
}

