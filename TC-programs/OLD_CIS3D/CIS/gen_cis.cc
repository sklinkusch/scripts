/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: gen_cis.cc                                                             *
 *                                                                              *
 * RHF program                                                                  *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
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
	    ofstream *outf, int nprint);


//Extern Functions
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void status(ofstream* outf);
extern void get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void read_sys(char* sysfile, double* coord, double* charges, double* mass, 
		     double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
		     double *Dz, long long int* sortcount, double* intval, 
		     unsigned short* intnums);
extern double calc_ion_rep(int nroa, double* coord, double* charges);
extern void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
			 double* mu_core);
extern double calc_S12(int nroao, double* Smat, double* Som12, double* tmpmat, 
		       double* tmpvecs, double* tmpvals);
extern void diag_Fmat(int nroao, double* Fmat, double* MOs, 
		      double* MOens, double* Som12, double* tmpmat);
extern double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old, 
			      double* MOs, double damp);
extern void build_Fmat(int nroao, double* Fmat, double* Pmat, double* Hmat,
		       double* intval, unsigned short* intnums,
		       long long int* sortcount, long long int nrofint);
extern double Calc_e_el(int nroao, double* Fmat, double* Pmat, double* Hmat);
extern double calc_op_1el(int nroao, double* opmat, double* Pmat);
extern void write_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void  diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
			    double *cismat, double* mumat, double mucore, 
			    double* MOs, double* Pmat, double* tmpvec_ao);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
extern void pmv(double* mat, double* vi, double* vo, int nroao);

//2EL INTEGRALS
extern double calc_mo2int_bf(int i, int j, int k, int l, int nroao, double* MOs,
			     long long int nrofint, long long int* sortcount, double* intval,
			     unsigned short* intnums);
extern double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
			long long int nrofint, long long int* sortcount,  double* intval,
			unsigned short* intnums, ofstream* outf);
extern double* precalc_ints_sd(int omo, int umo, int llim,
			       int nroao, double* MOs,
			       long long int nrofint, long long int* sortcount,  double* intval,
			       unsigned short* intnums, ofstream* outf);
extern double get_precalc_ints_sd(int i, int j, int k, int l,
				  int omo, int umo, int llim, double* prec_ints);



int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  int     nroe;           //Nr of electrons  (if negative, read in center of mass)
  int     llim;           //first MO used for correlation
  int     ulim;           //last  MO used for correlation
  char sysfile[256];      //binary system file
  char wavfile[256];      //HF-Wavefunction file
  int   calc_d=0;         //if 1 -> calc duplett correction PT
  int calc_COM = 0;       // 0 -> calc COM otherwise read in

  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "GCIS [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  
  ist >> sysfile >> nroe;

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
  
  ist >> llim >> ulim >> wavfile >> calc_d;
  outf << "System data \nNr of electrons: " << nroe << "\nsys-file: " << sysfile << "\n\n";
  outf << "CIS data\nLower limit of MOs used for correlation: " << llim << "\n";
  outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
  outf << "Reading HF-wave function from  " << wavfile << "\n";
  if(calc_d == 1) 
    outf << "Calculating 2nd order PT corrections [MP2, CIS(D)]\n";
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

  //atoms
  double*        coord;         //atomic coordinats               3*nroa
  double*        charges;       //atomic charges                    nroa
  double*        mass;          //atomic masses                     nroa

  //one electron mat&vecs
  double*        Smat;          //Overlap matrix S                  nroao*nroao
  double*        Hmat;          //one electron Hamiltionian         nroao*nroao
  double*        Tmat;          //Kinetic energy operator           nroao*nroao
  
  double*        Som12;         //S^-1/2                            nroao*nroao
  double*        Pmat;          //density matrix                    nroao*nroao
  double*        Pmat_old;      //Old density matrix                nroao*nroao

  double*        Fmat;          //Fock matrix                       nroao*nroao
  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
  double*        MOens_old;     //Old MO energies;                  nroao


  //Temporary memory spaces
  double*  tmpmat1;                //                               nroao*nroao
  double*  tmpmat2;                //                               nroao*nroao
  double*  tmpvecs;                //                               nroao*nroao

  double*  tmpvals;                //                               nroao
 
  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+14*nroao*nroao+3*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]); inc += 3*nroa; charges = &(dumd[inc]); inc += nroa; mass  = &(dumd[inc]); inc += nroa;
  Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao; Tmat = &(dumd[inc]); inc += nroao*nroao;
  Som12 = &(dumd[inc]); inc += nroao*nroao; Pmat = &(dumd[inc]); inc += nroao*nroao;Pmat_old = &(dumd[inc]); inc += nroao*nroao;
  Fmat = &(dumd[inc]); inc += nroao*nroao; MOs = &(dumd[inc]); inc += nroao*nroao; Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao; Dz = &(dumd[inc]); inc += nroao*nroao; 
  
  MOens = &(dumd[inc]); inc += nroao; MOens_old = &(dumd[inc]); inc += nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; tmpmat2 = &(dumd[inc]); inc += nroao*nroao; tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
  tmpvals = &(dumd[inc]); inc += nroao;

  //MEMORY ALLOCATION for two electron values 
  outf << "Need " << nrofint*(sizeof(double)+sizeof(unsigned short)*4) << " bytes for two electron data\n";
  outf.flush();
  
  double*         intval        = new double[nrofint];                     //two electron integrals
  unsigned short* intnums       = new unsigned short[nrofint*4];           //two electron indices
  long long int   sortcount[4];                                            //num of two electron Integrals in each perm. type

  
  //Read system data 
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, Smat, Dx, Dy,  Dz, sortcount, intval, intnums);
  

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
  double ion_rep =  calc_ion_rep( nroa, coord, charges);
  
  if(calc_COM == 0)
    calc_center_of_mass( nroa, coord,  mass,  center_of_mass);
  calc_mu_core( nroa, coord, charges, center_of_mass, mu_core);
  
  outf << "\n\nIon repulsion is: " << ion_rep  << "\n";
  outf << "Center of mass (x,y,z): " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << "\n";
  outf << "Core dipole moment (at C.o.M.): " << mu_core[0] << " " << mu_core[1] << " " << mu_core[2] << "\n";
  outf << "...............................................................................\n";
  outf << "Init HF operators:\n";
  
  outf << "Calculating S^-1/2\n";
  outf << "Minimal eigenvalue of S is " << calc_S12( nroao, Smat, Som12, tmpmat1, tmpvecs, tmpvals) <<"\n";

  read_wav_HF(wavfile, nroao, MOens, MOs);
  outf << "HF-wave function  read from " << wavfile << "\n";

  outf << "Testing for orthonormal MOs\n";
  for(int x = 0; x < nroao; x++)
    pmv(Smat, &(MOs[x*nroao]), &(tmpvecs[x*nroao]), nroao);
  
  double max_err = 0.;
  int mx, my;
  for(int x = 0; x < nroao; x++){
    for(int y = x+1; y < nroao; y++){
      double ov = 0.;
      for(int z = 0; z < nroao; z++) ov += MOs[x*nroao+z]*tmpvecs[y*nroao+z];
      if(fabs(ov) > max_err){
	max_err = fabs(ov);
	mx = x;
	my = y;
      }
    }
  }
  outf << "Maximal overlapp is " << max_err << " for " << mx << " , " << my << "\n";

  max_err = 0.;
  for(int x = 0; x < nroao; x++){
    double ov = 0.;
    for(int z = 0; z < nroao; z++) ov += MOs[x*nroao+z]*tmpvecs[x*nroao+z];
    if(fabs(ov-1.) > max_err) max_err = fabs(ov-1.);
  }
  outf << "Maximal norm deviation is " << max_err << "\n";

  outf << "Calculating one particle energies\n";
  for(int x = 0; x < nroao; x++){
    pmv(Hmat, &(MOs[x*nroao]), &(tmpvecs[x*nroao]), nroao);
    MOens_old[x] = 0.;
    for(int z = 0; z < nroao; z++) MOens_old[x] += MOs[x*nroao+z]*tmpvecs[x*nroao+z];
  }
  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++){
      tmpmat2[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++) tmpmat2[x*nroao+y] += MOs[x*nroao+z]*tmpvecs[y*nroao+z];
    }
  }
  double* e1pmat = tmpmat2;   //of diagonals of h1
  
  outf << "-------------------------------------------------------------------------------\n";
  

  outf << "Starting Gen-CIS calculation:\n";   
  sprintf(dumc,"%s.bcs",argv[2]);
  outf << "Writing binary data to " << dumc << "\n";
  ofstream datf(dumc);
  

  outf << "\n\n";
  

  //CHECK limits
  if(llim < 0)      { outf << "Invalid lower limit in MO-range!\n"; outf.flush(); exit(4);}
  if(ulim >= nroao) { outf << "Invalid upper limit in MO-range!\n"; outf.flush();exit(4);}
  if(llim >= nroe/2){ outf << "No occupied orbitals in MO-range!\n"; outf.flush(); exit(4);}

  int nrof = ulim - llim +1;
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  int cis_size = omo*umo+1;
  long long int Cis_size = cis_size;
  outf << "nr of CSF               : " << cis_size << "\n";


  //Write data to data-file
  datf.write((char *) &nroao, sizeof(int));
  datf.write((char *) &nroe, sizeof(int));
  datf.write((char *) &llim, sizeof(int));
  datf.write((char *) &ulim, sizeof(int));
  

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Allocating Memory for CIS:\n";
 
  //CIS Variables
  outf << "Need " << (3*cis_size*cis_size+cis_size)*sizeof(double) << " bytes of Memomry for CIS calculation\n";
  outf.flush();

  double*  cismat    = new double[Cis_size*Cis_size];     //Matrix for Operator representations in CIS-Basis
  double*  cistmpmat = new double[Cis_size*Cis_size];     //Temp space
  double*  cisvals   = new double[cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors

  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Doing cis calculation\n";
  outf.flush();

  //CALC MORE INTEGRALS

  //## START EXTRA CODE FOR PRECALCULATED MOINTS
  // we need <(0-nroe/2)(0-nroao)|(0-nroao)(0-nroao)
  outf << "Precalculating <oo|oo> up to <ov|vv>\n";
  int prec_mem = nroe/2*nroao*nroao*nroao;
  outf << "Need " << prec_mem*sizeof(double) << " bytes (" << prec_mem*sizeof(double)/1024/1024 << " MB) for precalculation\n";
  double* prec_ints =  new double[prec_mem];
    
  long long int prec_count = 0;

  int kstep = nroao;
  int jstep = nroao*kstep;
  int istep = nroao*jstep;

  for(int i = 0; i < nroe/2; i++){
    for(int j = 0; j < nroao; j++){
      for(int k = 0; k < nroao; k++){
        for(int l = 0; l <  nroao;l++){
          prec_ints[prec_count++] = mo2int_op(i, j, k, l,
                                              nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
        }
      }
    }
    outf << i << "\t";   
    if((i+1)%10==0) outf << "\n";
    outf.flush();
  }
  outf << "\nInts claculated: " << prec_count  << "\n";
  outf << "\nDone\n";
  //Print out stat 
  mo2int_op(-1, -1, -1, -1,
            nroao, MOs, nrofint, sortcount, intval, intnums, &outf); 
  

  //## END OF EXTRA CODE
  
  //BUILD Hamiltionian in CIS-Basis 
  for(int x = 0; x < cis_size*cis_size; x++)
    cismat[x] = 0.;
  
  outf << "Calculating upper triangle of CIS-matrix\n";

  //// CODE FOR GENERAL LOOP FOR GS-EXCITED MATRIX ELEMENTS
  
  //DIAGONAL ELEMENTS 0-0
  for(int x = 0; x < nroe/2; x++) cismat[0] += 2.*MOens_old[x];
  for(int x = 0; x < nroe/2; x++){
    for(int y = 0; y < nroe/2; y++)
      cismat[0] += (2.*prec_ints[x *istep+ x *jstep+ y *kstep+ y ] 
	            -  prec_ints[x *istep+ y *jstep+ x *kstep+ y ]);
  }
  cismat[0] += ion_rep;


  //DIAGONAL ELEMENTS EX-EX
  for(int x = 1; x < cis_size; x++){
    //V_kk (ion rep)
    cismat[x*cis_size+x] += ion_rep;
    int i1 = (x-1)/umo+llim;
    int f1 = (x-1)%umo+omo+llim;
    
    //one particle energie
    for(int y = 0; y < nroe/2; y++)
      cismat[x*cis_size+x] += 2.*MOens_old[y];
    cismat[x*cis_size+x] -=    MOens_old[i1];
    cismat[x*cis_size+x] +=    MOens_old[f1];
    //e-e interaction
    for(int C = 0; C < nroe/2; C++){
      if(C != i1){
	for(int D = 0; D < nroe/2; D++){
	  if(D != i1){
	    cismat[x*cis_size+x] += (2.*prec_ints[C *istep+ C *jstep+ D *kstep+ D ] 
				     -  prec_ints[C *istep+ D *jstep+ C *kstep+ D ] );
	  }
	}
      
	cismat[x*cis_size+x] += (2.*prec_ints[C *istep+ C  *jstep+ i1 *kstep+ i1 ] 
				 -  prec_ints[C *istep+ i1 *jstep+ C  *kstep+ i1 ] );
      
	cismat[x*cis_size+x] += (2.*prec_ints[C *istep+ C  *jstep+ f1 *kstep+ f1 ] 
				 -  prec_ints[C *istep+ f1 *jstep+ C  *kstep+ f1 ] );
      }
      
    }
    cismat[x*cis_size+x] +=         prec_ints[i1 *istep+ i1 *jstep+ f1  *kstep+ f1 ];   
    cismat[x*cis_size+x] +=         prec_ints[i1 *istep+ f1 *jstep+ i1  *kstep+ f1 ];  
  }

  //OFF-DIAGONAL ELEMENTS 0-EX
  for(int x = 1; x < cis_size; x++){
    int i1 = (x-1)/umo+llim;
    int f1 = (x-1)%umo+omo+llim;
    cismat[x] += -sqrt(2.)*e1pmat[i1*nroao+f1];
    for(int C = 0; C < nroe/2; C++){
	cismat[x] += -sqrt(2.)* (2.*prec_ints[C *istep+ C  *jstep+ i1 *kstep+ f1 ] 
				 -  prec_ints[C *istep+ i1 *jstep+ C  *kstep+ f1 ] );
    }
  }
  

  ////   //OFF-DIAGONAL ELEMENTS EX-EX
  for(int x = 1 ; x < cis_size; x++){
    int i1 = (x-1)/umo+llim; 
    int f1 = (x-1)%umo+omo+llim; 
    for(int y = x +1; y < cis_size; y++){
      int i2 = (y-1)/umo+llim;
      int f2 = (y-1)%umo+omo+llim;
      //case IV
      if(i1 == i2){
	cismat[x*Cis_size+y] += e1pmat[f1*nroao+f2];
	for(int C = 0; C < nroe/2; C++){
	  if(C != i1)
	    cismat[x*Cis_size+y] += (2.*prec_ints[C *istep+ C  *jstep+ f1 *kstep+ f2 ] 
				     - prec_ints[C *istep+ f2  *jstep+ C *kstep+ f1 ] );
	}
	cismat[x*Cis_size+y] += prec_ints[i1 *istep+ i1  *jstep+ f1 *kstep+ f2 ] + prec_ints[i1 *istep+ f2  *jstep+ i1 *kstep+ f1 ];
      }
      // case V
      if(f1 == f2){
	cismat[x*Cis_size+y] -= e1pmat[i1*nroao+i2];
	for(int C = 0; C < nroe/2; C++){
	  if(C != i1 && C != i2)
	    cismat[x*Cis_size+y] -= (2.*prec_ints[C *istep+ C  *jstep+ i1 *kstep+ i2 ] 
				     - prec_ints[C *istep+ i1  *jstep+ C *kstep+ i2 ] );
	}
	cismat[x*Cis_size+y] += (-prec_ints[i1 *istep+ i2  *jstep+ f1 *kstep+ f1 ] 
				 +prec_ints[i1 *istep+ f1  *jstep+ f1 *kstep+ i2 ] 
				 +prec_ints[i2 *istep+ f1  *jstep+ f1 *kstep+ i1 ]
				 -prec_ints[i1 *istep+ i2  *jstep+ i1 *kstep+ i1 ]
				 -prec_ints[i1 *istep+ i2  *jstep+ i2 *kstep+ i2 ]); 
      }
      //case VI
      if(i1 != i2 && f1 != f2){
	cismat[x*Cis_size+y] +=   -prec_ints[i1 *istep+ i2 *jstep+ f1 *kstep+ f2 ];
	cismat[x*Cis_size+y] += 2.*prec_ints[i1 *istep+ f1 *jstep+ i2  *kstep+ f2];
      }
    }
    outf  << x << "\t";   
    if(x%10==0) outf << "\n";
    outf.flush();
  }

  


  outf << "Done !!!\n";
  outf << "Calling diag!\n"; outf.flush();
  

  for(int x = 0; x < cis_size; x++){
    for(int y = x; y < cis_size; y++){
      cismat[y*cis_size+x] = cismat[x*cis_size+y];
    }
  }

  diag_mat(cis_size, cismat, cisvals, cisvecs);

  outf << "Ground state energy is: "<< cisvals[0] << "\n";
  

  //WRITE DATA TO FILE
  datf.write((char *) cisvals, Cis_size*(long long int) sizeof(double));
  datf.write((char *) cisvecs, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();
  

  int nprint =  8; 
  if(nprint > cis_size) nprint = cis_size;

  outf << "First " << nprint-1 << " exc. states\n [energies:  (Hartree) (eV),  (LC  x->y)  coeff :\n";
  outf << "...............................................................................\n";
  
  for(int x = 1; x < nprint; x++){
    double max_coeff = 0.;
    int mi = 0, mf = 0;
    for(int y = 0; y < cis_size; y++){
      if( fabs(cisvecs[x*cis_size+y]) > fabs(max_coeff)){
	max_coeff = cisvecs[x*cis_size+y];
	mi = (y-1)/umo+llim;
	mf = (y-1)%umo+omo+llim;
      }
    }
    sprintf(dumc,"%.8f\t%.8f,\t %.3i->%.3i\t%+.8f\n",(cisvals[x]-cisvals[0]),(cisvals[x]-cisvals[0])*27.211383,mi,mf,max_coeff/sqrt(2.));
    outf << x << "          "<< dumc;
  }

  outf << "-------------------------------------------------------------------------------\n";




  
  outf << "Building dipole transition marticies!\n";

   outf << "...............................................................................\n";
   outf << "X\n";
   outf.flush();
   cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dx,
	  mu_core[0], MOs, Pmat, tmpmat1, cistmpmat, &outf, nprint);
   
  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  outf << "Values ranging from " << cisvals[0] << " to " << cisvals[cis_size-1] << "\n";
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  outf << "...............................................................................\n";
  outf << "Y\n";
  outf.flush();
  cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dy,
	 mu_core[1], MOs, Pmat, tmpmat1, cistmpmat, &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  outf << "Values ranging from " << cisvals[0] << " to " << cisvals[cis_size-1] << "\n";
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();


  outf << "...............................................................................\n";
  outf << "Z\n";
  outf.flush();
  cp_mus(cis_size,  nroao,  omo,  umo, llim, nroe,  cismat,  cisvecs, Dz,
	 mu_core[2], MOs, Pmat, tmpmat1, cistmpmat, &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, cismat, cisvals, cistmpmat);
  outf << "Values ranging from " << cisvals[0] << " to " << cisvals[cis_size-1] << "\n";
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();



  outf << "-------------------------------------------------------------------------------\n";
  outf << "Execution ended on/at:\n";
  status(&outf);
  outf.flush();
}

void cp_mus(int cis_size, int nroao, int omo, int umo, int llim, int nroe, double* cismat, double* cisvecs, double* D,
	    double mu_core, double* MOs, double* Pmat, double* tmpvecs, double* cistmpmat,
	    ofstream *outf, int nprint){
  char dumc[512];
  long long int Cis_size = cis_size;
  

  outf->flush();
  *outf << "Calculating 1p matrix elements\n";

  //ERASE 1P MATRIX
  for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;

  for(int x = 0; x < nroao; x++)
    pmv(D, &(MOs[x*nroao]), &(tmpvecs[x*nroao]), nroao);
  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++){
      for(int z = 0; z < nroao; z++) Pmat[x*nroao+y] += -MOs[x*nroao+z]*tmpvecs[y*nroao+z]; //negative position !!
    }
  }
  
  //ERASE MATRIX
  for(int x = 0; x < cis_size*cis_size; x++) cismat[x] = 0.;
  //ADD CORE DIPOLE
  for(int x = 0; x < cis_size; x++) cismat[x*cis_size+x] +=  mu_core;

  //0-0
  for(int C = 0; C < nroe/2; C++)  cismat[0] += 2.*Pmat[C*nroao+C];
  
  //0-EX
  for(int x = 1; x < cis_size; x++){
    int i1 = (x-1)/umo+llim;
    int f1 = (x-1)%umo+omo+llim;
    cismat[x] += -sqrt(2.)*Pmat[i1*nroao+f1];
  }

  //EX-EX DIAG
  for(int x = 1; x < cis_size; x++){
    int i1 = (x-1)/umo+llim;
    int f1 = (x-1)%umo+omo+llim;
    for(int y = 0; y < nroe/2; y++)
      cismat[x*cis_size+x] +=  2.*Pmat[y*nroao+y];
    cismat[x*cis_size+x] -=     Pmat[i1*nroao+i1];
    cismat[x*cis_size+x] +=     Pmat[f1*nroao+f1];
  }
  
  ////   //OFF-DIAGONAL ELEMENTS EX-EX
  for(int x = 1 ; x < cis_size; x++){
    int i1 = (x-1)/umo+llim; 
    int f1 = (x-1)%umo+omo+llim; 
    for(int y = x +1; y < cis_size; y++){
      int i2 = (y-1)/umo+llim;
      int f2 = (y-1)%umo+omo+llim;
      //case IV
      if(i1 == i2){
	cismat[x*Cis_size+y] +=  Pmat[f1*nroao+f2];
      }
      // case V
      if(f1 == f2){
	cismat[x*Cis_size+y] -=  Pmat[i1*nroao+i2];
      }
    }
  }
  
  for(int x = 0; x < cis_size; x++){
    for(int y = x; y < cis_size; y++){
      cismat[y*cis_size+x] = cismat[x*cis_size+y];
    }
  }

  *outf << "Transforming Dipole matrix\n";
  outf -> flush();

  if(cis_size < 10000)
    trans_mat(cis_size, cismat, cisvecs, cistmpmat, 1); //Forward transfornation
  else
    trans_mat(Cis_size, cismat, cisvecs, cistmpmat, 1); //Forward transfornation
   

  *outf << "Transition Moments\n";
  for(int x = 0; x < nprint; x++) *outf << "        " << x;
  *outf << "\n";
  
  for(long long int x = 0; x < nprint ; x++ ){
    char tmpc[12];
    sprintf(dumc,"%i   ",(int) x);
    for(long long int y = 0; y < nprint; y++) {
      sprintf(tmpc,"%+.4f  ", cismat[x*cis_size+y]);
      strcat(dumc,tmpc);
    }
    *outf << dumc << "\n";
  }
  *outf << "\n\n";
  outf->flush();
}
