/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: cis.cc                                                                 *
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


//MP2
extern double SO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);
extern double SO_calc_ijIIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints);
extern double SO_calc_ijIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints);
extern double SO_calc_a_ij_kl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints, double* MOens);
extern double SO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);
extern double SO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens);
extern double SO_calc_t1( int omo, int umo, int llim, double* cis_vec, double cis_en, double* prec_ints, double* MOens);
extern double SO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens);
extern double SO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
			     double* cis_vec, double cis_en, double* t1, double* t2);



//RO MP2 (speed test so far)
extern double RO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens);
extern double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double cis_en, double* t1, double* t2);
extern void SL_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
			   double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2, 
			   int cis_size, ofstream* outf);

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
  
  outf << "CIS [CIS3(D) Suite]\n";
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
  

  outf << "Building density and Fock matrix\n";
  
  build_Pmat_dscf( nroao, nroe,  Pmat,  Pmat_old,  MOs, 0.);
  build_Fmat(nroao,  Fmat, Pmat, Hmat, intval, intnums, sortcount, nrofint);

  double currE =  Calc_e_el( nroao, Fmat, Pmat, Hmat);
  double ekin = calc_op_1el( nroao, Tmat, Pmat);

  outf << "Final electronic energy: " << currE << "\n";
  outf << "Kinetic energy: " << ekin  << "\n";
  outf << "Total energy: " << currE+ion_rep << "\n\n";
  outf << "-V/T : " << -(currE+ion_rep-ekin)/ekin << "\n\n";

    

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
  

  outf << "Starting CIS calculation:\n";   
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
  outf << "Need " << (3*Cis_size*Cis_size+Cis_size)*sizeof(double) << " bytes of Memomry for CIS calculation\n";
  outf.flush();

  double*  cismat    = new double[Cis_size*Cis_size];     //Matrix for Operator representations in CIS-Basis
  double*  cistmpmat = new double[Cis_size*Cis_size];     //Temp space
  double*  cisvals   = new double[cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors
  double*  cisd_corr = new double[cis_size];

  //Singel loop D-CORR
  double*  C_t1      = new double[cis_size];
  double*  C_t2      = new double[cis_size];
  

  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Doing cis calculation\n";
  outf.flush();


  //## START EXTRA CODE FOR PRECALCULATED MOINTS
  double* prec_ints = precalc_ints_sd( omo,  umo,  llim,
				       nroao,  MOs,
				       nrofint,  sortcount, intval,
				       intnums, &outf);
  //## END OF EXTRA CODE
  
  //BUILD Hamiltionian in CIS-Basis 
  for(long long int x = 0; x < Cis_size*Cis_size; x++)
    cismat[x] = 0.;
  
  outf << "Calculating upper triangle of CIS-matrix\n";



  //// CODE FOR NORMAL LOOPS
  for(long long int x = 1 ; x < Cis_size; x++){
    int i1 = (x-1)/umo+llim;
    int f1 = (x-1)%umo+omo+llim;
    for(long long int y = x ; y < Cis_size; y++){
      int i2 = (y-1)/umo+llim;
      int f2 = (y-1)%umo+omo+llim;

      //## PRECALCULATED INTEGRAL --> biggest memory requierd
      cismat[x*Cis_size+y] +=   -get_precalc_ints_sd(i1, i2, f1, f2,
						     omo, umo, llim, prec_ints);
      cismat[x*Cis_size+y] += 2.*get_precalc_ints_sd(i1, f1, i2, f2,
						     omo, umo, llim, prec_ints);

      //## SEMI-DIRECT TRANSFORMATION
      //      cismat[x*cis_size+y] +=   -mo2int_op(i1, i2, f1, f2,
      // 					   nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
      //       cismat[x*cis_size+y] += 2.*mo2int_op(i1, f1, i2, f2,
      // 					   nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
      
      //## BRUT FORCE CALCULATIONS  --> lowest memory requierd
      //      cismat[x*cis_size+y] +=   -calc_mo2int_bf(i1, i2, f1, f2,
      // 						nroao, MOs, nrofint, sortcount, intval, intnums);
      //      cismat[x*cis_size+y] += 2.*calc_mo2int_bf(i1, f1, i2, f2,
      //  						nroao, MOs, nrofint, sortcount, intval, intnums);
      
      if(i1==i2 && f1==f2)
        cismat[x*Cis_size+y] += MOens[f1]-MOens[i1];
    }
    outf  << x << "\t";   
    if(x%10==0) outf << "\n";
    outf.flush();
  }


  //## ONLY MEANING FULL FOR SEMI DIRECT TRANS
  //Print out stat 
  // mo2int_op(-1, -1, -1, -1,
  // 	    nroao, MOs, nrofint, sortcount, intval, intnums, &outf);


  outf << "Done !!!\n";
  outf << "Calling diag!\n"; outf.flush();
  

  for(long long int x = 0; x < Cis_size; x++){
    for(long long int y = x; y < Cis_size; y++){
      cismat[y*Cis_size+x] = cismat[x*Cis_size+y];
    }
  }

  diag_mat(cis_size, cismat, cisvals, cisvecs);
  

  //WRITE DATA TO FILE
  datf.write((char *) cisvals, Cis_size*(long long int) sizeof(double));
  datf.write((char *) cisvecs, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  int nprint =  26; 
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
    sprintf(dumc,"%.8f\t%.8f,\t %.3i->%.3i\t%+.8f\n",cisvals[x],cisvals[x]*27.211383,mi,mf,max_coeff/sqrt(2.));
    outf << x << "          "<< dumc;
  }

  outf << "-------------------------------------------------------------------------------\n";




  
  if(calc_d == 1) {
    //## PERTURBATION THEORY PART -> ONLY IMPLEMENTED FOR PRECALCULATED INTEGRALS SOFAR
    outf << "Calculating doubles corrections\n";
    for(int x = 0; x < cis_size; x++) cisd_corr[x] = 0.;

    double E2=   RO_calc_MP2( omo,  umo,  llim,  prec_ints,  MOens);

    outf << "Ground state MP2 energy (RO-OPS):  T2 " << E2 << " TOT "<< E2+currE+ion_rep<< " \n";  
    outf << "...............................................................................\n";
   
    //SINGLE LOOP CODE 
    SL_calc_cis_dc( omo,  umo,  llim, prec_ints,  MOens,
		    cisvecs, cisvals,  cistmpmat, C_t1, C_t2, 
		    cis_size, &outf);
    
    for(int x = 0; x < cis_size; x++) cisd_corr[x] = C_t1[x]+C_t2[x];
    
    for(int x = 1; x < nprint; x++){
      outf << "State " << x << " : \n";
      outf << "doubles correction is " << C_t1[x]  << "\n";
      outf << "triples correction is " << C_t2[x] << "\n";
      outf << "Total corr " <<  cisd_corr[x]  << " ex en [Eh,ev] " << cisvals[x]+cisd_corr[x]  <<  " "<< (cisvals[x]+cisd_corr[x] )*27.211383 << "\n";
      outf << "...............................................................................\n";
    }
    outf.flush();

    //OLD D CORRECTIONS PART
//     outf << "D corrections for the first excited states !!\n";
//     outf.flush();
//     double T1_corr, T2_corr;
//     for(int x = 1; x < nprint; x++){
//       outf << "State " << x << " : \n";
//       double tot_corr = RO_calc_cis_dc (omo,  umo,  llim,  prec_ints, MOens,
// 					&(cisvecs[cis_size*x]), cisvals[x], &T1_corr, &T2_corr);
      
//       outf << "doubles correction is " << T1_corr << "\n";
//       outf << "triples correction is " << T2_corr << "\n";
//       outf << "Total corr " <<  tot_corr << " ex en [Eh,ev] " << cisvals[x]+tot_corr <<  " "<< (cisvals[x]+tot_corr)*27.211383 << "\n";
//       cisd_corr[x] = tot_corr;
//       outf << "...............................................................................\n";
//       outf.flush();
//     }


//     outf << "Calculating remaining corrections\n";
//     outf.flush();
//     for(int x = nprint; x < cis_size; x ++){
//       cisd_corr[x]  = RO_calc_cis_dc (omo,  umo,  llim,  prec_ints, MOens,
// 				      &(cisvecs[cis_size*x]), cisvals[x], &T1_corr, &T2_corr);
      
//       outf  << x << "\t";   
//       if((x-nprint+1)%10==0) outf << "\n";
//       outf.flush();
//     }
//     outf << "\n";
    outf << "-------------------------------------------------------------------------------\n";
    //## END OF PT-2nd Order part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  }
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


  if(calc_d == 1) {
    //## write out PERTURBATION THEORY PART
    outf << "\n\nWriting PT-corrections to bcs-file\n";
    
    datf.write((char *) cisd_corr, cis_size*sizeof(double));
    datf.flush();
  }

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Execution started on/at:\n";
  status(&outf);
  outf.flush();
}

void cp_mus(int cis_size, int nroao, int omo, int umo, int llim, int nroe, double* cismat, double* cisvecs, double* D,
	    double mu_core, double* MOs, double* Pmat, double* tmpvecs, double* cistmpmat,
	    ofstream *outf, int nprint){
  char dumc[512];
  long long int Cis_size = cis_size;
  

  outf->flush();
  calc_mu_mat_cis( cis_size,  nroao, omo,  umo, llim, nroe, cismat,  D,  mu_core,  MOs,  Pmat,  tmpvecs);
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
