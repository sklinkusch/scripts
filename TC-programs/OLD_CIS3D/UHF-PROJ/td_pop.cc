/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: td_pop.cc                                                              *
 *                                                                              *
 * td  population analysis   program                                            *
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
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
#include <complex>
#include <stdlib.h>

using namespace std;
#define Complex complex<double>


//Functions
void check_range(int tomuch, int current_number, ofstream *outf);

//Extern Functions
extern int    rem_com(char* filename, char* streamstring, int string_length);
extern void   status(ofstream* outf);
extern void   get_sys_size(char* sysfile, int* nroao, int* nroa, long long int* nrofint);
extern void   read_sys_1el(char* sysfile, double* coord, double* charges, double* mass, 
			   double* Hmat, double* Tmat, double* Smat,  double* Dx, double* Dy, 
			   double *Dz);
extern void   read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern double calc_ion_rep(int nroa, double* coord, double* charges);
extern void   calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void   calc_mu_core(int nroa, double* coord, double* charges, double* point, 
			   double* mu_core);
extern double calc_op_1el(int nroao, double* opmat, double* Pmat);
extern void   transformationCIS_SE(Complex* wav, Complex* SEwav, double* cisvecs, int cis_size);
extern void transformationCIS_SE(Complex* wav, Complex* SEwav, double* cisvecs, long long int cis_size); //for large matrices
extern void   calc_MOpops(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* pops);
extern void   calc_dens(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double* Pmat, double* MOs);
extern void   mat_mat(int np, double* mat_i1, double* mat_i2, double* mat_f);
extern void   output_matrix(double* mat,int nr_of_col, int nrop, ofstream *outf);

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "TDPOP [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";

  //INPUT VARIABLES
  
  //SYSTEM
  char    bcsfile[256];         //binary system file (CIS space)
  char    sysfile[256];         //binary system file (AO  space)
  char    hfwfile[256];         //Hartree-Fock wave function file
  char    tdwfile[256];         //time dependent wave function file
  int     use_d = 0;            //Use D correstion (0 false 1 true)
  
  //ANALYSIS
  int     nrots;                //Nr of time steps to analyse
  //print our intervalls
  int     main_mw;              //Analyse every main_mw-th time step
  int     mw_exp;               //Print out intervall for expectation values E(exc), T, mu_x, mu_y, mu_z, Mulliken charges and bond orders 
  int     mw_mopop;             //Print out intervall for MO-populations
  int     mw_mat;               //Print out intervall for matricies (desnity, full Mulliken)

  //population analysis
  int     nroc;                 //nr of mulliken charges to analyse
  int     nrob;                 //nr of bonds to analyse
  int*    pop_centers;          //wich charges to compute
  int*    pop_bonds_atom1;      //first atoms of bonds to compute        
  int*    pop_bonds_atom2;      //last  atoms of bonds to compute
  int*    first_basis;          //first basis function on a center
  int*    last_basis;           //last  basis function on a center
  int     nrof_basis_input;     //nr of basis ranges defined in the input file
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  
  ist >> bcsfile >> sysfile >> hfwfile >> tdwfile >> use_d;
  outf << "Binary system file (CIS space):  " << bcsfile << "\n";
  outf << "Binary system file (AO  space):  " << sysfile << "\n";
  outf << "HF wave function file:           " << hfwfile << "\n";
  outf << "TD-CIS wav file:                 " << tdwfile << "\n";
  if(use_d == 1)
    outf << "USING CIS(D) energy corrections.\n";
  else
    outf << "NO CIS(D) energy corrections.\n";
    
  outf << "\nReading HF-data.\n";
  
  //system size
  int            nroao;         //Nr of basis functions
  int            nroa;          //Nr of atoms 
  long long  int nrofint;       //Nr of two electron Integrals (not used here)
  
  get_sys_size( sysfile, &nroao, &nroa,  &nrofint);

  outf << "System sizes read from " << sysfile << "\n";
  outf << "Nr of basis functions: " << nroao   << "\n";
  outf << "Nr of atoms:           " << nroa << "\n";
  outf << "Allocating Memory\n";

  first_basis = new  int[nroa];         
  last_basis  = new  int[nroa];         
  for(int x = 0; x < nroa; x++)
    first_basis[x] = last_basis[x] = 0;
  
  
  //atoms
  double*        coord;         //atomic coordinats               3*nroa
  double*        charges;       //atomic charges                    nroa
  double*        mass;          //atomic masses                     nroa

  //one electron mat&vecs
  double*        Smat;          //Overlap matrix S                  nroao*nroao
  double*        Hmat;          //one electron Hamiltionian         nroao*nroao
  double*        Tmat;          //Kinetic energy operator           nroao*nroao
  
  double*        Pmat;          //density matrix                    nroao*nroao
  double*        Pmat_Mul;      //Mulliken matrix                   nroao*nroao

  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
  double*        tdPOP_mo;      //td MO-populations                 nroao

  //Temporary memory spaces
  double*  tmpmat1;             //                                  nroao*nroao
  double*  tmpmat2;             //                                  nroao*nroao
  double*  tmpvecs;             //                                  nroao*nroao

  double*  tmpvals;             //                                  nroao

  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+12*nroao*nroao+3*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]); inc += 3*nroa;     charges = &(dumd[inc]); inc += nroa;      mass  = &(dumd[inc]); inc += nroa;
  Smat = &(dumd[inc]); inc += nroao*nroao; Hmat = &(dumd[inc]); inc += nroao*nroao;  Tmat = &(dumd[inc]); inc += nroao*nroao;
  Pmat = &(dumd[inc]); inc += nroao*nroao; Pmat_Mul = &(dumd[inc]); inc += nroao*nroao;
  MOs = &(dumd[inc]); inc += nroao*nroao;  Dx = &(dumd[inc]); inc += nroao*nroao;
  Dy = &(dumd[inc]); inc += nroao*nroao;   Dz = &(dumd[inc]); inc += nroao*nroao; 
  
  MOens = &(dumd[inc]); inc += nroao;      tdPOP_mo = &(dumd[inc]); inc += nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; 
  tmpmat2 = &(dumd[inc]); inc += nroao*nroao; 
  tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
  tmpvals = &(dumd[inc]); inc += nroao;
 
  read_sys_1el(sysfile,  coord, charges, mass, Hmat,  Tmat,  Smat,   Dx,  Dy,  Dz);
  read_wav_HF(hfwfile, nroao, MOens, MOs);

 
  outf << "HF reference data:\n";
  outf << "...............................................................................\n";
  outf << "Ion Cores: coordinates mass charge\n";
  int count = 0;
  for(int x = 0; x < nroa; x++){
    outf << x << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%+.5f", coord[count++]);    outf  << dumc << "\t";
    sprintf(dumc,"%.5f", mass[x]);            outf  << dumc << "\t";
    sprintf(dumc,"%.0f", charges[x]);         outf  << dumc << "\n";     
  }
  
  outf.precision(10);

  double center_of_mass[3]; 
  double mu_core[3];
  double ion_rep =  calc_ion_rep( nroa, coord, charges);

  calc_center_of_mass( nroa, coord,  mass,  center_of_mass);
  calc_mu_core( nroa, coord, charges, center_of_mass, mu_core);
  
  outf << "\n\nIon repulsion is: " << ion_rep  << "\n";
  outf << "Center of mass (x,y,z): " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << "\n";
  outf << "Core dipole moment (at C.o.M.): " << mu_core[0] << " " << mu_core[1] << " " << mu_core[2] << "\n";
  outf << "\nMO-energies:\n";
  for(int x = 0; x < nroao; x++){
    sprintf(dumc,"%+.5f ",MOens[x]);
    outf << dumc << "\t";
    if((x+1)%5==0) outf << "\n";
  }
  outf << "...............................................................................\n";

  
  ist  >>  nrots >> main_mw >> mw_exp >> mw_mopop >> mw_mat;
  outf << "\nPRINT OUT INTERVALS\n";
  outf << "Main analysis intervall:         " << main_mw << "\n";
  outf << "Expectation values:              " << mw_exp  << "\n";
  outf << "MO-populations:                  " << mw_mopop<< "\n";
  outf << "Density Matrix:                  " << mw_mat  << "\n";

  outf << "\nANALYSIS DETAILS\n";
  ist  >> nroc;
  pop_centers = new int[nroc];
  outf << "Charges (Mulliken) computed for centers:\n";
  for(int x = 0; x < nroc; x++){
    ist >> pop_centers[x];
    check_range(nroa ,pop_centers[x], &outf);
    outf << "\t" << pop_centers[x] << "\n";
  }
  
  ist >> nrob;
  pop_bonds_atom1 = new int[nrob];
  pop_bonds_atom2 = new int[nrob];
  outf << "\nBond orders (according to I. Mayer) computed for bonds between:\n";
  for(int x = 0; x < nrob; x++){
    ist >> pop_bonds_atom1[x] >> pop_bonds_atom2[x];
    check_range(nroa ,pop_bonds_atom1[x], &outf);
    check_range(nroa ,pop_bonds_atom2[x], &outf);
    outf << "\t" <<pop_bonds_atom1[x] << "-" << pop_bonds_atom2[x] << "\n";
  }

  outf << "\nBASIS FUNCTION DISTRIBUTION\n";
  ist >> nrof_basis_input;
  outf << "Reading basis localizytion for " << nrof_basis_input << " input\n";
  for(int x = 0; x < nrof_basis_input; x++){
    int center;
    ist >> center;
    ist >> first_basis[center] >> last_basis[center];
    check_range(nroa, center, &outf);
    check_range(nroao, first_basis[center], &outf);
    check_range(nroao, last_basis[center],  &outf);
    outf << "Center: " << center << ", first: " << first_basis[center] << ", last: " << last_basis[center] << "\n";
  }
  
  
  //CIS-VARIABLES
  outf << "\nREADING CIS-SPACE\n";
  int     nroe;                 //Nr of electrons
  int     llim;                 //first MO used for correlation
  int     ulim;                 //last  MO used for correlation
  
  ifstream datf(bcsfile);
  datf.read((char*) &nroao, sizeof(int));
  datf.read((char*) &nroe,  sizeof(int));
  datf.read((char*) &llim,  sizeof(int));
  datf.read((char*) &ulim,  sizeof(int));

  outf << "CIS data read: \n";
  int nrof = ulim - llim +1;
  int omo  = nroe/2 - llim;
  int umo  = ulim - nroe/2+1;
  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  check_range(nroao, llim, &outf);
  check_range(nroao, ulim, &outf);
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  int cis_size = omo*umo+1;
  long long int Cis_size =  cis_size;
  outf << "nr of CSF               : " << cis_size << "\n";
  outf << "Allocating Memory\n";

  double* cisvecs     = new double[(long long int) cis_size* (long long int) cis_size];
  double* cis_vals    = new double[cis_size];    //CIS (or CIS(D) eigenvalues)
  double* dumv        = new double[cis_size];    //temporary space

  Complex* wav        = new Complex[cis_size];   //Time dependent wave function (eigen state basis)
  Complex* SEwav      = new Complex[cis_size];   //Time dependent wave function (CSF basis)

  outf << "Reading data:\n";    
  datf.read((char *) cis_vals, cis_size * sizeof(double));
  datf.read((char *) cisvecs, (long long int) cis_size * (long long int) cis_size * sizeof(double));
  outf << "First uncorrected excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
  outf << "Highest uncorrected excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";
  //READ D CORRECTIONS
  if(use_d == 1){
    for(int x = 0; x < 4*cis_size+4; x++)
      datf.read((char *) dumv, cis_size * sizeof(double));
    for(int x = 0; x < cis_size; x++)
      cis_vals[x] += dumv[x];
    outf << "First CORRECTED excitation energies are: " << cis_vals[1] << " " << cis_vals[2] << "  ....\n";
    outf << "Highest CORRECTED excitation energy  is: " <<  cis_vals[cis_size-1] << "\n";
  }
  outf << "-------------------------------------------------------------------------------\n";
  
  datf.close();
  
  //OPEN NEEDED FILE STREAMS
  sprintf(dumc,"%s.exp",argv[2]);
  ofstream expf(dumc);
  outf << "Writing expectation values\n    [t,E,T,mu_x,mu_y,mu_z,";
  for(int x = 0; x < nroc; x++)
    outf << "MC(" << pop_centers[x] << "),";
  for(int x = 0; x < nrob; x++)
    outf << "MBO(" << pop_bonds_atom1[x] << "-" << pop_bonds_atom2[x] <<"),";
  outf << "]\nto " << dumc << "\n";


  sprintf(dumc,"%s.mop",argv[2]);
  ofstream mopf(dumc);
  outf << "\nWriting MO-pops to " << dumc << "\n";

  sprintf(dumc,"%s.mat",argv[2]);
  ofstream matf(dumc);
  outf << "\nWriting density matrices  to " << dumc << "\n";
  
  double   curr_time;
  
  /////////////INITIAL ANALYSIS////////////////////
  ifstream wvf(tdwfile);
  //Read first wave function for inital analysis
  wvf.read((char *) &curr_time, sizeof(double));
  wvf.read((char *)  wav, sizeof(Complex)*cis_size);

  outf << "Detailed analysis of the inital state:\n";
  
  calc_dens(nroe, nroao, llim, ulim, SEwav, Pmat, MOs);
  if(cis_size < 10000)
    transformationCIS_SE(wav, SEwav, cisvecs, cis_size);
  else
    transformationCIS_SE(wav, SEwav, cisvecs, Cis_size);
  calc_dens(nroe, nroao, llim, ulim, SEwav, Pmat, MOs);
  
  double E = 0.;
  for(int x = 0; x < cis_size; x++)
    E+= norm(wav[x])*cis_vals[x];

  
  outf << "Inital Energy: " << E << "\n";
  outf << "Kinetic energy: " << calc_op_1el(nroao, Tmat, Pmat) << "\n";
  outf << "Dipole moment (x,y,z) : " 
       << -calc_op_1el(nroao, Dx,   Pmat) +mu_core[0] << " "
       << -calc_op_1el(nroao, Dy,   Pmat) +mu_core[1] << " "
       << -calc_op_1el(nroao, Dz,   Pmat) +mu_core[2] << "\n";
  //Calc Mulliken matrix
  mat_mat(nroao, Pmat, Smat, Pmat_Mul);
  
  outf << "\nMulliken charges:\n";
  for(int c = 0; c < nroc; c++){
    double MullC = 0.;
    for(int b = first_basis[pop_centers[c]]; b <= last_basis[pop_centers[c]]; b++)
      MullC +=  Pmat_Mul[b*nroao+b];
    outf << "Center " << pop_centers[c] <<  " : " << -MullC+charges[pop_centers[c]] << "\n";
  }
  
  outf << "\nBond orders (Mayer):\n";
  for(int b = 0; b < nrob; b++){
    double MayerB = 0.;
    for(int c1 =  first_basis[pop_bonds_atom1[b]]; c1 <= last_basis[pop_bonds_atom1[b]]; c1++){
      for(int c2 =  first_basis[pop_bonds_atom2[b]]; c2 <= last_basis[pop_bonds_atom2[b]]; c2++){
	MayerB += Pmat_Mul[c1*nroao+c2]*Pmat_Mul[c2*nroao+c1];
      }
    }
    outf << "Center " << pop_bonds_atom1[b] << "----Center " << pop_bonds_atom2[b] << " : "  << MayerB << "\n";
  }
  
  calc_MOpops(nroe, nroao, llim, ulim,  SEwav,  tdPOP_mo);
  outf << "\nMO_populations:\n";
  for(int x = 0; x < nroao; x++){
    outf << "MO:\t" << x << "\tEn\t" << MOens[x] << "\t" << tdPOP_mo[x] << "\n";
  }
  
  outf << "\nHF-Ground state populations: " <<  norm(SEwav[0]) << "\n";
  outf << "\nLeading CSF populations (> 0.01):\n";
  for(int x = 1; x < cis_size; x++){
    int i = (x-1)/umo+llim;
    int f = (x-1)%umo+omo+llim;
    if(norm(SEwav[x]) > 0.01)
      outf << i << "-->" << f << ": " << norm(SEwav[x]) << " " << SEwav[x] << "\n";
  }
  
  

  outf << "\n*******************************************************************************\n";
  outf.flush();

  //////////////LOOP OVER TIME STEPS///////////////////
  for(int x = 0; x < nrots; x++){
    if(x != 0){
      wvf.read((char *) &curr_time, sizeof(double));
      wvf.read((char *)  wav, sizeof(Complex)*cis_size);
    }
    if(x%main_mw == 0){
      outf << "Current analysis time: " << curr_time << "\n";
      
      //EXPECTATION VALUES
      if((x/main_mw)%mw_exp == 0){
	if(cis_size < 10000)
	  transformationCIS_SE(wav, SEwav, cisvecs, cis_size);
	else
	  transformationCIS_SE(wav, SEwav, cisvecs, Cis_size);
	calc_dens(nroe, nroao, llim, ulim, SEwav, Pmat, MOs);
	E = 0.;
	for(int x = 0; x < cis_size; x++)
	  E+= norm(wav[x])*cis_vals[x];
	expf << curr_time << " " << E << " "
	     <<  calc_op_1el(nroao, Tmat, Pmat) << " "
	     << -calc_op_1el(nroao, Dx,   Pmat) +mu_core[0] << " "
	     << -calc_op_1el(nroao, Dy,   Pmat) +mu_core[1] << " "
	     << -calc_op_1el(nroao, Dz,   Pmat) +mu_core[2] << " ";
	//Calc Mulliken matrix
	mat_mat(nroao, Pmat, Smat, Pmat_Mul);
	
	//Calc charges
	for(int c = 0; c < nroc; c++){
	  double MullC = 0.;
	  for(int b = first_basis[pop_centers[c]]; b <= last_basis[pop_centers[c]]; b++)
	    MullC +=  Pmat_Mul[b*nroao+b];
	  expf << -MullC+charges[pop_centers[c]] << " ";
	}
	//Calc bond orders
	for(int b = 0; b < nrob; b++){
	  double MayerB = 0.;
	  for(int c1 =  first_basis[pop_bonds_atom1[b]]; c1 <= last_basis[pop_bonds_atom1[b]]; c1++){
	    for(int c2 =  first_basis[pop_bonds_atom2[b]]; c2 <= last_basis[pop_bonds_atom2[b]]; c2++){
	      MayerB += Pmat_Mul[c1*nroao+c2]*Pmat_Mul[c2*nroao+c1];
	    }
	  }
	  expf << MayerB << " ";
	}
	expf << "\n"; expf.flush();
      }
      
      outf.flush();
    }
    //MO POPS
    if((x/main_mw)%mw_mopop == 0){
      transformationCIS_SE(wav, SEwav, cisvecs, cis_size);
      calc_MOpops(nroe, nroao, llim, ulim,  SEwav,  tdPOP_mo);
      mopf << curr_time << " ";
      for(int x = 0; x < nroao; x++)
	mopf << tdPOP_mo[x]  << " ";
      mopf << "\n";
      mopf.flush();
    }
    //DENSITY MATRIX
    
    if((x/main_mw)%mw_mat == 0){
      if((x/main_mw)%mw_exp != 0){
	transformationCIS_SE(wav, SEwav, cisvecs, cis_size);
	calc_dens(nroe, nroao, llim, ulim, SEwav, Pmat, MOs);
      }
      matf << "Current time " << curr_time << "\n";
      matf << "---------------------------------------------------------------------------------------------------\n";
      output_matrix(Pmat,5, nroao, &matf);    
    }
  }
}

//extern void   calc_MOpops(int nroe, int nroao, int llim, int ulim, Complex* SEwav, double pops);

void check_range(int tomuch, int current_number, ofstream *outf){
  if( current_number  < 0 ||  current_number >= tomuch){
    *outf << "Invalid input: " << current_number  << "\n";
    outf->flush(); exit(1);
  }
}
