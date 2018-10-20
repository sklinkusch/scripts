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
void cp_mus(int cisd_size, int nroao, int llim, int ulim,  int nroe, double* cisdmat, double* cisdvecs, double* D,
            double mu_core, double* MOs, double* Mat_m11, double* tmpvecs, double* cisdtmpmat,
            int* bra_dets, double* bra_facs, int bra_nr,
            int* ket_dets, double* ket_facs, int ket_nr,
            int* tmp_det_bra, int* tmp_det_ket,
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
extern void   calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void   calc_mu_core(int nroa, double* coord, double* charges, double* point, 
                         double* mu_core);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern void diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void pmv(double* mat, double* vi, double* vo, int nroao);
extern double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
                        long long int nrofint, long long int* sortcount,  double* intval,
                        unsigned short* intnums, ofstream* outf);

extern void calc_det(int i, int nroe, int llim, int ulim, int* dets, double* facs, int* nr);
extern double calc_op1el_csf(int nroe, int nroao, 
                             int* bra_dets, double* bra_facs, int bra_nr,
                             int* ket_dets, double* ket_facs, int ket_nr,
                             int* tmp_det_bra, int* tmp_det_ket, 
                             double* opmat);
extern double calc_op2el_csf(int nroe, int nroao, 
                             int* bra_dets, double* bra_facs, int bra_nr,
                             int* ket_dets, double* ket_facs, int ket_nr,
                             int* tmp_det_bra, int* tmp_det_ket, 
                             double* prec_ints);

//DEBUG
extern void print_dets_ket(int nroe, int* det2);

int main(int argc, char* argv[]){

  //INPUT VARIABLES
   int    nroe;            //Nr of electrons  (if negative, read in center of mass)
   int    llim;            //first MO used for correlation
   int    ulim;            //last  MO used for correlation
   char sysfile[256];      //binary system file
   char wavfile[256];      //HF-Wavefunction file
   int calc_COM = 0;       // 0 -> calc COM otherwise read in

   if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }

  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  

  outf << "CID \n";
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

  ist >> llim >> ulim >> wavfile;
  outf << "System data \nNr of electrons: " << nroe << "\nsys-file: " << sysfile << "\n\n";
  outf << "CID data\nLower limit of MOs used for correlation: " << llim << "\n";
  outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
  outf << "Reading HF-wave function from  " << wavfile << "\n";

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
  double*        Hmat;          //one electron Hamiltionian         nroao*nroao
  double*        Tmat;          //Kinetic energy operator           nroao*nroao
  

  double*        MOs;           //MO coeffs                         nroao*nroao
  double*        Dx;            //Dipole X                          nroao*nroao

  double*        Dy;            //Dipole Y                          nroao*nroao
  double*        Dz;            //Dipole Z                          nroao*nroao

  double*        MOens;         //MO Energies                       nroao
  double*        Mat_h11;       //Energie of 1p hamiltionian        nroao*nroa
  double*        Mat_m11;       //Dipole operator in MO spac        nroao*nroa


  //Temporary memory spaces
  double*  tmpmat1;                //                               nroao*nroao
  double*  tmpmat2;                //                               nroao*nroao
  double*  tmpvecs;                //                               nroao*nroao

  double*  tmpvals;                //                               nroao
 
  //MEMORY ALLOCATION for one electron atoms, mat&vecs
  int atom_ao_mem = 5*nroa+11*nroao*nroao+2*nroao;
  outf << "Need " << atom_ao_mem*sizeof(double) << " bytes for atomic + one electron data\n";
  outf.flush();
  
  double* dumd  = new double[atom_ao_mem]; int inc = 0;

  coord = &(dumd[inc]);   inc += 3*nroa; 
  charges = &(dumd[inc]); inc += nroa; 
  mass  = &(dumd[inc]);   inc += nroa;
  Hmat = &(dumd[inc]);    inc += nroao*nroao; 
  Tmat = &(dumd[inc]);    inc += nroao*nroao;
  MOs = &(dumd[inc]);     inc += nroao*nroao; 
  Dx = &(dumd[inc]);      inc += nroao*nroao;
  Dy = &(dumd[inc]);      inc += nroao*nroao; 
  Dz = &(dumd[inc]);      inc += nroao*nroao;
  MOens = &(dumd[inc]);   inc += nroao; 
  Mat_h11 = &(dumd[inc]); inc += nroao*nroao; 
  Mat_m11 = &(dumd[inc]); inc += nroao*nroao;
  
  tmpmat1 = &(dumd[inc]); inc += nroao*nroao; 
  tmpmat2 = &(dumd[inc]); inc += nroao*nroao; 
  tmpvecs = &(dumd[inc]); inc += nroao*nroao;
  
  tmpvals = &(dumd[inc]); inc += nroao;

  //MEMORY ALLOCATION for two electron values 
  outf << "Need " << nrofint*(sizeof(double)+sizeof(unsigned short)*4) << " bytes for two electron data\n";
  outf.flush();
  
  double*         intval        = new double[nrofint];                     //two electron integrals
  unsigned short* intnums       = new unsigned short[nrofint*4];           //two electron indices
  long long int   sortcount[4];                                            //num of two electron Integrals in each perm. type

  
  //Read system data 
  read_sys(sysfile, coord, charges, mass, Hmat, Tmat, tmpmat1, Dx, Dy,  Dz, sortcount, intval, intnums);

    //System Properties for HF
  outf << "-------------------------------------------------------------------------------\n";
  outf << "System properties:\n";
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

  read_wav_HF(wavfile, nroao, MOens, MOs);
  outf << "HF-wave function  read from " << wavfile << "\n";
  
  outf << "-------------------------------------------------------------------------------\n";
  

  outf << "Starting CID calculation with Determinant procedure:\n";   
  sprintf(dumc,"%s.bcs",argv[2]);
  outf << "Writing binary data to " << dumc << "\n";
  ofstream datf(dumc);
  
  outf << "\n\n";
    
 //CHECK limits
  if(llim < 0)      { outf << "Invalid lower limit in MO-range!\n"; outf.flush(); exit(4);}
  if(ulim >= nroao) { outf << "Invalid upper limit in MO-range!\n"; outf.flush();exit(4);}
  if(llim >= nroe/2){ outf << "No occupied orbitals in MO-range!\n"; outf.flush(); exit(4);}

  
  int omo      = nroe/2 - llim;
  int umo      = ulim - nroe/2+1;
  int nrof = ulim - llim +1;
  int cid_size;


  // Loops to determine cisd_size
  int sum2 = 0;
  for (int occ = 1; occ <= omo-1; occ++){
    sum2 += occ; 
  }

  int sum1 = 0;
  for (int virt = 1; virt <= umo-1; virt++){
    sum1 += virt; 
  }

  cid_size = 1 + omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;
  long long int Cid_size = cid_size;


  outf << "# MOs for correlation   : " << nrof << "\n";
  outf << "limits (l,u)            : " << llim << " , " << ulim << "\n";
  outf << "occupied                : " << omo << "\n";
  outf << "virtual                 : " << umo  << "\n";
  outf << "Nr of CSF               : " << cid_size << "\n";

  //Write data to data-file
  datf.write((char *) &nroao, sizeof(int));
  datf.write((char *) &nroe, sizeof(int));
  datf.write((char *) &llim, sizeof(int));
  datf.write((char *) &ulim, sizeof(int));
  


  outf << "-------------------------------------------------------------------------------\n";
  outf << "Allocating Memory for CID:\n";
  
  //CIS Variables
  outf << "Need " << (3*Cid_size*Cid_size+Cid_size)*sizeof(double) << " bytes of Memomry for CISD calculation\n";
  outf.flush();

  double*  cidmat    = new double[Cid_size*Cid_size];     //Matrix for Operator representations in CID-Basis
  double*  cidtmpmat = new double[Cid_size*Cid_size];     //Temp space
  double*  cidvals   = new double[cid_size];               //Space for excpetation values
  double*  cidvecs   = new double[Cid_size*Cid_size];     //Space for eigen vectors



  //PRECAL ONE ELECTRON HAMILTIONIAN IN AO BASIS
  outf << "Calculating one particle energies\n";
  //calc h(1) C_x
  for(int x = 0; x < nroao; x++){
    pmv(Hmat, &(MOs[x*nroao]), &(tmpvecs[x*nroao]), nroao);
  }
  //calc c_y h(1) c_x
  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++){
      Mat_h11[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++) Mat_h11[x*nroao+y] += MOs[x*nroao+z]*tmpvecs[y*nroao+z];
    }
  }

  //Do complete for index transformation 
  //VERY INEFFECTIVE CAN BE IMPROVED !!!!
  long long size_2el_ints = (long long int) nroao * (long long int) nroao * (long long int) nroao * (long long int) nroao;
  outf << "Allocationg " << size_2el_ints/131072 << " MB for 2 el integrals\n";
  double* prec_ints = new double[size_2el_ints];
  long long int prec_count = 0;
  for(int i = 0; i < nroao; i++){
    for(int j = 0; j < nroao; j++){
      for(int k = 0; k < nroao; k++){
        for(int l = 0; l < nroao;l++){
          prec_ints[prec_count++] = mo2int_op(i, j, k, l,
                                              nroao, MOs, nrofint, sortcount, intval, intnums, &outf);
        }
      }
    }
    outf << i << "\t";   
    if((i+1)%10==0) outf << "\n";
    outf.flush();
  }
  outf << "\nInts claculated: " << prec_count << "\n";
  outf << "\nDone\n";
  outf.flush();
  //Print out stat 
  mo2int_op(-1, -1, -1, -1,
            nroao, MOs, nrofint, sortcount, intval, intnums, &outf); 


  //SPACE FOR DETERMINANTS
  int max_nr_of_dets = 6;        //maximal 6 determinats in singlet CSFs

  int     bra_nr;
  double* bra_facs = new double[max_nr_of_dets];
  int*    bra_dets = new int[max_nr_of_dets*nroe];

  int     ket_nr;
  double* ket_facs = new double[max_nr_of_dets];
  int*    ket_dets = new int[max_nr_of_dets*nroe];

  int* tmp_det_bra = new int[nroe];
  int* tmp_det_ket = new int[nroe];
  
  //BUILD Hamiltionian in CID-Basis
  for(int x = 0; x < Cid_size*Cid_size; x++)
    cidmat[x] = 0.;

  for(int x = 0; x < Cid_size; x++)
    cidmat[x*Cid_size+x] =ion_rep;
  
  for (int i = 0; i < cid_size; i++){
    calc_det(i, nroe, llim, ulim, bra_dets, bra_facs, &bra_nr);
//     clog << i << ":\n";
//     for(int x = 0; x < bra_nr; x++){
//       clog << bra_facs[x] << "*"; print_dets_ket(nroe, &(bra_dets[x*nroe]));
//     }
    for(int j = i; j < cid_size; j++){
      calc_det(j, nroe, llim, ulim, ket_dets, ket_facs, &ket_nr);
      cidmat[i*Cid_size+j] +=  (calc_op1el_csf(nroe, nroao,  
						 bra_dets, bra_facs, bra_nr,
						 ket_dets, ket_facs, ket_nr,
						 tmp_det_bra, tmp_det_ket,
						 Mat_h11)+
				  calc_op2el_csf(nroe, nroao, 
						 bra_dets, bra_facs, bra_nr,
						 ket_dets, ket_facs, ket_nr,
						 tmp_det_bra, tmp_det_ket,
						 prec_ints));
      
    }
  }
  
  for(int x = 0; x < Cid_size; x++){
    for(int y = x; y < Cid_size; y++){
      cidmat[y*Cid_size+x] = cidmat[x*Cid_size+y];
    }
  }
  
  outf.precision(12);
  //  outf << cidmat[0] << "\n";
  

  diag_mat(cid_size, cidmat, cidvals, cidvecs); 
  
  //WRITE DATA TO FILE
  datf.write((char *) cidvals, Cid_size*(long long int) sizeof(double));
  datf.write((char *) cidvecs, Cid_size*Cid_size*(long long int) sizeof(double));
  datf.flush();

  outf << "Ground state energy: " << cidvals[0] << "\n";

  
  int nprint =  8;
  if(nprint > cid_size) nprint = cid_size;
  
  outf << "First " << nprint-1 << " exc. states\n [energies:  (Hartree) (eV),  (LC  x->y)  coeff :\n";
  outf << "...............................................................................\n";

  for(int x = 1; x < nprint; x++){
    double max_coeff = 0.;
    int mi = 0, mf = 0;
    for(int y = 0; y < cid_size; y++){
      if( fabs(cidvecs[x*cid_size+y]) > fabs(max_coeff)){
        max_coeff = cidvecs[x*cid_size+y];
        mi = (y-1)/umo+llim;
        mf = (y-1)%umo+omo+llim;
      }
    }
    sprintf(dumc,"%.8f\t%.8f,\t %.3i->%.3i\t%+.8f\n",(cidvals[x]-cidvals[0]),(cidvals[x]-cidvals[0])*27.211383,mi,mf,max_coeff/sqrt(2.));
    outf << x << "          "<< dumc;
  }
  
  outf << "-------------------------------------------------------------------------------\n";
 


  outf << "Building dipole transition marticies!\n";

  outf << "...............................................................................\n";
  outf << "X\n";
  outf.flush();
  cp_mus(cid_size,  nroao, llim, ulim, nroe,  cidmat,  cidvecs, Dx,
         mu_core[0], MOs, Mat_m11, tmpmat1, cidtmpmat,
         bra_dets, bra_facs, bra_nr,
         ket_dets, ket_facs, ket_nr,
         tmp_det_bra, tmp_det_ket,
         &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cid_size, cidmat, cidvals, cidtmpmat);
  outf << "Values ranging from " << cidvals[0] << " to " << cidvals[cid_size-1] << "\n";
  datf.write((char *) cidvals, cid_size*sizeof(double));
  datf.write((char *) cidtmpmat, Cid_size*Cid_size*(long long int) sizeof(double));
  datf.flush();

  outf << "...............................................................................\n";
  outf << "Y\n";
  outf.flush();
  cp_mus(cid_size,  nroao, llim, ulim, nroe,  cidmat,  cidvecs, Dy,
         mu_core[1], MOs, Mat_m11, tmpmat1, cidtmpmat,
         bra_dets, bra_facs, bra_nr,
         ket_dets, ket_facs, ket_nr,
         tmp_det_bra, tmp_det_ket,
         &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cid_size, cidmat, cidvals, cidtmpmat);
  outf << "Values ranging from " << cidvals[0] << " to " << cidvals[cid_size-1] << "\n";
  datf.write((char *) cidvals, cid_size*sizeof(double));
  datf.write((char *) cidtmpmat, Cid_size*Cid_size*(long long int) sizeof(double));
  datf.flush();


  outf << "...............................................................................\n";
  outf << "Z\n";
  outf.flush();
  cp_mus(cid_size,  nroao, llim, ulim, nroe,  cidmat,  cidvecs, Dz,
         mu_core[2], MOs, Mat_m11, tmpmat1, cidtmpmat,
         bra_dets, bra_facs, bra_nr,
         ket_dets, ket_facs, ket_nr,
         tmp_det_bra, tmp_det_ket,
         &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cid_size, cidmat, cidvals, cidtmpmat);
  outf << "Values ranging from " << cidvals[0] << " to " << cidvals[cid_size-1] << "\n";
  datf.write((char *) cidvals, cid_size*sizeof(double));
  datf.write((char *) cidtmpmat, Cid_size*Cid_size*(long long int) sizeof(double));
  datf.flush();

  outf << "\n" << "Alle Hoellenhunde!!" << "\n";
  //outf << "\n" << "Mille sabords!!" << "\n";

}

void cp_mus(int cid_size, int nroao, int llim, int ulim,  int nroe, double* cidmat, double* cidvecs, double* D,
            double mu_core, double* MOs, double* Mat_m11, double* tmpvecs, double* cidtmpmat,
            int* bra_dets, double* bra_facs, int bra_nr,
            int* ket_dets, double* ket_facs, int ket_nr,
            int* tmp_det_bra, int* tmp_det_ket,
            ofstream *outf, int nprint){
  char dumc[512];
  long long int Cid_size = cid_size;

  *outf << "Claculating mu op in MO-basis\n";
  for(int x = 0; x < nroao; x++){
    pmv(D, &(MOs[x*nroao]), &(tmpvecs[x*nroao]), nroao);
  }
  //calc c_y h(1) c_x
  for(int x = 0; x < nroao; x++){
    for(int y = 0; y < nroao; y++){
      Mat_m11[x*nroao+y] = 0.;
      for(int z = 0; z < nroao; z++) Mat_m11[x*nroao+y] -= MOs[x*nroao+z]*tmpvecs[y*nroao+z];
    }
  }
  
  outf->flush();

  //CIS loop in determinant space
  *outf << "Calculating upper triangle of CID-matrix\n";
  for(int x = 0; x < cid_size*cid_size; x++)
    cidmat[x] = 0.;

  //ADD MU CORE
  for(int x = 0; x < cid_size; x++)
    cidmat[x*cid_size+x] += mu_core;
  
  for (int i = 0; i < cid_size; i++){
    calc_det(i, nroe, llim, ulim, bra_dets, bra_facs, &bra_nr);
    for(int j = i; j < cid_size; j++){
      calc_det(j, nroe, llim, ulim, ket_dets, ket_facs, &ket_nr);      
      cidmat[i*cid_size+j] +=  calc_op1el_csf(nroe, nroao,
						bra_dets, bra_facs, bra_nr,
						ket_dets, ket_facs, ket_nr,
						tmp_det_bra, tmp_det_ket,
						Mat_m11);
    }
  }
    
  for(int x = 0; x < cid_size; x++){
    for(int y = x; y < cid_size; y++){
      cidmat[y*cid_size+x] = cidmat[x*cid_size+y];
    }
  }
  
  *outf << "Transforming Dipole matrix\n";
  outf -> flush();
  
  if(cid_size < 10000)
    trans_mat(cid_size, cidmat, cidvecs, cidtmpmat, 1); //Forward transfornation
  else
    trans_mat(Cid_size, cidmat, cidvecs, cidtmpmat, 1); //Forward transfornation
  
  
  *outf << "Transition Moments\n";
  for(int x = 0; x < nprint; x++) *outf << "        " << x;
  *outf << "\n";

  for(long long int x = 0; x < nprint ; x++ ){
    char tmpc[12];
    sprintf(dumc,"%i   ",(int) x);
    for(long long int y = 0; y < nprint; y++) {
      sprintf(tmpc,"%+.4f  ", cidmat[x*cid_size+y]);
      strcat(dumc,tmpc);
    }
    *outf << dumc << "\n";
  }
  *outf << "\n\n";
  outf->flush();
} 





