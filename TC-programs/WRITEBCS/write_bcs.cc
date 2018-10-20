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
void cp_mus(int cis_size, double* mumat, ofstream *outf, int nprint);
void read_gamess(double* data, int recnr, int nromo, int type, char* filename);


//Extern Functions
extern int  rem_com(char* filename, char* streamstring, int string_length);
extern void status(ofstream* outf);
extern void  diag_mat(int nroao, double* mat, double* vals, double* vecs);
extern void calc_mu_mat_cis(int cis_size, int nroao, int omo, int umo, int llim, int nroe, 
			    double *cismat, double* mumat, double mucore, 
			    double* MOs, double* Pmat, double* tmpvec_ao);
extern void trans_mat(long long int np, double* mat, double* trans, double* tmpmat, int dir);
extern void trans_mat(int np, double* mat, double* trans, double* tmpmat, int dir);
extern void symmortho_mat(int nroao, double *mat, double* tmat, double* dummat);
extern void pmv(double* mat, double* vi, double* vo, int nroao);
extern void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);
extern double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old, 
		       double* MOs, double damp);
extern void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
extern void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
		  double* mu_core);
extern double calc_ion_rep(int nroa, double* coord, double* charges);
extern void build_cismat(long long int Cis_size, double* cismat, double* cisvals, double* cisvecs);
extern void mat_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);         //for large matrices
extern void mat_mat_T(long long int np, double* mat_i1, double* mat_i2, double* mat_f);       //for large matrices
extern void mat_T_mat(long long int np, double* mat_i1, double* mat_i2, double* mat_f);
extern "C" void dafrd_(double* V, long long int* LEN, long long int* RECN, char* FNAME);

int main(int argc, char* argv[]){
  
  //INPUT VARIABLES
  char enfile[256];
  char vecfile[256];
  char dipfile[256];
  int nroao;      // Nr of atomic orbitals
  int nroe;      // Nr of electrons
  int nroa;
  int llim;
  int ulim;
  int calc_d=0;         //if 1 -> calc duplett correction PT

  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "WRITE_BCS [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

  ist >> enfile >> vecfile >> dipfile >> nroe >> nroao >> nroa >> llim >> ulim >> calc_d;

  if(nroe%2 != 0){
    cerr << "Nr of electrons not even !!!\n";
    exit(2);
  }
 
  outf << "File containing CIS energies and CIS(D) corrections: " << enfile << "\n";
  outf << "File containing CIS eigenvectors: " << vecfile << "\n";
  outf << "File containing dipole moments: " << dipfile << "\n";
  outf << "System data \nNr of electrons: " << nroe << "\n";
  outf << "CIS data\nLower limit of MOs used for correlation: " << llim << "\n";
  outf << "Upper limit of MOs used for correlation: " << ulim << "\n";
  outf << "Nr of basis functions: " << nroao   << "\n";
  outf << "Nr of atoms: " << nroa << "\n";
  outf.flush();

  sprintf(dumc,"%s.bcs",argv[2]);
  outf << "Writing binary data to " << dumc << "\n";
  ofstream datf(dumc);

  outf << "\n\n";
  outf.flush();

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
  outf.flush();

  outf << "nroao, nroe, llim, and ulim written to bcs-file\n";
  outf.flush();

  double*  cistmpmat = new double[Cis_size*Cis_size];     //Temp space
  double*  cisvals   = new double[Cis_size];              //Space for excpetation values
  double*  cisvecs   = new double[Cis_size*Cis_size];     //Space for eigen vectors
  double*  cisd_corr = new double[Cis_size];

// Read CIS Energies and CIS(D) corrections
  ifstream ensf;
  ensf.open(enfile);
  for(int x = 1; x < cis_size; x++)   ensf >> cisvals[x] >> cisd_corr[x];
  ensf.close();
  outf << "Eigenenergies and corrections: \n";
  for(int x = 0; x < cis_size; x++){
      outf << x << " " << cisvals[x] << " " << cisd_corr[x] << "\n";
      outf.flush();
  }
  
  // Read CI vectors
  outf << "Eigenvectors: \n";
  ifstream vecf;
  vecf.open(vecfile);
  int st, mi, mf;
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  if(x == 0 && y == 0){
	      cisvecs[x*cis_size+y] = 1.;
	      continue;
	  }else if(x == 0 &&  y != 0){
	      cisvecs[x*cis_size+y] = 0.;
	      continue;
	  }else if(x != 0 && y == 0){
	      cisvecs[x*cis_size+y] = 0.;
	      continue;
	  }else{
	      vecf >> st >> mi >> mf >> cisvecs[x*cis_size+y];
	  }
      }
  }
  vecf.close();
  for(int x = 0; x < cis_size; x++){
      outf << x << ": ";
      for(int y = 0; y < cis_size; y++){
	  outf << cisvecs[x*cis_size+y] << " ";
	  if(y%5 == 4) outf << "\n";
      }
      outf << "\n";
      outf.flush();
  }

  //WRITE DATA TO FILE
  datf.write((char *) &nroao, sizeof(int));
  datf.write((char *) &nroe, sizeof(int));
  datf.write((char *) &llim, sizeof(int));
  datf.write((char *) &ulim, sizeof(int));
  datf.write((char *) cisvals, Cis_size*(long long int) sizeof(double));
  datf.write((char *) cisvecs, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();
  int nprint = 8;

  double* unitmat = new double[Cis_size*Cis_size];
  mat_mat_T(Cis_size, cisvecs, cisvecs, unitmat); 
  double dtrace = 0.;
  double doffdiagonal = 0.;
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < cis_size; y++){
	  if (x == y){
	      dtrace += fabs(unitmat[x*cis_size+y] - 1.);
	  }else{
	      doffdiagonal += fabs(unitmat[x*cis_size+y]);
	  }
      }
  }
  outf << "Cisvecs Test Trace: " << dtrace << " " << dtrace/cis_size << "\n";
  outf << "Cisvecs Test Off-diagonal: " << doffdiagonal << " " << doffdiagonal/(cis_size*cis_size-cis_size) << "\n";
  outf.flush();

  double* dx = new double[Cis_size*Cis_size];
  double* dy = new double[Cis_size*Cis_size];
  double* dz = new double[Cis_size*Cis_size];
  ifstream dipf;
  dipf.open(dipfile);
  for(int x = 0; x < cis_size; x++){
      for(int y = x; y < cis_size; y++){
	  dipf >> dx[x*cis_size+y] >> dy[x*cis_size+y] >> dz[x*cis_size+y];
      }
  }
  dipf.close();
  for(int x = 0; x < cis_size; x++){
      for(int y = 0; y < x; y++){
	  dx[x*cis_size+y] = dx[y*cis_size+x];
	  dy[x*cis_size+y] = dy[y*cis_size+x];
	  dz[x*cis_size+y] = dz[y*cis_size+x];
      }
  }


  outf << "-------------------------------------------------------------------------------\n";

  outf << "Building dipole transition marticies!\n";

  outf << "...............................................................................\n";
  outf << "X\n";
  outf.flush();
  cp_mus(cis_size,  dx, &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, dx, cisvals, cistmpmat);
  outf << "Values ranging from " << cisvals[0] << " to " << cisvals[cis_size-1] << "\n";
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();

  outf << "...............................................................................\n";
  outf << "Y\n";
  outf.flush();
  cp_mus(cis_size,  dy, &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, dy, cisvals, cistmpmat);
  outf << "Values ranging from " << cisvals[0] << " to " << cisvals[cis_size-1] << "\n";
  datf.write((char *) cisvals, cis_size*sizeof(double));
  datf.write((char *) cistmpmat, Cis_size*Cis_size*(long long int) sizeof(double));
  datf.flush();


  outf << "...............................................................................\n";
  outf << "Z\n";
  outf.flush();
  cp_mus(cis_size,  dz, &outf, nprint);

  outf << "Calculating dipole eigenvectors \n";
  diag_mat(cis_size, dz, cisvals, cistmpmat);
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

void cp_mus(int cis_size, double* mumat, ofstream *outf, int nprint){
  char dumc[512];
  long long int Cis_size = cis_size;
  

  outf->flush();

  *outf << "Transition Moments\n";
  for(int x = 0; x < nprint; x++) *outf << "        " << x;
  *outf << "\n";
  
  for(long long int x = 0; x < nprint ; x++ ){
    char tmpc[12];
    sprintf(dumc,"%i   ",(int) x);
    for(long long int y = 0; y < nprint; y++) {
      sprintf(tmpc,"%+.4f  ", mumat[x*Cis_size+y]);
      strcat(dumc,tmpc);
    }
    *outf << dumc << "\n";
  }
  *outf << "\n\n";
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

