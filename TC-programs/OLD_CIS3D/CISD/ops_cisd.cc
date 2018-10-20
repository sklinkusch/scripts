/********************************************************************************
 * Functions for CIS example program with determinant expressions               *
 *                                                                              *
 * file: ops_cis.cc                                                             *
 *                                                                              *
 * contains the 4 index transformation                                          * 
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;


//Functions
double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
		 long long int nrofint, long long int* sortcount,  double* intval,
		 unsigned short* intnums, ofstream* outf);

		      
/*******************************************************************************
 *                                                                             *
 * calc_mo2int_op                                                              *
 *                                                                             *
 *                                                                             *
 * optimized     4 index transformation                                        *
 * should be called with i being the slowest varying index                     *
 * (4 index transformation)  (chemist notation)                                *
 *******************************************************************************/


//TRANSFORM THE FIRST INDEX AND BUILD ALL PERMUTATIONS OF THE REMAINING INDICES => AVOIDS K^4 MEMORY STORAGE
void trans_I(int i, double* i____Trans, int nroao, double* MOs,
	     long long int nrofint, long long int* sortcount,  double* intval,
	     unsigned short* intnums){

  //STATIC VARIABLES
  static long long int sM3 = (long long int) nroao *(long long int) nroao *(long long int) nroao;
  static long long int fac_d = nroao * nroao;   //!!!!MAKE b the fast running index => will be transfromed next
  static long long int fac_c = nroao;

  //erase old I;
  for(long long int x = 0; x < sM3; x++) i____Trans[x] = 0.;
  
  int off_i = nroao*i;


  ///!!!!!!!!!!!DO TRANSFORMATION and PERMUTATIONS!!!!!!!!!!!!!!!
  
  //PERM_1
  for(long long int x = 0; x < sortcount[0]; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!) 
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
  }
  
  //PERM_15
  for(long long int x = sortcount[0]; x < sortcount[1]; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!) 
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
  }

  //PERM_1234
  for(long long int x = sortcount[1]; x < sortcount[2]; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!) 
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc) 
    i____Trans[intnums[x*4+0]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
    i____Trans[intnums[x*4+0]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
  }

  
  //PERM_1256
  for(long long int x = sortcount[2]; x < sortcount[3]; x++) {
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!) 
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc) 
    i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
    i____Trans[intnums[x*4+2]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(6)    (dc|ab)
  }

 
  //PERM_ALL
  for(long long int x = sortcount[3]; x < nrofint; x++){
    //                     b=1            c=2                  d=3                                a=0   (Transfromation over a!!) 
    i____Trans[intnums[x*4+1]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(1)    (ab|cd)
    i____Trans[intnums[x*4+1]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+0]]*intval[x]; //(2)    (ab|dc) 
    i____Trans[intnums[x*4+0]+intnums[x*4+2]*fac_c+intnums[x*4+3]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(3)    (ba|cd)
    i____Trans[intnums[x*4+0]+intnums[x*4+3]*fac_c+intnums[x*4+2]*fac_d] += MOs[off_i+intnums[x*4+1]]*intval[x]; //(4)    (ba|dc)
    i____Trans[intnums[x*4+3]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(5)    (cd|ab)
    i____Trans[intnums[x*4+2]+intnums[x*4+0]*fac_c+intnums[x*4+1]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(6)    (dc|ab)
    i____Trans[intnums[x*4+3]+intnums[x*4+1]*fac_c+intnums[x*4+0]*fac_d] += MOs[off_i+intnums[x*4+2]]*intval[x]; //(7)    (cd|ba)
    i____Trans[intnums[x*4+2]+intnums[x*4+1]*fac_c+intnums[x*4+0]*fac_d] += MOs[off_i+intnums[x*4+3]]*intval[x]; //(8)    (dc|ba)
  }
}



double mo2int_op(int i, int j, int k, int l, int nroao, double* MOs,
		 long long int nrofint, long long int* sortcount,  double* intval,
		 unsigned short* intnums, ofstream* outf){
  
  static double* dmem;        //pointer for memory allocation
  static double* i____Trans;  //transformed ints for only i       (nroao^3)
  static double* ij___Trans;  //transformed ints for only ij      (nroao^2) 
  static double* ijk__Trans;  //transfodmed ints for only ijk     (nroao  )
  
  static int alt_i = -1;
  static int alt_j = -1;
  static int alt_k = -1;
 
  static long long int stati = 0, statj=0, statk = 0;


  static long long int fac_d = nroao * nroao;   //!!!!MAKE j the fast running index => will be transfromed next
  static long long int fac_c = nroao;
  
  
  if(dmem == NULL){
    *outf << "...............................................................................\n";
    *outf << "Memory allocation for improved 4 index transformation\n";
    long long int M = nroao;
    long long int dM = M*M*M+M*M+M;
    *outf << "Need " << dM*8 << " bytes of Memory for one index  at a time transformation\n";
    outf->flush();
    dmem = new double[dM];
    long long int inc = 0;
    i____Trans = &(dmem[inc]); inc += M*M*M;
    ij___Trans = &(dmem[inc]); inc += M*M;
    ijk__Trans = &(dmem[inc]); inc += M;    
    *outf << "Done \n";
    outf->flush();
    *outf << "...............................................................................\n";
  }

  //Print statistcs 
  if(i == -1){
    *outf << "\n...............................................................................\n";
    *outf << "mo2int_op transformation statistics: \n";
    *outf << "\nmo2int_op: # i trans so far " << stati << "\n";
    *outf <<   "mo2int_op: # j trans so far " << statj << "\n";
    *outf <<   "mo2int_op: # k trans so far " << statk << "\n"; 
    *outf << "...............................................................................\n";
    outf->flush();
    return(0);
  }



  if(i != alt_i){
    //Make transformation for i
    trans_I( i,  i____Trans,  nroao,  MOs, nrofint, sortcount,   intval, intnums);
    alt_i = i;
    //Statistics
    stati++;
  }
  
  if(j != alt_j){
    //Make transformation    j
    for(int d = 0; d < nroao; d++){
      for(int c = 0; c < nroao; c++){
	ij___Trans[d*nroao+c] = 0.;
	long long int ioff_b = fac_c*c +  fac_d*d;
	int off_j = nroao*j;
	for(int b = 0; b < nroao; b++)
	  ij___Trans[d*nroao+c] += i____Trans[ioff_b+b]*MOs[off_j+b];
      }
    }
    alt_j = j;
    //Statistics
    statj++;
  }
  
  if(k != alt_k){
    //Make transformation for k
    for(int d = 0; d < nroao; d++){
      int off_k = k*nroao;
      ijk__Trans[d] = 0.;
      for(int c = 0; c < nroao; c++)
	ijk__Trans[d] += ij___Trans[d*nroao+c]*MOs[off_k+c];
    }
    alt_k = k;
    //Statistics
    statk++;
  }

  
  //Make transfromation for l
  double inte = 0.;
  int off_l = nroao*l;
  for(int d = 0; d < nroao; d++)
    inte += ijk__Trans[d]*MOs[off_l+d];
  
  return(inte);
}




