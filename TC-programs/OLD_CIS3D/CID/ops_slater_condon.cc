/********************************************************************************
 * Functions to calculate matrix elements between diffenrent Slater             *
 * determinants according to the Slater-Condon-Rules                            *
 *                                                                              *
 * Determinants (Configurations) are stored as integer fields, which            *
 * describe the Hartree product used to from the determinant.                   *
 * The index, i, stands for the position/spin  variable of the i th electron    *
 * the integer det[i] stands for the det[i] th spin orbital.                    *
 * Even numbers det[i] stand for the (det[i]/2)th spatial orbital times         *
 * an alpha spin function, det[i]+1 for the corresponding beta spin orbital.    *
 *  e.g.                                                                        *
 *                                                                              *
 * det[] = {0,1,2,7,4,5}                                                        *
 *                                                                              *
 * means the determinant is formed from permutations of the Hartree             * 
 * product                                                                      *
 *                                                                              *
 * chi_0(x_0)*chi_1(x_1)*chi_2(x_2)*chi_7(x_3)*chi_4(x_4)*chi_5(x_5)            *
 *                                                                              *
 * which is in terms of spatial orbitals psi                                    *
 *                                                                              *
 * psi_0(r_0)*alpha(0) * psi_0(r_1)*beta(1)  * psi_1(r_2)*alpha(1) *            *
 * psi_3(r_3)*beta(1)  * psi_2(r_4)*alpha(4) * psi_2(r_5)*beta(5)               *
 *                                                                              *
 *                                                                              *
 *                                                                              *
 * Slow rountines no assumption made, only that spin orbitals are               *
 * orthonormal.                                                                 *
 *                                                    T. Klamroth 2005          *
 ********************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

using namespace std;


//FUNCTIONS
int check_nr_diff(int nroe, int* det1, int* det2);
int bring_max_co(int nroe, int* det1, int* det2);
double calc_op1el_det(int nroe, int nroao, int* det1, int* det2, double* opmat);

//DEBUG
void print_dets_bra_ket(int nroe, int* det1, int* det2){
    clog << "<";
    for(int x = 0; x < nroe-1; x++)
      clog << det1[x] << ",";
    clog << det1[nroe-1] << "|O";

    clog << "|";
    for(int x = 0; x < nroe-1; x++)
      clog << det2[x] << ",";
    clog << det2[nroe-1] << ">\n";   
}

//DEBUG
void print_dets_ket(int nroe, int* det2){
    clog << "|";
    for(int x = 0; x < nroe-1; x++)
      clog << det2[x] << ",";
    clog << det2[nroe-1] << ">\n";   
}


/********************************************************************************
 *                                                                              *
 * Checks in how many spin orbitals two determinants are different              *
 *                                                                              *
 ********************************************************************************/

int check_nr_diff(int nroe, int* det1, int* det2){
  int nr_diff = 0;
  
  for(int x = 0; x < nroe; x++){
    int is_diff = 1;
    for(int y = 0; y < nroe; y++){
      is_diff *= (det1[x]+1)-(det2[y]+1);   //if ever zero, is_diff will be zero
    }
    if(is_diff != 0)
      nr_diff += 1;
  }
  return(nr_diff);
}

/********************************************************************************
 *                                                                              *
 * Perform one permutation between a and b                                      *
 *                                                                              *
 ********************************************************************************/

inline void swap_inds(int a, int b, int *det){
  int tmp = det[a];
  det[a]  = det[b];
  det[b]  = tmp;
}


/********************************************************************************
 *                                                                              *
 * Bring det2  in maximal conicidence with respect to det1                      *
 * Put the all the differences to the front.                                    *  
 *                                                                              *
 ********************************************************************************/

int bring_max_co(int nroe, int* det1, int* det2){
  int nr_perm = 0;
  //bring in maximal conincidence 
  for(int x = 0; x < nroe; x++){
    if(det1[x] != det2[x]){
      for(int y = 0; y < nroe; y++){
	if(det1[x] == det2[y]){ 
	  swap_inds(x,y, det2);
	  nr_perm++;
	  y = nroe;
	}
      }
    }
  }
  //put differences to the front
  int diff_count = 0;
  for(int x = 0; x < nroe; x++){
    if(det1[x] != det2[x]){
      swap_inds(x,diff_count, det1);
      swap_inds(x,diff_count, det2);
      diff_count++;
    }
  }
  return(nr_perm);
}

/********************************************************************************
 *                                                                              *
 * Calc one electron operator between 2 determinants                            *
 * opmat contains matrix elements in !!spatial!! MO basis                       *  
 *                                                                              *
 ********************************************************************************/

double calc_op1el_det(int nroe, int nroao, int* det1, int* det2, double* opmat){
  int nr_diff = check_nr_diff(nroe,  det1,  det2);
  
  double opval = 0.;
  
  if( nr_diff > 1){

    return(0.);   //Zero if different in more than 1 SO

  }
  
  //bring to maximal coincidence, save prefac
  int    nr_perm = bring_max_co(nroe, det1, det2);
  double fac     = pow(-1.,nr_perm);

  if( nr_diff == 1){  

    //Different in one SO 
    if((det1[0]+det2[0])%2 != 0) 
      return(0);  //if different spins everthing is zero
    //other wise one integral left
    opval =  opmat[(int)(det1[0]/2)*nroao+(int)(det2[0]/2)];

  }else{                         
    
    // all SO identical
    for(int x = 0; x < nroe; x++)
      opval += opmat[(int)(det1[x]/2)*nroao+(int)(det2[x]/2)];
    
  }
  return(fac*opval);   //give result back
  
}

/*******************************************************************************
 *                                                                              *
 * SO_calc_ijIIkl                                                               *
 *                                                                              *
 * Calculates <ij||kl> = <ij|kl>-<ik|jl>  Phyiscist NOTATION !!!!!!!!           *
 *******************************************************************************/

double SO_calc_ijIIkl(int i, int j, int k, int l, int nroao, double* prec_ints){
  static int off_i = nroao*nroao*nroao;
  static int off_j = nroao*nroao;
  static int off_k = nroao;

  int I = i/2;
  int J = j/2;
  int K = k/2;
  int L = l/2;
  
  double inte = 0.;
  

  //Calc <ij|kl> [ => (ik|jl)]
  // check spin
  if( (i+k)%2 == 0 && (j+l)%2 == 0)
    inte += +prec_ints[I*off_i + K*off_j + J*off_k + L];
    
  //Calc <ij|lk> [ => (il|kj)]
  // check spin
  if( (i+l)%2 == 0 && (k+j)%2 == 0) 
    inte += -prec_ints[I*off_i + L*off_j + K*off_k + J];

  return(inte);
}


/********************************************************************************
 *                                                                              *
 * Calc two electron operator between 2 determinants                            *
 * opmat contains matrix elements in !!spatial!! MO basis                       *  
 *                                                                              *
 ********************************************************************************/

double calc_op2el_det(int nroe, int nroao, int* det1, int* det2, double* prec_ints){
  int nr_diff = check_nr_diff(nroe,  det1,  det2);
  
  double opval = 0.;
  
  if( nr_diff > 2){

    return(0.);   //Zero if different in more than 2 SO
  }
  
  //bring to maximal coincidence, save prefac
  int    nr_perm = bring_max_co(nroe, det1, det2);
  double fac     = pow(-1.,nr_perm);


  if( nr_diff == 2){              
    //Different in two SO 
    opval += SO_calc_ijIIkl(det1[0], det1[1], det2[0], det2[1],  nroao, prec_ints);
  }
  
  if(nr_diff == 1){
    //Different in one SO
    for(int x = 1; x < nroe; x++)
      opval += SO_calc_ijIIkl(det1[0], det1[x], det2[0], det2[x], nroao,  prec_ints);
  }

  if(nr_diff == 0){
    //Different in 0 SO
    for(int x = 0; x < nroe; x++){
      for(int y = 0; y < nroe; y++)
	opval += SO_calc_ijIIkl(det1[x], det1[y], det2[x], det2[y], nroao, prec_ints);
    }
    opval *= 0.5;
  }
  return(fac*opval);   //give result back
  
}

/********************************************************************************
 *                                                                              *
 * Calc one electron operator between 2 linear combinations of determinants(CSF)*
 * opmat contains matrix elements in !!spatial!! MO basis                       *  
 *                                                                              *
 ********************************************************************************/

double calc_op1el_csf(int nroe, int nroao, 
		      int* bra_dets, double* bra_facs, int bra_nr,
		      int* ket_dets, double* ket_facs, int ket_nr,
		      int* tmp_det_bra, int* tmp_det_ket, 
		      double* opmat){

  double opval = 0.;
  size_t meml  = sizeof(int)*nroe;


  //loops over all simple determinant matrix elements 
  for(int x = 0; x < bra_nr; x++){
    for(int y = 0; y < ket_nr; y++){
      //Copy original dets, because they might be permuted for max coincidence.
      memcpy(tmp_det_bra, &(bra_dets[x*nroe]), meml);
      memcpy(tmp_det_ket, &(ket_dets[y*nroe]), meml);
      
      opval +=  bra_facs[x]*ket_facs[y]*
	calc_op1el_det(nroe, nroao, tmp_det_bra,  tmp_det_ket, opmat);
    }
  }
  return(opval);   
}

/********************************************************************************
 *                                                                              *
 * Calc one electron operator between 2 linear combinations of determinants(CSF)*
 * opmat contains matrix elements in !!spatial!! MO basis                       *  
 *                                                                              *
 ********************************************************************************/

double calc_op2el_csf(int nroe, int nroao, 
		      int* bra_dets, double* bra_facs, int bra_nr,
		      int* ket_dets, double* ket_facs, int ket_nr,
		      int* tmp_det_bra, int* tmp_det_ket, 
		      double* prec_ints){

  double opval = 0.;
  size_t meml  = sizeof(int)*nroe;


  //loops over all simple determinant matrix elements (Upper triangular)
  for(int x = 0; x < bra_nr; x++){
    for(int y = 0; y < ket_nr; y++){
      //Copy original dets, because they might be permuted for max coincidence.
      memcpy(tmp_det_bra, &(bra_dets[x*nroe]), meml);
      memcpy(tmp_det_ket, &(ket_dets[y*nroe]), meml);
      
      opval +=  bra_facs[x]*ket_facs[y]*
	calc_op2el_det(nroe, nroao, tmp_det_bra, tmp_det_ket, prec_ints);
    }
  }
  return(opval);   
}

