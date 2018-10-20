/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_so_corr.cc                                                         *
 *                                                                              *
 * operarotrs in spin orbital sapce  needed for cis(d) mp2 clacluations         * 
 *                                                                              *
 * NOTATION:                                                                    *
 *  The Ith spin orbital is (I integer >= 0)                                    *
 *                                                                              *
 *    if I%2 == 0                                                               * 
 *       the I/2 spatial orbtital times a alpha spin function                   *
 *                                                                              *
 *    if I%2 == 1                                                               *
 *       the I/2 spatial orbtital times a beta  spin function                   *
 *                                                                              *
 *                                                                              *
 *                                                                              *
 *  2 ELECTRON INTEGRALS IN PHYSICIST NOTATION IN THIS FILE !!!!!!!!!!!         *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/ 


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;


//Functions
inline int SO_to_MO(int SO);
double SO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);
double SO_calc_ijIIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints);
double SO_calc_ijIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints);
double SO_calc_a_ij_kl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints, double* MOens);
double SO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);
double SO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens);
double SO_calc_u_ij_ab(int i, int j, int a, int b,  int omo, int umo, int llim, double* cis_vec, double* prec_ints);
double SO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en, double* prec_ints, double* MOens);
double SO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens);
double SO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double cis_en, double* t1, double* t2);
//Extern Functions
extern double get_precalc_ints_sd(int i, int j, int k, int l,
				  int omo, int umo, int llim, double* prec_ints);

/*******************************************************************************
 *                                                                             *
 *  SO_to_MO                                                                   *
 *                                                                             *
 *                                                                             *
 * returns the corresponding spatiol orbital +1 !!!!!!!! corresponding to SO ! *
 * value is negative for beta spin                                             *
 *******************************************************************************/

inline int SO_to_MO(int SO){
  int mo_nr = SO/2+1;
  if(SO%2 != 0)  mo_nr *= -1;
  return(mo_nr);
}

/*******************************************************************************
 *                                                                              *
 * SO_calc_delta_ij_ab                                                          *
 *                                                                              *
 * 4 mo energie difference needed for MP2 and CIS(D)                            *
 *******************************************************************************/

double SO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens){
  return(-MOens[abs(SO_to_MO(i))-1]-MOens[abs(SO_to_MO(j))-1]
	 +MOens[abs(SO_to_MO(a))-1]+MOens[abs(SO_to_MO(b))-1]);
}
			   
/*******************************************************************************
 *                                                                              *
 * SO_calc_ijIIkl                                                               *
 *                                                                              *
 * Calculates <ij||kl> = <ij|kl>-<ik|jl>  Phyiscist NOTATION !!!!!!!!           *
 *******************************************************************************/

double SO_calc_ijIIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints){
  int I = SO_to_MO(i);
  int J = SO_to_MO(j);
  int K = SO_to_MO(k);
  int L = SO_to_MO(l);
  
  double inte = 0.;
  
  //Calc <ij|kl> [ => (ik|jl)]
  if(I*K >= 0 && J*L >= 0) inte += +get_precalc_ints_sd(abs(I)-1, abs(K)-1, abs(J)-1, abs(L)-1,
						       omo, umo, llim, prec_ints);

  //Calc <ij|lk> [ => (il|kj)]
  if(I*L >= 0 && K*J >= 0) inte += -get_precalc_ints_sd(abs(I)-1, abs(L)-1, abs(K)-1, abs(J)-1,
						       omo, umo, llim, prec_ints);

  return(inte);
}
						
/*******************************************************************************
 *                                                                              *
 * SO_calc_ijIkl                                                                *
 *                                                                              *
 * Calculates <ij|kl>                                                           *
 *******************************************************************************/

double   SO_calc_ijIkl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints){
  int I = SO_to_MO(i);
  int J = SO_to_MO(j);
  int K = SO_to_MO(k);
  int L = SO_to_MO(l);
  
  double inte = 0.;
  
  //calc <ij|kl> [ => (ik|jl)
  if(I*K>= 0 && J*L >= 0) inte += +get_precalc_ints_sd(abs(I)-1, abs(K)-1, abs(J)-1, abs(L)-1,
						       omo, umo, llim, prec_ints);

  return(inte);
}
						

/*******************************************************************************
 *                                                                              *
 * SO_calc_a_ij_kl                                                              *
 *                                                                              *
 * Calculates a_ij_kl = -1/DELTA_ij^kl * <kl||ij>                               *
 *                                                                              *
 * NOTE <ij||kl> used instead of <kl||ij> since all MOs are real and it is      *
 * faster this way !!!!!!!!!!!!!!!!!                                            *
 *******************************************************************************/

double SO_calc_a_ij_kl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints, double* MOens){
  
  return(-SO_calc_ijIIkl( i, j, k, l, omo, umo, llim, prec_ints)/SO_calc_delta_ij_ab( i,  j,  k,  l, MOens));
}

/*******************************************************************************
 *                                                                              *
 * SO_b_ia                                                                      *
 *                                                                              *
 * Calculates b_i^a from cis_vec , i and a are Spin orbitals                    *
 * cis_vec  is the cis vector in CSF-Basis (for one specific eigenvalue)        *
 *                                                                              *
 *******************************************************************************/


double SO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec){
  int I =  SO_to_MO(i);
  int A =  SO_to_MO(a);

  //CHECK whether coefficent belongs to a pure triplett -> if yes return zero
  if(I*A < 0) return(0.);
  
  //Start with sqrt(.5) for splitting up the CSF  
  //NO MINUS sign here !!!
  double coeff = sqrt(.5);

  //Reshift MOs
  I = abs(I)-1;
  A = abs(A)-1;
  
  //Calculate position in CIS-VEC
  int pos = 1 +                   // skip HF-ground state
           (I-llim)*umo  +        // set zero to lower limit used for correlation multiply with umo 
            A-llim-omo;           // add postion for virt orbs

  return(coeff*cis_vec[pos]);
}


/*******************************************************************************
 *                                                                              *
 * SO_v_ia                                                                      *
 *                                                                              *
 * Calculates v_i^a                                                             *
 *                                                                              *
 *******************************************************************************/

double SO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens){
  double v_ia = 0.;
  
  //Loop over k j b c (eq 13. CPL 219 (1994) 21-29
  for(int j = 2*llim; j < 2*(llim+omo); j++){
    for(int k = 2*llim; k < 2*(llim+omo); k++){
      for(int b =  2*(llim+omo); b < 2*(llim+omo+umo); b++){
	for(int c =  2*(llim+omo); c < 2*(llim+omo+umo); c++){
	
	  //clac leading factor
	  double jkbc = SO_calc_ijIIkl(j, k, b, c,
				       omo, umo, llim, prec_ints);
	  
	  //only proceed if != zero
	  double sum = 0.;
	  double coeff = 0.;
	  double a_term = 0.;
	  if(jkbc != 0.){	  
	    //!!!!! b_i_b * a_jk^ca
	    //first coeff
	    coeff = SO_b_ia(i, b,
			    omo,  umo,  llim, cis_vec);
	    // check before calc a
	    if(coeff != 0.){
	      a_term = SO_calc_a_ij_kl(j, k, c, a, 
				       omo, umo, llim, prec_ints, MOens);
	      sum += coeff*a_term;
	    }

	    //!!!!! b_j_a * a_ik^cb
	    //first coeff
	    coeff = SO_b_ia(j, a,  
			    omo,  umo,  llim, cis_vec);
	    // check before calc a
	    if(coeff != 0.){
	      a_term = SO_calc_a_ij_kl(i, k, c, b, 
				       omo, umo, llim, prec_ints, MOens);
	      sum += coeff*a_term;
	    }
	    
	    //!!!!! 2*b_j_b * a_jk^ac
	    //first coeff
	    coeff = SO_b_ia(j, b,  
			    omo,  umo,  llim, cis_vec);
	    // check before calc a
	    if(coeff != 0.){
	      a_term = SO_calc_a_ij_kl(i, k, a, c,
				       omo, umo, llim, prec_ints, MOens);
	      sum += 2*coeff*a_term;
	    }
	  }
	  
	  v_ia += 0.5*jkbc*sum;
	}
      }
    }
  }  
  return(v_ia);
}

/*******************************************************************************
 *                                                                              *
 * SO_calc_u_ij_ab                                                              *
 *                                                                              *
 * Calculates u_ij_ab from  eq 8    CPL 219 (1994) 21-29                        *
 *                                                                              *
 *******************************************************************************/

double SO_calc_u_ij_ab(int i, int j, int a, int b,  int omo, int umo, int llim, 
		       double* cis_vec, double* prec_ints){
  double u_ij_ab = 0.;
  
  //***********SUM OVER c*****************************************
  for(int c =  2*(llim+omo); c < 2*(llim+omo+umo); c++){
    //+ <ab||cj> b_i^c
    double b_ic = SO_b_ia(i, c,  
			  omo,  umo,  llim, cis_vec);
    if(b_ic != 0.)
      u_ij_ab += (SO_calc_ijIIkl( a, b, c, j,
				 omo, umo, llim, prec_ints)
		  *b_ic); 

    
    //- <ab||ci> b_j^c 
    double b_jc = SO_b_ia(j, c,  
			  omo,  umo,  llim, cis_vec);
    if(b_jc != 0.)
      u_ij_ab -= (SO_calc_ijIIkl( a, b, c, i, 
				  omo, umo, llim, prec_ints)
		  *b_jc);

 
  }
  
  //***********SUM OVER k*****************************************
  for(int k = 2*llim; k < 2*(llim+omo); k++){
    //+ <ka||ij> b_k^b
    double b_kb = SO_b_ia(k, b,  
			  omo,  umo,  llim, cis_vec);
    if(b_kb != 0.)
      u_ij_ab += (SO_calc_ijIIkl( k, a, i, j, 
				 omo, umo, llim, prec_ints)
		  *b_kb);

    
    //- <kb||ij> b_k^a
    double b_ka = SO_b_ia(k, a,  
			  omo,  umo,  llim, cis_vec);
    if(b_ka != 0.)
      u_ij_ab -= (SO_calc_ijIIkl( k, b, i, j, 
				 omo, umo, llim, prec_ints)
		  *b_ka);

  }
  
  return(u_ij_ab);
}


/*******************************************************************************
 *                                                                              *
 * SO_calc_t1                                                                   *
 *                                                                              *
 * Calculates firts sum in eq 14    CPL 219 (1994) 21-29                        *
 *                                                                              *
 *******************************************************************************/

double SO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en,
		  double* prec_ints, double* MOens){
  double t1 = 0.;
  //Loop over i j a b   eq 14. CPL 219 (1994) 21-29

//   for(int i = 2*llim; i < 2*(llim+omo); i++){
//     for(int j = 2*llim; j < 2*(llim+omo); j++){
//       for(int a =  2*(llim+omo); a < 2*(llim+omo+umo); a++){
// 	for(int b =  2*(llim+omo); b < 2*(llim+omo+umo); b++){
// 	  t1 += pow(SO_calc_u_ij_ab( i, j, a,  b,  omo, umo, llim, cis_vec, prec_ints),2)/
//                    (SO_calc_delta_ij_ab( i,  j,  a,  b, MOens)-cis_en);
//         }
//       }
//     }
//   }  
//   return(-0.25*t1);

  //MAKE IT FASTER
  for(int i = 2*llim; i < 2*(llim+omo); i++){
    for(int j = i; j < 2*(llim+omo); j++){
      for(int a =  2*(llim+omo); a < 2*(llim+omo+umo); a++){
	for(int b =  a; b < 2*(llim+omo+umo); b++){
	  t1 += pow(SO_calc_u_ij_ab( i, j, a,  b,  omo, umo, llim, cis_vec, prec_ints),2)/
                   (SO_calc_delta_ij_ab( i,  j,  a,  b, MOens)-cis_en);
	}
      }
    }
  }  
  return(-t1);

}



/*******************************************************************************
 *                                                                              *
 * CALC E2                                                                      *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the MP2 energy correction to the ground state                     *
 *                                                                              *
 *******************************************************************************/

double SO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens){
  double E2=0.;
  for(int i = 2*llim; i < 2*(llim+omo); i++){
    for(int j = 2*llim; j < 2*(llim+omo); j++){
      for(int a = 2*(omo+llim); a < 2*(llim+omo+umo); a++){
	for(int b = 2*(omo+llim); b < 2*(llim+omo+umo); b++){
	  E2 += (0.25*SO_calc_a_ij_kl( i, j, a, b, omo, umo, llim, prec_ints, MOens)*
		   SO_calc_ijIIkl(i, j, a, b, omo, umo, llim, prec_ints));
	}
      }
    }
  }
  return(E2);
}



/*******************************************************************************
 *                                                                              *
 * CALC D-correction                                                            *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the CIS(D) correction to CIS excitations energies                 *
 *                                                                              *
 *******************************************************************************/

double SO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double cis_en, double* t1, double* t2){
  
  *t1 = SO_calc_t1(omo, umo, llim, cis_vec, cis_en, prec_ints, MOens);
			      
  clog << "Hallo  " << *t1 << "\n";


  *t2 = 0.;
  for(int i = 2*llim; i < 2*(llim+omo); i++){
    for(int a = 2*(llim+omo); a < 2*(llim+omo+umo); a++){
      double coeff =  SO_b_ia( i,  a,  omo,   umo,   llim,  cis_vec);
      if(coeff != 0.){
   	coeff *= SO_calc_v_ia(i, a,  omo,  umo, llim, cis_vec, prec_ints,  MOens);
      }
      *t2 += coeff;
    }
  }
  return(*t1+*t2);
}
