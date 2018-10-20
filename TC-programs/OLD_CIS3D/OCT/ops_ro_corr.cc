/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_ro_corr.cc                                                         *
 *                                                                              *
 * operarotrs in "room"  orbital sapce  needed for cis(d) mp2 clacluations      * 
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
double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);
double RO_calc_a_ij_kl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_ints, double* MOens);
inline double RO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);
double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens);
double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim, double* cis_vec, double* prec_ints);
double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en, double* prec_ints, double* MOens);
double RO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens);
double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double cis_en, double* t1, double* t2);
void   SL_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2, 
		      int cis_size, ofstream* outf);
//Extern Functions
extern double get_precalc_ints_sd(int i, int j, int k, int l,
				  int omo, int umo, int llim, double* prec_ints);
extern double get_precalc_ints_ovov(int i, int j, int k, int l,
				    int omo, int umo, int llim, double* prec_ints);
extern void status(ofstream* outf);

/*******************************************************************************
 *                                                                              *
 * RO_calc_delta_ij_ab                                                          *
 *                                                                              *
 * 4 mo energie difference needed for MP2 and CIS(D)                            *
 *******************************************************************************/

double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens){
  return(-MOens[i]-MOens[j]
	 +MOens[a]+MOens[b]);
}
			   
						
					       

/*******************************************************************************
 *                                                                              *
 * RO_b_ia                                                                      *
 *                                                                              *
 * Calculates b_i^a from cis_vec , i and a are "room" orbitals                  *
 * cis_vec  is the cis vector in CSF-Basis (for one specific eigenvalue)        *
 *                                                                              *
 *******************************************************************************/


inline double RO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec){
  //Calculate position in CIS-VEC
  int pos = 1 +                   // skip HF-ground state
            (i-llim)*umo  +        // set zero to lower limit used for correlation multiply with umo 
            a-llim-omo;           // add postion for virt orbs

  return(cis_vec[pos]);
}


/******************************************************************************** 
 *                                                                              *
 * pos_i a                                                                      *
 *                                                                              *
 * Calculates pos of i^a in  cis_vec , i and a are "room" orbitals              *
 *                                                                              *
 *******************************************************************************/

int pos_i_a(int i, int a, int omo, int umo, int llim){
  int pos = 1 +                   // skip HF-ground state
           (i-llim)*umo  +        // set zero to lower limit used for correlation multiply with umo 
            a-llim-omo;           // add postion for virt orbs
  return(pos);
}


/*******************************************************************************
 *                                                                              *
 * RO_v_ia                                                                      *
 *                                                                              *
 * Calculates v_i^a                                                             *
 *                                                                              *
 *******************************************************************************/


double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* prec_ints, double* MOens){
  
  double X_IB = 0.;
  double X_JA = 0.;
  double X_JB = 0.;

  //Loop over k j b c (eq 13. CPL 219 (1994) 21-29
  //RO FORMULATION see DOC
  for(int j = llim; j < (llim+omo); j++){
    for(int k = llim; k < (llim+omo); k++){
      for(int b =  (llim+omo); b < (llim+omo+umo); b++){
	for(int c =  (llim+omo); c < (llim+omo+umo); c++){
	  //INIT prefacs
	  double F_IB = -RO_b_ia(i, b, omo,  umo,  llim, cis_vec)/RO_calc_delta_ij_ab( j,  k,  c,  a, MOens)/sqrt(2.);
	  double F_JA = -RO_b_ia(j, a, omo,  umo,  llim, cis_vec)/RO_calc_delta_ij_ab( i,  k,  c,  b, MOens)/sqrt(2.);
	  double F_JB = -RO_b_ia(j, b, omo,  umo,  llim, cis_vec)/RO_calc_delta_ij_ab( i,  k,  a,  c, MOens)*sqrt(2.);
	  
	  
	  double JKBC =  get_precalc_ints_ovov(j, b, k, c, omo, umo, llim, prec_ints);     //<JK|BC> = (JB|KC) => (JB|KC) -> J<B;K<C
	  double JKCB =  get_precalc_ints_ovov(j, c, k, b, omo, umo, llim, prec_ints);     //<JK|CB> = (JC|KB) => (JC|KB) -> J<C;K<B

	  double CAJK =  get_precalc_ints_ovov(j, c, k, a, omo, umo, llim, prec_ints);     //<CA|JK> = (CJ|AK) => (JC|KA) -> J<C;K<A
	  double CAKJ =  get_precalc_ints_ovov(k, c, j, a, omo, umo, llim, prec_ints);     //<CA|KJ> = (CK|AJ) => (KC|JA) -> K<C;J<A

	  double CBIK =  get_precalc_ints_ovov(i, c, k, b, omo, umo, llim, prec_ints);     //<CB|IK> = (CI|BK) => (IC|KB) -> I<C;K<B
	  double CBKI =  get_precalc_ints_ovov(k, c, i, b, omo, umo, llim, prec_ints);     //<CB|KI> = (CK|BI) => (KC|IB) -> K<C;I<B

	  double ACIK =  get_precalc_ints_ovov(i, a, k, c, omo, umo, llim, prec_ints);     //<AC|IK> = (AI|CK) => (IA|KC) -> I<A;K<C
 	  double ACKI =  get_precalc_ints_ovov(k, a, i, c, omo, umo, llim, prec_ints);     //<AC|KI> = (AK|CI) => (KA|IC) -> K<A;I<C

	  X_IB += F_IB*(JKBC*(CAJK-2.*CAKJ)+JKCB*(CAKJ-2.*CAJK));
	  X_JA += F_JA*(JKBC*(CBIK-2.*CBKI)+JKCB*(CBKI-2.*CBIK));
	  
	  X_JB += F_JB*(2*JKBC-JKCB)*(2*ACIK-ACKI);

	}
      }
    }
  }  
  return(0.5*(X_IB+X_JA+X_JB));
}

/*******************************************************************************
 *                                                                              *
 * RO_calc_u_ij_ab                                                              *
 *                                                                              *
 * Calculates u_ij_ab from  eq 8    CPL 219 (1994) 21-29                        *
 * Closed shell formulations see DOC 
 *                                                                              *
 *******************************************************************************/

double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim, 
		      double* cis_vec, double* prec_ints){

  double Cup = 0.;
  double Cum = 0.;
  
  double Kup = 0.;
  double Kum = 0.;

  
  //***********SUM OVER c*****************************************
  for(int c = (llim+omo); c < (llim+omo+umo); c++){
    double b_IC = RO_b_ia(i, c, omo,  umo,  llim, cis_vec);    
    double b_JC = RO_b_ia(j, c, omo,  umo,  llim, cis_vec);

    double ABIC = get_precalc_ints_sd(a, i, b, c, omo, umo, llim, prec_ints);     //<AB|IC> = (AI|BC)
    double ABCI = get_precalc_ints_sd(a, c, b, i, omo, umo, llim, prec_ints);     //<AB|CI> = (AC|BI)
    double ABJC = get_precalc_ints_sd(a, j, b, c, omo, umo, llim, prec_ints);     //<AB|JC> = (AJ|BC)
    double ABCJ = get_precalc_ints_sd(a, c, b, j, omo, umo, llim, prec_ints);     //<AB|CJ> = (AC|BJ)
    
    Cup += ABIC*b_JC+ABCJ*b_IC;
    Cum += ABCI*b_JC+ABJC*b_IC;
  }

  Cup *= 1./sqrt(2.);
  Cum *= 1./sqrt(2.);
  
  //***********SUM OVER k*****************************************
  for(int k = llim; k < (llim+omo); k++){
    double b_KB = RO_b_ia(k, b, omo,  umo,  llim, cis_vec);
    double b_KA = RO_b_ia(k, a, omo,  umo,  llim, cis_vec);

    double KAIJ = get_precalc_ints_sd(k, i, a, j, omo, umo, llim, prec_ints);     //<KA|IJ> = (KI|AJ)
    double KAJI = get_precalc_ints_sd(k, j, a, i, omo, umo, llim, prec_ints);     //<KA|JI> = (KJ|AI)
    double KBIJ = get_precalc_ints_sd(k, i, b, j, omo, umo, llim, prec_ints);     //<KB|IJ> = (KI|BJ)
    double KBJI = get_precalc_ints_sd(k, j, b, i, omo, umo, llim, prec_ints);     //<KB|JI> = (KJ|BI)

    Kup += KAIJ*b_KB+KBJI*b_KA;
    Kum += KAJI*b_KB+KBIJ*b_KA;
  }
  
  Kup *= 1./sqrt(2.);
  Kum *= 1./sqrt(2.);
  
  double u_ro = 2.*(pow(Cup-Cum+Kup-Kum,2)+pow(Cup-Kum,2)+pow(Kup-Cum,2));
  
  return(u_ro);
}


/*******************************************************************************
 *                                                                              *
 * RO_calc_t1                                                                   *
 *                                                                              *
 * Calculates firts sum in eq 14    CPL 219 (1994) 21-29                        *
 *                                                                              *
 *******************************************************************************/

double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en,
		  double* prec_ints, double* MOens){
  double t1 = 0.;

  //MAKE IT FASTER
  
  //i,i a,a terms unique
  for(int i = llim; i < (llim+omo); i++){
    for(int a =  (llim+omo); a < (llim+omo+umo); a++){
      t1 += 0.25*RO_calc_u_term( i, i, a,  a,  omo, umo, llim, cis_vec, prec_ints)/
	(RO_calc_delta_ij_ab( i,  i,  a,  a, MOens)-cis_en);
    }
  }

  //i,i,a,b terms = i,i,b,a =        x2
  for(int i = llim; i < (llim+omo); i++){
    for(int a =  (llim+omo); a < (llim+omo+umo); a++){
      for(int b =  a+1; b < (llim+omo+umo); b++){
	t1 += 0.5*RO_calc_u_term( i, i, a,  b,  omo, umo, llim, cis_vec, prec_ints)/
	  (RO_calc_delta_ij_ab( i,  i,  a,  b, MOens)-cis_en);
      }
    }
  }

  //i,j,a,a terms = j,i,a,a           => x2
  for(int i = llim; i < (llim+omo); i++){
    for(int j = i+1; j < (llim+omo); j++){
      for(int a =  (llim+omo); a < (llim+omo+umo); a++){
	t1 += 0.5*RO_calc_u_term( i, j, a,  a,  omo, umo, llim, cis_vec, prec_ints)/
	  (RO_calc_delta_ij_ab( i,  j,  a,  a, MOens)-cis_en);
      }
    }
  }


  //i,j,a,b terms = i,j,b,a ; j,i,a,b ; j,i,b,a    => x4
  for(int i = llim; i < (llim+omo); i++){
    for(int j = i+1; j < (llim+omo); j++){
      for(int a =  (llim+omo); a < (llim+omo+umo); a++){
	for(int b =  a+1; b < (llim+omo+umo); b++){
	  t1 +=   RO_calc_u_term( i, j, a,  b,  omo, umo, llim, cis_vec, prec_ints)/
	    (RO_calc_delta_ij_ab( i,  j,  a,  b, MOens)-cis_en);
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

double RO_calc_MP2(int omo, int umo, int llim, double* prec_ints, double* MOens){
  double E2=0.;
  for(int i = llim; i < (llim+omo); i++){
    for(int j = llim; j < (llim+omo); j++){
      for(int a = (omo+llim); a < (llim+omo+umo); a++){
	for(int b = (omo+llim); b < (llim+omo+umo); b++){
	  double J = get_precalc_ints_sd(i, a, j, b,
					 omo, umo, llim, prec_ints);
	  double K = get_precalc_ints_sd(i, b, j, a,
					 omo, umo, llim, prec_ints);
	  E2 += -0.25/RO_calc_delta_ij_ab( i,  j,  a,  b, MOens)*
	         2.*(pow(J-K,2)+J*J+K*K);
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

double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		      double* cis_vec, double cis_en, double* t1, double* t2){
  
  *t1 = RO_calc_t1(omo, umo, llim, cis_vec, cis_en, prec_ints, MOens);

  *t2 = 0.;
  for(int i = llim; i < (llim+omo); i++){
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
      double coeff =  RO_b_ia( i,  a,  omo,   umo,   llim,  cis_vec);
      coeff *= sqrt(2.)*RO_calc_v_ia(i, a,  omo,  umo, llim, cis_vec, prec_ints,  MOens);
      *t2 += coeff;
    }
  }
  return(*t1+*t2);
}



/*******************************************************************************
 *                                                                              *
 * CALC D-correction single loop                                                *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the CIS(D) correction to CIS excitations energies in one singel   *
 * loop for all states ( hopefully faster!!)                                    *
 *******************************************************************************/

void SL_calc_cis_dc(int omo, int umo, int llim, double* prec_ints, double* MOens,
		    double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2, 
		    int cis_size, ofstream* outf){
  
  
  //distribute tmp memory
  double* coeff = &(cistmpmat[0*cis_size]);
  
  double* X_IB  = &(cistmpmat[1*cis_size]);
  double* X_JA  = &(cistmpmat[2*cis_size]);
  double* X_JB  = &(cistmpmat[3*cis_size]);

  double* F_IB  = &(cistmpmat[4*cis_size]);
  double* F_JA  = &(cistmpmat[5*cis_size]);
  double* F_JB  = &(cistmpmat[6*cis_size]);
 

  //NO SQRT
  const double sq2 = sqrt(2.);



  *outf << "Single loop D-correction clacluation\n";
  status(outf);
  *outf << "Calculating doubles corrections\n";
  *outf << "...............................................................................\n";
  outf->flush();
  t1[0] = 0.;
  for(int x = 1; x < cis_size; x++){
    t1[x] = RO_calc_t1(omo, umo, llim, &(cis_vec[x*cis_size]), cis_en[x], prec_ints, MOens);
    *outf  << x << "\t";   
    if((x)%10==0) *outf << "\n";
    outf->flush(); 
  }
  
  //Make t2 zero
  for(int x = 0; x < cis_size; x++)
    t2[x] = 0.;
  
  *outf << "\n";
  status(outf);
  *outf << "Calculating tripels corrections\n";
  *outf << "...............................................................................\n";
  outf->flush();
  
  int count_det = 1;
  for(int i = llim; i < (llim+omo); i++){
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
      *outf  << count_det++ << "\t";   
      if((count_det)%10==0) *outf << "\n";
      outf->flush(); 
      
      int pos_ia =  pos_i_a(i, a, omo,  umo,  llim);
      
      for(int x = 1; x < cis_size; x++){
	coeff[x] = cis_vec[x*cis_size+pos_ia];
	X_IB[x]  = 0.;
	X_JA[x]  = 0.;
	X_JB[x] = 0.;
      }
      
      //Loop over k j b c (eq 13. CPL 219 (1994) 21-29
      //RO FORMULATION see DOC
      for(int j = llim; j < (llim+omo); j++){
	int pos_ja =  pos_i_a(j, a, omo,  umo,  llim);
	for(int k = llim; k < (llim+omo); k++){
	  for(int b =  (llim+omo); b < (llim+omo+umo); b++){
	    int pos_ib =  pos_i_a(i, b, omo,  umo,  llim);
	    int pos_jb =  pos_i_a(j, b, omo,  umo,  llim);
	    for(int c =  (llim+omo); c < (llim+omo+umo); c++){

	      //INIT prefacs
	      for(int x = 1; x < cis_size; x++){
		F_IB[x] = -cis_vec[x*cis_size+pos_ib]/RO_calc_delta_ij_ab( j,  k,  c,  a, MOens)/sq2;
		F_JA[x] = -cis_vec[x*cis_size+pos_ja]/RO_calc_delta_ij_ab( i,  k,  c,  b, MOens)/sq2;
		F_JB[x] = -cis_vec[x*cis_size+pos_jb]/RO_calc_delta_ij_ab( i,  k,  a,  c, MOens)*sq2;
	      }
	  
	      double JKBC =  get_precalc_ints_ovov(j, b, k, c, omo, umo, llim, prec_ints);     //<JK|BC> = (JB|KC) => (JB|KC) -> J<B;K<C
	      double JKCB =  get_precalc_ints_ovov(j, c, k, b, omo, umo, llim, prec_ints);     //<JK|CB> = (JC|KB) => (JC|KB) -> J<C;K<B

	      double CAJK =  get_precalc_ints_ovov(j, c, k, a, omo, umo, llim, prec_ints);     //<CA|JK> = (CJ|AK) => (JC|KA) -> J<C;K<A
	      double CAKJ =  get_precalc_ints_ovov(k, c, j, a, omo, umo, llim, prec_ints);     //<CA|KJ> = (CK|AJ) => (KC|JA) -> K<C;J<A
	      
	      double CBIK =  get_precalc_ints_ovov(i, c, k, b, omo, umo, llim, prec_ints);     //<CB|IK> = (CI|BK) => (IC|KB) -> I<C;K<B
	      double CBKI =  get_precalc_ints_ovov(k, c, i, b, omo, umo, llim, prec_ints);     //<CB|KI> = (CK|BI) => (KC|IB) -> K<C;I<B
	      
	      double ACIK =  get_precalc_ints_ovov(i, a, k, c, omo, umo, llim, prec_ints);     //<AC|IK> = (AI|CK) => (IA|KC) -> I<A;K<C
	      double ACKI =  get_precalc_ints_ovov(k, a, i, c, omo, umo, llim, prec_ints);     //<AC|KI> = (AK|CI) => (KA|IC) -> K<A;I<C
	      
	      for(int x = 1; x < cis_size; x++){
		X_IB[x] += F_IB[x]*(JKBC*(CAJK-2.*CAKJ)+JKCB*(CAKJ-2.*CAJK));
		X_JA[x] += F_JA[x]*(JKBC*(CBIK-2.*CBKI)+JKCB*(CBKI-2.*CBIK));
		
		X_JB[x] += F_JB[x]*(2*JKBC-JKCB)*(2*ACIK-ACKI);
	      }
	      
	    }
	  }
	}
      }

      for(int x = 1; x < cis_size; x++){      
	coeff[x] *= sqrt(2.)*.5*(X_IB[x]+X_JA[x]+X_JB[x]);
	t2[x] += coeff[x];
      }
 
    } 
  }
}
