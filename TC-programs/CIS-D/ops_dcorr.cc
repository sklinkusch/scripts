/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_ro_corr.cc                                                         *
 *                                                                              *
 * operators in "room"  orbital sapce  needed for cis(d) mp2 calculations       * 
 *                                                                              *
 *  2 ELECTRON INTEGRALS IN PHYSICIST NOTATION IN THIS FILE !!!!!!!!!!!         *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>

using namespace std;


//Functions
double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens);
double RO_calc_a_ij_kl(int i, int j, int k, int l, int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens);
inline double RO_b_ia(int i, int a, int omo, int umo, int llim, double* cis_vec);
double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* twoel, double* MOdelta);
double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim, double* cis_vec, double* ints_ooov, double* ints_ovvv);
double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en, double* ints_ooov, double* ints_ovvv, double* MOdelta);
double RO_calc_MP2(int omo, int umo, int llim, double* ints_ovov, double* MOdelta);
//double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
//                      double* cis_vec, double cis_en, double* t1, double* t2);
/*void   SL_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
                      double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2,
                      int cis_size, ofstream* outf);*/
void calc_MOdelta(int omo, int umo, int llim, double* MOdelta, double* MOens);
long long int map_oovv(int i, int j, int a, int b, int omo, int umo, int llim);
long long int map_ovov(int i, int a, int j, int b, int omo, int umo, int llim);
long long int map_ooov(int i, int j, int k, int a, int omo, int umo, int llim);
void map_ints_ooov(int omo, int umo, int llim, double* ints_ooov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
long long int map_ovvv(int i, int a, int b, int c, int omo, int umo, int llim);
void map_ints_ovvv(int omo, int umo, int llim, double* ints_ovvv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
void map_ints_ovov(int omo, int umo, int llim, double* ints_ovov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
void map_ints_oovv(int omo, int umo, int llim, double* ints_oovv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
void SK_calc_cis_dc(int x, int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOdelta, double* cisvec, double* cisen, int cis_size, double* t1);
void SK_calc_cis_tc(int x, int omo, int umo, int llim, double* cisvec, int cis_size, double *t2, double* twoel, double* MOdelta);
//Extern Functions
extern double get_precalc_ints_sd(int i, int j, int k, int l,
                                  int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
extern double get_precalc_ints_ovov(int i, int j, int k, int l,
                                    int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows);
extern void status(ofstream* outf);
extern double get_matrix_element(long long int row, long long int column, double* nonzero, uint32_t* col_ind, uint32_t* row_ptr);

/*******************************************************************************
 *                                                                              *
 * RO_calc_delta_ij_ab                                                          *
 *                                                                              *
 * 4 mo energie difference needed for MP2 and CIS(D)                            *
  *******************************************************************************/

double RO_calc_delta_ij_ab(int i, int j, int a, int b, double* MOens){
 return(-MOens[i]-MOens[j]+MOens[a]+MOens[b]);
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
 int pos = 1 +                    // skip HF-ground state
           (i-llim)*umo  +        // set zero to lower limit used for correlation multiply with umo 
            a-llim-omo;           // add postion for virt orbs
 return(pos);
}

/********************************************************************************
 *                                                                              *
 * RO_v_ia                                                                      *
 *                                                                              *
 * Calculates v_i^a                                                             *
 *                                                                              *
 ********************************************************************************/

double RO_calc_v_ia(int i, int a, int omo, int umo, int llim, double* cis_vec, double* ints_ovov, double* MOdelta){
 // define individual terms of v_i^a and set them to 0.
 double X_IB = 0.;
 double X_JA = 0.;
 double X_JB = 0.;
 int mapvar;
 //Loop over k j b c (eq 13. CPL 219 (1994) 21-29
 //RO FORMULATION see DOC
#pragma omp parallel for reduction (+: X_IB, X_JA, X_JB)
 for(int c = (llim+omo); c < (llim+omo+umo); c++){
  for(int k =  llim; k < (llim+omo); k++){
   mapvar = map_ovov(i, a, k, c, omo, umo, llim);
   double ACIK =  ints_ovov[mapvar];     //<AC|IK> = (AI|CK) => (IA|KC) -> I<A;K<C
   mapvar = map_ovov(k, a, i, c, omo, umo, llim);
   double ACKI =  ints_ovov[mapvar];     //<AC|KI> = (AK|CI) => (KA|IC) -> K<A;I<C
   for(int j = llim; j < (llim+omo); j++){                                          // four fold loop
    mapvar = map_ovov(j, c, k, a, omo, umo, llim);
    double CAJK =  ints_ovov[mapvar];     //<CA|JK> = (CJ|AK) => (JC|KA) -> J<C;K<A
    mapvar = map_ovov(k, c, j, a, omo, umo, llim);
    double CAKJ =  ints_ovov[mapvar];     //<CA|KJ> = (CK|AJ) => (KC|JA) -> K<C;J<A
    for(int b =  (llim+omo); b < (llim+omo+umo); b++){
     //calculate prefactors of two-electron integrals
     mapvar = map_oovv(j, k, c, a, omo, umo, llim);
     double F_IB = -RO_b_ia(i, b, omo,  umo,  llim, cis_vec)/MOdelta[mapvar]/sqrt(2.);
     mapvar = map_oovv(i, k, c, b, omo, umo, llim);
     double F_JA = -RO_b_ia(j, a, omo,  umo,  llim, cis_vec)/MOdelta[mapvar]/sqrt(2.);
     mapvar = map_oovv(i, k, a, c, omo, umo, llim);
     double F_JB = -RO_b_ia(j, b, omo,  umo,  llim, cis_vec)/MOdelta[mapvar]*sqrt(2.);

     // get precalculated two-electron integrals
     mapvar = map_ovov(j, b, k, c, omo, umo, llim);
     double JKBC =  ints_ovov[mapvar];     //<JK|BC> = (JB|KC) => (JB|KC) -> J<B;K<C
     mapvar = map_ovov(j, c, k, b, omo, umo, llim);
     double JKCB =  ints_ovov[mapvar];     //<JK|CB> = (JC|KB) => (JC|KB) -> J<C;K<B

     mapvar = map_ovov(i, c, k, b, omo, umo, llim);
     double CBIK =  ints_ovov[mapvar];     //<CB|IK> = (CI|BK) => (IC|KB) -> I<C;K<B
     mapvar = map_ovov(k, c, i, b, omo, umo, llim);
     double CBKI =  ints_ovov[mapvar];     //<CB|KI> = (CK|BI) => (KC|IB) -> K<C;I<B

     // calculate individual terms of v_i^a as the sum of the antisymmetrized two-electron integrals weighed by the prefactors
     X_IB += F_IB*(JKBC*(CAJK-2.*CAKJ)+JKCB*(CAKJ-2.*CAJK));
     X_JA += F_JA*(JKBC*(CBIK-2.*CBKI)+JKCB*(CBKI-2.*CBIK));
     X_JB += F_JB*(2*JKBC-JKCB)*(2*ACIK-ACKI);
    }
   }
  }
 }
 return(0.5*(X_IB+X_JA+X_JB));
}

/********************************************************************************
 *                                                                              *
 * RO_calc_u_ij_ab                                                              *
 *                                                                              *
 * Calculates u_ij_ab from  eq 8    CPL 219 (1994) 21-29                        *
 * Closed shell formulations see DOC                                            *
 *                                                                              *
 ********************************************************************************/

double RO_calc_u_term(int i, int j, int a, int b,  int omo, int umo, int llim,
                      double* cis_vec, double* ints_ooov, double* ints_ovvv){

 double Cup = 0.;
 double Cum = 0.;
 double Kup = 0.;
 double Kum = 0.;
 long long int mapvar;

 //***********SUM OVER c*****************************************
 for(int c = (llim+omo); c < (llim+omo+umo); c++){
  // calculate prefactors
  double b_IC = RO_b_ia(i, c, omo,  umo,  llim, cis_vec);
  double b_JC = RO_b_ia(j, c, omo,  umo,  llim, cis_vec);
 // get precalculated integrals 
  mapvar = map_ovvv(i, a, b, c, omo, umo, llim);
  double ABIC = ints_ovvv[mapvar]; // get_precalc_ints_sd(a, i, b, c, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<AB|IC> = (AI|BC)
  mapvar = map_ovvv(i, b, a, c, omo, umo, llim);
  double ABCI = ints_ovvv[mapvar]; // get_precalc_ints_sd(a, c, b, i, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<AB|CI> = (AC|BI)
  mapvar = map_ovvv(j, a, b, c, omo, umo, llim);
  double ABJC = ints_ovvv[mapvar]; // get_precalc_ints_sd(a, j, b, c, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<AB|JC> = (AJ|BC)
  mapvar = map_ovvv(j, b, a, c, omo, umo, llim);
  double ABCJ = ints_ovvv[mapvar]; // get_precalc_ints_sd(a, c, b, j, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<AB|CJ> = (AC|BJ)
 // sum up the two terms of the "c" sum
  Cup += ABIC*b_JC+ABCJ*b_IC;
  Cum += ABCI*b_JC+ABJC*b_IC;
 }
 // calculate final values for the "c" sum terms
 Cup *= 1./sqrt(2.);
 Cum *= 1./sqrt(2.);

 //***********SUM OVER k*****************************************
 for(int k = llim; k < (llim+omo); k++){
  // calculate prefactors
  double b_KB = RO_b_ia(k, b, omo,  umo,  llim, cis_vec);
  double b_KA = RO_b_ia(k, a, omo,  umo,  llim, cis_vec);
  // get precalculated integrals
  mapvar = map_ooov(k, i, j, a, omo, umo, llim);
  double KAIJ = ints_ooov[mapvar]; // get_precalc_ints_sd(k, i, a, j, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<KA|IJ> = (KI|AJ)
  mapvar = map_ooov(k, j, i, a, omo, umo, llim);
  double KAJI = ints_ooov[mapvar]; // get_precalc_ints_sd(k, j, a, i, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<KA|JI> = (KJ|AI)
  mapvar = map_ooov(k, i, j, b, omo, umo, llim);
  double KBIJ = ints_ooov[mapvar]; // get_precalc_ints_sd(k, i, b, j, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<KB|IJ> = (KI|BJ)
  mapvar = map_ooov(k, j, i, b, omo, umo, llim);
  double KBJI = ints_ooov[mapvar]; // get_precalc_ints_sd(k, j, b, i, omo, umo, llim, prec_vals, prec_cols, prec_rows);     //<KB|JI> = (KJ|BI)
  // sum up the two terms of the "k" sum
  Kup += KAIJ*b_KB+KBJI*b_KA;
  Kum += KAJI*b_KB+KBIJ*b_KA;
 }
 // calculate final values for the "k" sum terms
 Kup *= 1./sqrt(2.);
 Kum *= 1./sqrt(2.);
 // calculate doubles corrections
 double u_ro = 2.*(pow(Cup-Cum+Kup-Kum,2)+pow(Cup-Kum,2)+pow(Kup-Cum,2));
 return(u_ro);
}

/********************************************************************************
 *                                                                              *
 * RO_calc_t1                                                                   *
 *                                                                              *
 * Calculates firts sum in eq 14    CPL 219 (1994) 21-29                        *
 *                                                                              *
 ********************************************************************************/

double RO_calc_t1(int omo, int umo, int llim, double* cis_vec, double cis_en,
	                  double* ints_ooov, double* ints_ovvv, double* MOdelta){
 int mapvar;
 double t1 = 0.;
 //MAKE IT FASTER
 //i,i a,a terms unique and i,i,a,b = x2
#pragma omp parallel for reduction (+:t1)
 for(int i = llim; i < (llim+omo); i++){
  for(int a =  (llim+omo); a < (llim+omo+umo); a++){
      mapvar = map_oovv(i, i, a, a, omo, umo, llim);
    t1 += 0.25*RO_calc_u_term( i, i, a,  a,  omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
          (MOdelta[mapvar]-cis_en);
    for(int b = a+1; b < (llim+omo+umo); b++){
	mapvar = map_oovv(i, i, a, b, omo, umo, llim);
	t1 += 0.5*RO_calc_u_term(i, i, a, b, omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
	    (MOdelta[mapvar]-cis_en);
    }
  }
 }
 /*
 //i,i,a,b terms = i,i,b,a =        x2
 for(int i = llim; i < (llim+omo); i++){
  for(int a =  (llim+omo); a < (llim+omo+umo); a++){
   for(int b =  a+1; b < (llim+omo+umo); b++){
       mapvar = map_oovv(i, i, a, b, omo, umo, llim);
    t1 += 0.5*RO_calc_u_term( i, i, a,  b,  omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
          (MOdelta[mapvar]-cis_en);
   }
  }
 }
 */
 //i,j,a,a terms = j,i,a,a           => x2 and i, j, a, b => x4
#pragma omp parallel for reduction (+:t1)
 for(int i = llim; i < (llim+omo); i++){
  for(int j = i+1; j < (llim+omo); j++){
   for(int a =  (llim+omo); a < (llim+omo+umo); a++){
       mapvar = map_oovv(i, j, a, a, omo, umo, llim);
    t1 += 0.5*RO_calc_u_term( i, j, a,  a,  omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
          (MOdelta[mapvar]-cis_en);
    for(int b = a+1; b < (llim+omo+umo); b++){
	mapvar = map_oovv(i, j, a, b, omo, umo, llim);
	t1 += RO_calc_u_term(i, j, a, b, omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
	    (MOdelta[mapvar]-cis_en);
    }
   }
  }
 }
 /*
 //i,j,a,b terms = i,j,b,a ; j,i,a,b ; j,i,b,a    => x4
 for(int i = llim; i < (llim+omo); i++){
  for(int j = i+1; j < (llim+omo); j++){
   for(int a =  (llim+omo); a < (llim+omo+umo); a++){
    for(int b =  a+1; b < (llim+omo+umo); b++){
	mapvar = map_oovv(i, j, a, b, omo, umo, llim);
     t1 +=   RO_calc_u_term( i, j, a,  b,  omo, umo, llim, cis_vec, ints_ooov, ints_ovvv)/
             (MOdelta[mapvar]-cis_en);
    }
   }
  }
 }
 */
 return(-t1);
}

/********************************************************************************
 *                                                                              *
 * CALC E2                                                                      *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the MP2 energy correction to the ground state                     *
 *                                                                              *
 ********************************************************************************/

double RO_calc_MP2(int omo, int umo, int llim, double* ints_ovov, double* MOdelta){
 double E2=0.;
 long long int mapvar;
#pragma omp parallel for reduction(+:E2)
 for(int i = llim; i < (llim+omo); i++){
  for(int j = llim; j < (llim+omo); j++){
   for(int a = (omo+llim); a < (llim+omo+umo); a++){
    for(int b = (omo+llim); b < (llim+omo+umo); b++){
     mapvar = map_ovov(i, a, j, b, omo, umo, llim);
     double J = ints_ovov[mapvar]; // get_precalc_ints_sd(i, a, j, b, omo, umo, llim, prec_vals, prec_cols, prec_rows);
     mapvar = map_ovov(i, b, j, a, omo, umo, llim);
     double K = ints_ovov[mapvar]; // get_precalc_ints_sd(i, b, j, a, omo, umo, llim, prec_vals, prec_cols, prec_rows);
     mapvar = map_oovv(i, j, a, b, omo, umo, llim);
     E2 += -0.25/MOdelta[mapvar]*2.*(pow(J-K,2)+J*J+K*K);
    }
   }
  }
 }
 return(E2);
}

/********************************************************************************
 *                                                                              *
 * CALC D-correction                                                            *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the CIS(D) correction to CIS excitations energies                 *
 *                                                                              *
 ********************************************************************************/
/*
double RO_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
                      double* cis_vec, double cis_en, double* t1, double* t2){
 *t1 = RO_calc_t1(omo, umo, llim, cis_vec, cis_en, prec_vals, prec_cols, prec_rows, MOens);
 *t2 = 0.;
 for(int i = llim; i < (llim+omo); i++){
  for(int a = (llim+omo); a < (llim+omo+umo); a++){
   double coeff =  RO_b_ia( i,  a,  omo,   umo,   llim,  cis_vec);
   coeff *= sqrt(2.)*RO_calc_v_ia(i, a,  omo,  umo, llim, cis_vec, prec_vals, prec_cols, prec_rows,  MOens);
   *t2 += coeff;
  }
 }
 return(*t1+*t2);
}
*/
/********************************************************************************
 *                                                                              *
 * CALC D-correction single loop                                                *
 *                                  CPL 219 (1994) 21-29                        *
 * Calculates the CIS(D) correction to CIS excitations energies in one single   *
 * loop for all states ( hopefully faster!!)                                    *
 ********************************************************************************/
/*
void SL_calc_cis_dc(int omo, int umo, int llim, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows, double* MOens,
                    double* cis_vec, double* cis_en, double* cistmpmat, double* t1, double* t2,
                    int cis_size, ofstream* outf){
 // distribute tmp memory
// double* coeff = &(cistmpmat[0*cis_size]);
 double* coeff = new double[cis_size];

 double* X_IB  = new double[cis_size]; //&(cistmpmat[1*cis_size]);
 double* X_JA  = new double[cis_size]; //&(cistmpmat[2*cis_size]);
 double* X_JB  = new double[cis_size]; //&(cistmpmat[3*cis_size]);

 double* F_IB  = new double[cis_size]; //&(cistmpmat[4*cis_size]);
 double* F_JA  = new double[cis_size]; //&(cistmpmat[5*cis_size]);
 double* F_JB  = new double[cis_size]; //&(cistmpmat[6*cis_size]);

 // NO SQRT
 const double sq2 = sqrt(2.);
 
 *outf << "Single loop D-correction calculation\n";
 status(outf);
 *outf << "Calculating doubles corrections\n";
 *outf << "...............................................................................\n";

 outf->flush();
// t1[0] = 0.;
 int count_doub = 1;
#pragma omp parallel for
 for(int x = 0; x < cis_size; x++){
  t1[x] = RO_calc_t1(omo, umo, llim, &(cis_vec[x*cis_size]), cis_en[x], prec_vals, prec_cols, prec_rows, MOens);
 *outf  << x << "\t";
 count_doub++;
  if((count_doub)%10==0) *outf << count_doub << "\n";
  outf->flush();
 }
 // Make t2 zero
#pragma omp parallel for
 for(int x = 0; x < cis_size; x++)
  t2[x] = 0.;
 *outf << "\n";
 status(outf);
 *outf << "Calculating triple corrections\n";
 *outf << "...............................................................................\n";
 outf->flush();
long long int dimens = calc_mapvar(llim+omo-1, llim+omo-1, llim+omo+umo-1, llim+omo+umo, omo, umo, llim);
double* twoel = new double[dimens];
double* MOdelta = new double[dimens];
calc_twoel_MOdelta(omo, umo, llim, twoel, MOdelta, MOens, prec_vals, prec_cols, prec_rows);
long long int mapvar;
 for(int a = (llim+omo); a < (llim+omo+umo); a++){
  for(int i = llim; i < (llim+omo); i++){
      int pos_ia = pos_i_a(i, a, omo, umo, llim);
#pragma omp parallel for
      for(int x = 0; x < cis_size; x++){
	  coeff[x] = cis_vec[x*cis_size+pos_ia];
	  X_IB[x] = 0.;
	  X_JA[x] = 0.;
	  X_JB[x] = 0.;
      }
      for(int c = (llim+omo); c < (llim+omo+umo); c++){
	  for(int k = llim; k < (llim+omo); k++){
	      mapvar = calc_mapvar(i, k, a, c, omo, umo, llim);
	      double ACIK = twoel[mapvar];
	      mapvar = calc_mapvar(k, i, a, c, omo, umo, llim);
	      double ACKI = twoel[mapvar];
	      for(int j = llim; j < (llim+omo); j++){
		  int pos_ja = pos_i_a(j, a, omo, umo, llim);
		  mapvar = calc_mapvar(j, k, c, a, omo, umo, llim);
		  double CAJK = twoel[mapvar];
		  mapvar = calc_mapvar(k, j, c, a, omo, umo, llim);
		  double CAKJ = twoel[mapvar];
		  for(int b = (llim+omo); b < (llim+omo+umo); b++){
		      int pos_ib = pos_i_a(i, b, omo, umo, llim);
		      int pos_jb = pos_i_a(j, b, omo, umo, llim);
		      mapvar = calc_mapvar(i, k, c, b, omo, umo, llim);
		      double CBIK = twoel[mapvar];
		      mapvar = calc_mapvar(k, i, c, b, omo, umo, llim);
		      double CBKI = twoel[mapvar];
		      mapvar = calc_mapvar(j, k, b, c, omo, umo, llim);
		      double JKBC = twoel[mapvar];
		      mapvar = calc_mapvar(j, k, c, b, omo, umo, llim);
		      double JKCB = twoel[mapvar];
		      for(int x = 0; x < cis_size; x++){
			  mapvar = calc_mapvar(j, k, c, a, omo, umo, llim);
			  F_IB[x] = -cis_vec[x*cis_size+pos_ib]/MOdelta[mapvar]/sq2;
			  mapvar = calc_mapvar(i, k, c, b, omo, umo, llim);
			  F_JA[x] = -cis_vec[x*cis_size+pos_ja]/MOdelta[mapvar]/sq2;
			  mapvar = calc_mapvar(i, k, a, c, omo, umo, llim);
			  F_JB[x] = -cis_vec[x*cis_size+pos_jb]/MOdelta[mapvar]*sq2;
			  X_IB[x] += F_IB[x]*(JKBC*(CAJK-2.*CAKJ)+JKCB*(CAKJ-2.*CAJK));
			  X_JA[x] += F_JA[x]*(JKBC*(CBIK-2.*CBKI)+JKCB*(CBKI-2.*CBIK));
			  X_JB[x] += F_JB[x]*(2.*JKBC-JKCB)*(2.*ACIK-ACKI);
			  }
		     }
		 }
	     }
	 }
      for(int x = 0; x < cis_size; x++){
	  coeff[x] *= 0.5*sq2*(X_IB[x]+X_JA[x]+X_JB[x]);
	  t2[x] += coeff[x];
      }
     }
 }
}*/

long long int map_oovv(int i, int j, int a, int b, int omo, int umo, int llim){
    long long int dumvar = ((long long int) i - (long long int) llim) * (long long int) omo * (long long int) umo * (long long int) umo
	+ ((long long int) j - (long long int) llim) * (long long int) umo * (long long int) umo
	+ ((long long int) a - (long long int) llim - (long long int) omo) * (long long int) umo
	+ ((long long int) b - (long long int) llim - (long long int) omo);
    return(dumvar);
}

long long int map_ovov(int i, int a, int j, int b, int omo, int umo, int llim){
    long long int dumvar = ((long long int) i - (long long int) llim) * (long long int) omo * (long long int) umo * (long long int) umo
	+ ((long long int) a - (long long int) llim - (long long int) omo) * (long long int) umo * (long long int) omo
	+ ((long long int) j - (long long int) llim) * (long long int) umo
	+ ((long long int) b - (long long int) llim - (long long int) omo);
    return(dumvar);
}

void calc_MOdelta(int omo, int umo, int llim, double* MOdelta, double* MOens){
    long long int mapvar;
#pragma omp parallel for
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
	for(int i = llim; i < (llim+omo); i++){
	    for(int j = llim; j < (llim+omo); j++){
		for(int b = (llim+omo); b < (llim+omo+umo); b++){
		    mapvar = map_oovv(i, j, a, b, omo, umo, llim);
		    MOdelta[mapvar] = RO_calc_delta_ij_ab(i, j, a, b, MOens);
		}
	    }
	}
    }
}

long long int map_ooov(int i, int j, int k, int a, int omo, int umo, int llim){
    long long int x = ((long long int) a - (long long int) llim - (long long int) omo)
	+ ((long long int) k - (long long int) llim) * (long long int) umo
	+ ((long long int) j - (long long int) llim) * (long long int) umo * (long long int) omo
	+ ((long long int) i - (long long int) llim) * (long long int) umo * (long long int) omo * (long long int) omo;
    return(x);
}

void map_ints_ooov(int omo, int umo, int llim, double* ints_ooov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows){
    long long int mapvar;
#pragma omp parallel for
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
	for(int i = llim; i < (llim+omo); i++){
	    for(int j = llim; j < (llim+omo); j++){
		for(int k = llim; k < (llim+omo); k++){
		    mapvar = map_ooov(i, j, k, a, omo, umo, llim);
		    ints_ooov[mapvar] = get_precalc_ints_sd(i, j, k, a, omo, umo, llim, prec_vals, prec_cols, prec_rows);
		}
	    }
	}
    }
}

long long int map_ovvv(int i, int a, int b, int c, int omo, int umo, int llim){
    long long int x = ((long long int) c - (long long int) llim - (long long int) omo) +
	((long long int) b - (long long int) llim - (long long int) omo) * (long long int) umo +
	((long long int) a - (long long int) llim - (long long int) omo) * (long long int) umo * (long long int) umo +
	((long long int) i - (long long int) llim) * (long long int) umo * (long long int) umo * (long long int) umo;
    return(x);
}

void map_ints_ovvv(int omo, int umo, int llim, double* ints_ovvv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows){
    long long int mapvar;
#pragma omp parallel for
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
	for(int i = llim; i < (llim+omo); i++){
	    for(int b = (llim+omo); b < (llim+omo+umo); b++){
		for(int c = (llim+omo); c < (llim+omo+umo); c++){
		    mapvar = map_ovvv(i, a, b, c, omo, umo, llim);
		    ints_ovvv[mapvar] = get_precalc_ints_sd(i, a, b, c, omo, umo, llim, prec_vals, prec_cols, prec_rows);
		}
	    }
	}
    }
}

void map_ints_ovov(int omo, int umo, int llim, double* ints_ovov, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows){
    long long int mapvar;
#pragma omp parallel for
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
	for(int i = llim; i < (llim+omo); i++){
	    for(int j = llim; j < (llim+omo); j++){
		for(int b = (llim+omo); b < (llim+omo+umo); b++){
		    mapvar = map_ovov(i, a, j, b, omo, umo, llim);
		    ints_ovov[mapvar] = get_precalc_ints_ovov(i, a, j, b, omo, umo, llim, prec_vals, prec_cols, prec_rows);
		}
	    }
	}
    }
}

void map_ints_oovv(int omo, int umo, int llim, double* ints_oovv, double* prec_vals, uint32_t* prec_cols, uint32_t* prec_rows){
    long long int mapvar;
#pragma omp parallel for
    for(int a = (llim+omo); a < (llim+omo+umo); a++){
	for(int i = llim; i < (llim+omo); i++){
	    for(int j = llim; j < (llim+omo); j++){
		for(int b = (llim+omo); b < (llim+omo+umo); b++){
		    mapvar = map_oovv(i, j, a, b, omo, umo, llim);
		    ints_oovv[mapvar] = get_precalc_ints_sd(i, j, a, b, omo, umo, llim, prec_vals, prec_cols, prec_rows);
		}
	    }
	}
    }
}

void SK_calc_cis_dc(int x, int omo, int umo, int llim, double* ints_ooov, double* ints_ovvv, double* MOdelta, double* cisvec, double* cisen, int cis_size, double* t1){
    *t1 = RO_calc_t1(omo, umo, llim, &(cisvec[x*cis_size]), cisen[x], ints_ooov, ints_ovvv, MOdelta);
}

void SK_calc_cis_tc(int x, int omo, int umo, int llim, double* cisvec, int cis_size, double *t2, double* twoel, double* MOdelta){
    double sq2 = sqrt(2.);
    *t2 = 0.;
    for(int i = llim; i < (llim+omo); i++){
	for(int a = (llim+omo); a < (llim+omo+umo); a++){
	    int pos_ia = pos_i_a(i, a, omo, umo, llim);
	    double coeff = cisvec[x*cis_size+pos_ia];
	    coeff *= sq2*RO_calc_v_ia(i, a, omo, umo, llim, &(cisvec[x*cis_size]), twoel, MOdelta);
	    *t2 += coeff;
	}
    }
}

