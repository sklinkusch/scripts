/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_rhf.cc                                                             *
 *                                                                              *
 * contains functions for (time indepentend)  operations  needed for RHF        *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/ 
#include <math.h>

using namespace std;


//Functions
void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
		  double* mu_core);
double calc_op_1el(int nroao, double* opmat, double* Pmat);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_center_of_mass                                             */
/*                                                                               */
/* calculates the ion core center of mass                                         */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass){
  for(int x = 0; x < 3; x++) center_of_mass[x] = 0.;
  
  double mass_tot = 0.;
  
  for(int x = 0; x < nroa; x++){
    mass_tot += mass[x];
    for(int i = 0; i < 3; i++)
      center_of_mass[i] += mass[x] * coord[3*x+i];
  }
  
  for(int i = 0; i < 3; i++)
    center_of_mass[i] /= mass_tot;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_mu_core                                                    */
/*                                                                               */
/* calculates the ion core dipole moment  at point                               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void calc_mu_core(int nroa, double* coord, double* charges, double* point, 
		  double* mu_core){
  for(int x = 0; x < 3; x++) mu_core[x] = 0.;
  
  
  for(int x = 0; x < nroa; x++){
    for(int i = 0; i < 3; i++)
      mu_core[i] += charges[x] * (coord[3*x+i] - point[i]);
  }
  
}

/*******************************************************************************
 * This routine builds the new Pmat in a damped SCF                            *
 *                                                                             *
 ******************************************************************************/

double build_Pmat_dscf(int nroao, int nroe, double* Pmat, double* Pmat_old, 
		       double* MOs, double damp){
  double max_diff = 0.;
  
  for(int x = 0; x < nroao*nroao; x++) Pmat[x] = 0.;
  
  for(int e = 0; e < nroe/2; e++){
    for(int x = 0; x < nroao; x++){
      for(int y = 0; y < nroao; y++){
        Pmat[x*nroao+y] += 2.*MOs[e*nroao+y]*MOs[e*nroao+x];
      }
    }
  }

  for(int x = 0; x < nroao*nroao; x++){
    if(fabs(Pmat_old[x] - Pmat[x]) > max_diff) max_diff = fabs(Pmat_old[x] - Pmat[x]);
    Pmat[x] = (1.-damp)*Pmat[x] + damp*Pmat_old[x];
    Pmat_old[x] = Pmat[x];
  }
  
  return(max_diff);
}

/*******************************************************************************
 * calculation of one electrion expc. values                                   *
 *                                                                             *
 ******************************************************************************/

double   calc_op_1el(int nroao, double* opmat, double* Pmat){
  double op = 0.;
  for(int x = 0; x < nroao;x++){
    for(int y = 0; y < nroao;y++){
      op += Pmat[x*nroao+y]*opmat[y*nroao+x];
    }
  }
  return(op);
}
