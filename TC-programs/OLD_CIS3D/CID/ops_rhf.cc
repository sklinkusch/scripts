/********************************************************************************
 * Functions for CIS example program with determinant expressions               *
 *                                                                              *
 * file: ops_rhf.cc                                                             *
 *                                                                              *
 * contains functions for some nuclei properties                                *
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
 ********************************************************************************/ 
#include <math.h>

using namespace std;


//Functions
double calc_r_ab(int a, int b, double* coord);
double calc_ion_rep(int nroa, double* coord, double* charges);
void   calc_center_of_mass(int nroa, double* coord, double* mass, double* center_of_mass);
void   calc_mu_core(int nroa, double* coord, double* charges, double* point, 
		    double* mu_core);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_r_ab                                                       */
/*                                                                               */
/* caluculates the ion core coulomb repulsion                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_r_ab(int a, int b, double* coord){
  double r = sqrt(pow(coord[3*a+0]-coord[3*b+0],2)+
		  pow(coord[3*a+1]-coord[3*b+1],2)+
		  pow(coord[3*a+2]-coord[3*b+2],2));
  return(r);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*               calc_ion_rep                                                    */
/*                                                                               */
/* calculates the ion core coulomb repulsion                                     */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double calc_ion_rep(int nroa, double* coord, double* charges){
  double ion_rep = 0.;
  for(int x = 0; x < nroa; x++){
    for(int y = x+1; y < nroa; y++){
      ion_rep += charges[x]*charges[y]/calc_r_ab(x,y,coord);
    }
  }

  return(ion_rep);
}


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

