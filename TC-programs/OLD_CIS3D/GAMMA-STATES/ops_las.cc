/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: ops_las.cc                                                             *
 *                                                                              *
 * contains functions for laser pulses                                          *
 *                                                                              *
 *                                                    Tillmann Klamroth  2004   *
 ********************************************************************************/ 
#include <complex>
#include <math.h>
#include <fstream>

using namespace std;
#define Complex complex<double>

//Functions
void calc_pol(int nrolp, double*  Ao_p, double* phase_p, int* pol_type, double* pol_vec);
inline double vec_norm(double* v3);
double vec_pord(double* v1, double* v2, double* v3);
double sca_pord(double* v1, double* v2);
double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time);
void wirte_las(char* dumc, int plas, double omega, double tpuls, double width, double* Ao, double* phase, 
               double* pol_vec, int pol_type);

//Extern Functions

/*******************************************************************************
 * Analyze laser pulse params                                                  *
 *                                                                             *
 ******************************************************************************/

void calc_pol(int nrolp, double*  Ao_p, double* phase_p, int* pol_type, double* pol_vec){
  double v1[3], v2[3]; 
  for(int x = 0; x < nrolp; x++){

    //Build Polarisation at Peak_time (v1) and omega*\delta t = PI/2 later (v2)
    for(int y = 0; y < 3; y++){
      v1[y] = Ao_p[3*x+y]*cos(phase_p[3*x+y]);
      v2[y] = Ao_p[3*x+y]*cos(M_PI/2.+phase_p[3*x+y]);
    }
    

    //check if v1 is collinear with v2
    double norm_vecpr = vec_pord(v1,  v2, &(pol_vec[3*x]));

    //If yes -> linear polarized light, store polarization in pol_vec
    if(norm_vecpr/(vec_norm(v1)*vec_norm(v2)) < 1e-8 || (vec_norm(v1)*vec_norm(v2)) < 1e-24){
      pol_type[x] = 0;
      double normc = vec_norm(&(Ao_p[3*x]));
      if(normc == 0) normc = 1.;
      for(int y = 0; y < 3; y++) pol_vec[3*x+y] = Ao_p[3*x+y]/normc;
    }
    //if no -> eliptical polarized light, store pointing vector in pol_vec
    else{
      pol_type[x] = 1;
      for(int y = 0; y < 3; y++) pol_vec[3*x+y] = pol_vec[3*x+y]/norm_vecpr;
    }
  }
}

/*******************************************************************************
 * calc field in one dimension (x,y, or z)                                     *
 *                                                                             *
 ******************************************************************************/

double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time){
  if(fabs(curr_time-tpeak) > width) return(0.);
  else{
    return(Ao*pow(cos(M_PI/(2.*width)*(curr_time-tpeak)),2)*cos(omega*(curr_time-tpeak)+phase));
  }
}

/*******************************************************************************
 * write laser puls "human readable" to file                                   *
 *                                                                             *
 ******************************************************************************/

void wirte_las(char* dumc, int plas, double omega, double tpuls, double width, double* Ao, double* phase, 
               double* pol_vec, int pol_type){
  ofstream outf(dumc);

  double curr_t = tpuls-width;
  double dt = 2.*width/(plas-1.);
  
  if(pol_type == 0){
    for(int x = 0; x < plas; x++){
      double fv[3];
            
      for(int y = 0; y < 3; y++)
        fv[y] = calc_field( omega,  tpuls,  width, Ao[y],  phase[y], curr_t);
      
      outf << curr_t << " 0. " << sca_pord(fv,pol_vec) << "\n";
      curr_t += dt;
    }
  }else{
    double fv[3], v1[3], v2[3];
    for(int y = 0; y < 3; y++){
      v1[y] = Ao[y]*cos(phase[y]);
      v2[y] = Ao[y]*cos(M_PI/2.+phase[y]);
    }
    double nv1 = vec_norm(v1);
    double nv2 = vec_norm(v2);

    for(int y = 0; y < 3; y++){
      v1[y] /= nv1;
      v2[y] /= nv2;
    }
   
    for(int x = 0; x < plas; x++){
    
      for(int y = 0; y < 3; y++)
        fv[y] = calc_field( omega,  tpuls,  width, Ao[y],  phase[y], curr_t);
      outf <<  curr_t << " " << sca_pord(v1,fv) << " " << sca_pord(v2,fv) << "\n";
      curr_t += dt;
    }
  }
  
  outf.close();
}


/*******************************************************************************
 * simple 3d vector routines                                                   *
 *                                                                             *
 ******************************************************************************/

//calc vecnorm
inline double vec_norm(double* v3){
  return(sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]));
}


//clacs v1 x v2 = v3  outer produkt
double vec_pord(double* v1, double* v2, double* v3){
  v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
  
  return(vec_norm(v3));
}

//calc inner product 
double sca_pord(double* v1, double* v2){
  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}


