/********************************************************************************
 * CIS3(D) Program suite                                                        *
 *                                                                              *
 * file: mk_field.cc                                                            *
 *                                                                              *
 * computes initial field for OCT                                               *
 *                                                                              *
 *        1.    laser pulses                                                    *
 *        2.    random field frequency restrictions                            *
 *                                                                              *
 *                                                    Tillmann Klamroth  2005   *
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
#include <fftw.h>
#include <complex>
#include <stdlib.h>

using namespace std;
#define Complex complex<double>

double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time);
int    rem_com(char* filename, char* streamstring, int string_length);
void   status(ofstream* outf);

int main(int argc, char* argv[]){
  
  //INITIALIZE RANDOM NUMBER GENERATOR
  time_t start_time;
  time(&start_time);
  srandom((int)start_time);

  //INPUT DATA 
  int    nrots;                //nr of time steps
  double dt;                   //time step length
  int    fac[3];               //Polarization vector
  double sto;                  //mid time of shape function 
  double sw;                   //width of shape functions
  int    expo;                 //exponent for shape function 
  int    mode = 2;             //0 Laser pulses, 1 random field


  if(argc != 3 && argc != 2){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }

  if(argc == 2){
    clog << "Reading field and printing to std-out (t,ex,ey,ez,shape)\n";
    ifstream inf(argv[1]);
    inf.read((char *) &nrots, sizeof(int));
    inf.read((char *) &dt, sizeof(double));
    inf.read((char *) fac, 3*sizeof(int));
    double* data = new double[4*nrots];
    inf.read((char *) data, 4*nrots*sizeof(double));
    for(int x = 0; x < nrots; x++)
      cout << x*dt 
	   << " " << data[0*nrots+x] 
	   << " " << data[1*nrots+x] 
	   << " " << data[2*nrots+x] 
	   << " " << data[3*nrots+x] << "\n";
    return(0);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "MK_FIELD [CIS3(D) Suite]\n";
  outf << "Execution started on/at:\n";
  status(&outf);

  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";
 
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

  
  ist >> nrots >> dt >> fac[0] >> fac[1] >> fac[2] >> sto >> sw >> expo >> mode;
  

  outf << "Nr of time steps is " << nrots << "\n";   
  outf << "dt is " << dt << "\n";	
  outf << "Polarization vector: " << fac[0] << " " << fac[1] << " " << fac[2] << "\n";
  if((fac[0] != 0 && fac[0] != 1) || (fac[1] != 0 && fac[1] != 1) || (fac[2] != 0 && fac[2] != 1)){
    outf << "Only 1. or 0. allowed in polarization vector\n";
    outf.flush();
    exit(1);
  }
  outf << "Shape data:\n";
  outf << "mid time " << sto  << "\n"; 
  outf << "width "    << sw   << "\n";	 
  outf << "exponent " << expo << "\n";

  //Memory allocation
  double* field = new double[nrots*3];

  for(int x = 0; x < 3*nrots; x++) field[x] = 0.;

  double* shape = new double[nrots];
  for(int x = 0; x < nrots; x++){
    double logshape = pow(fabs(x*dt-sto)/(0.5*sw),expo);
    shape[x] = 0.;
    if(logshape < 1.e3)
      shape[x] = exp(-logshape);
  }

  outf << "-------------------------------------------------------------------------------\n";
  int nrolp;
  if(mode == 0){
    outf << "Creating inital pulse from Laser pulses\n";
    ist >> nrolp;
    outf << nrolp << " sets of laser pulse parameters will be read\n\n";
    outf << "Reading laserpulses\n";
    
    //laser 
    double     *omega_p;                //laser frequencies                        nrolp
    double     *t_puls;                 //laser peak times                         nrolp
    double     *width_p;                //laser pulse width                        nrolp
    double     *Ao_p;                   //Amplidute vector                       3*nrolp
    double     *phase_p;                //Phase shifts in x y z                  3*nrolp
    double     *pol_vec;                //Polarization vectros                   3*nrolp  
      
    int  doubmem = 12*nrolp;
    double* dumd = new double[doubmem]; int inc = 0;
    
    //laser
    omega_p = &(dumd[inc]); inc += nrolp; t_puls  = &(dumd[inc]); inc += nrolp; width_p = &(dumd[inc]); inc += nrolp;
    Ao_p    = &(dumd[inc]); inc += 3*nrolp; phase_p = &(dumd[inc]); inc += 3*nrolp; pol_vec = &(dumd[inc]); inc += 3*nrolp;
    
    int*    pol_type = new int[nrolp];
    
    for(int x = 0; x < nrolp; x++) {
      ist >> omega_p[x] >> t_puls[x] >> width_p[x]
	  >> Ao_p[3*x+0] >> phase_p[3*x+0] >> Ao_p[3*x+1] >> phase_p[3*x+1] >> Ao_p[3*x+2] >> phase_p[3*x+2];
    }
    outf << "Nr   Omega   t_peak     width    X (Ao,phase),  Y (Ao,phase),  Z (Ao,phase)\n";
    outf << "..............................................................................\n";
    for(int x = 0; x < nrolp; x++){
      sprintf(dumc,"%i %2.4f %9.2f  %9.2f %2.4f %+2.4f, %2.4f %+2.4f, %2.4f %+2.4f\n",x,  
	      omega_p[x],t_puls[x],width_p[x],Ao_p[3*x+0],phase_p[3*x+0],Ao_p[3*x+1],phase_p[3*x+1],Ao_p[3*x+2],phase_p[3*x+2]);
      outf << dumc;
    }
    

    for(int x = 0; x < nrots; x++){
      for(int l = 0; l < nrolp; l++){
	for(int y = 0; y < 3; y++){
	  field[y*nrots+x] += calc_field( omega_p[l], t_puls[l], width_p[l], Ao_p[3*l+y],  phase_p[3*l+y], x*dt);
	}
      }
    }  
  }
  
  if(mode == 1){
    outf << "Creating random inital field\n";
    
    double fmax; 
    double fmin; 
    double emax; 
    ist >> fmin >> fmax >> emax;
    
    outf << "fmin is " << fmin << "\n";		   
    outf << "fmax is " << fmax << "\n";		   
    outf << "emax is " << emax << "\n";      
    
    Complex* cfield = new Complex[nrots];
    Complex* cfreq  = new Complex[nrots];
    
    fftwnd_plan planf;
    fftwnd_plan planb;

    planf = fftwnd_create_plan(1, &nrots, FFTW_FORWARD,  FFTW_ESTIMATE);
    planb = fftwnd_create_plan(1, &nrots, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int y = 0; y < 3; y++){  
      for(int x = 0; x < nrots; x++){
	cfield[x] = ((double)rand()/(double) RAND_MAX -.5)*2;
      }
      
      fftwnd(planf,1,(FFTW_COMPLEX*)cfield,1,0,(FFTW_COMPLEX*)cfreq,1,0);
	
      for(int x = 0; x < nrots;x++){
	int kpunkt  = x; 
	if(kpunkt >= nrots/2) kpunkt = -(kpunkt - nrots); 
	double om = kpunkt*2*M_PI/dt/(double)nrots;
	if(om < fmin || om > fmax)
	  cfreq[x] = 0.;
      }

      fftwnd(planb,1,(FFTW_COMPLEX*)cfreq,1,0,(FFTW_COMPLEX*)cfield,1,0);

      double maxf = 0.;
  
      for(int x = 0; x < nrots;x++){
	if(maxf < fabs(real(cfield[x]))) maxf = fabs(real(cfield[x]));
      }
  
  
    
      for(int x = 0; x < nrots; x++)
	field[y*nrots+x] = real(cfield[x])/maxf*emax;
    }
    
  }   
  
  for(int y = 0; y < 3; y++){
    for(int x = 0; x < nrots; x++){
      field[y*nrots+x] *= (double) fac[y];
    }
  }
  sprintf(dumc,"%s.ef",argv[2]);
  ofstream datf(dumc);
  outf << "Writing field to " << dumc << "\n";
  datf.write((char *) &nrots, sizeof(int));
  datf.write((char *) &dt,    sizeof(double));
  datf.write((char *) fac,    3*sizeof(int));
  datf.write((char *) field,  3*nrots*sizeof(double));
  datf.write((char *) shape,  nrots*sizeof(double));
}

double calc_field(double omega, double tpeak, double width, double Ao, double phase, double curr_time){
  if(fabs(curr_time-tpeak) > width) return(0.);
  else{
    return(Ao*pow(cos(M_PI/(2.*width)*(curr_time-tpeak)),2)*cos(omega*(curr_time-tpeak)+phase));
  }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      STATUS                                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void status(ofstream* outf){
  char str1[1024];
  
  size_t len = 1024;
  time_t curr_time;
  time(&curr_time);
  
  gethostname(str1, len);
  
  *outf << "Host: " << str1 << ", Date: " << ctime(&curr_time);
  outf->flush();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      REM COM                                                  */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int rem_com(char* filename, char* streamstring, int string_length){
  const char com_B = '#';
  const char com_E = '\n';
  
  int pos = 0;
  char cc;

  ifstream inf(filename);
  
  while(inf.get(cc)&& pos < string_length-1){
    if(cc != com_B) 
      streamstring[pos++] = cc;
    else{
      while(cc != com_E && inf.get(cc));
      streamstring[pos++] = com_E;
    }
  }
  streamstring[pos] = 0;
  if(pos == string_length-1){
    cerr << "Buffer size exceeded !\n"; exit(0);
  }
  return(strlen(streamstring));
}
