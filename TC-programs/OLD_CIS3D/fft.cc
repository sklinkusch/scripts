#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include <fftw.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

#define Complex complex<double>

int rem_com(char* filename, char* streamstring, int string_length);


int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "need input outpref\n";
    exit(1);
  }
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  outf << "Calcuating fourier transformation \n";
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n\n";


  //INPUT DATA
  int    tpoints;          //nr of points to use for fourier transfromation should be minimal of nrots/2.
  int    nrots;            //nrots in input file
  int    nr_col;           //nr of collums in input file
  int    t_col;            //colum with t
  int    sig_col;          //colum with signal
  char   data_file[1024];
  
  //  double wo;
  //  double sig_w;
  
  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);

  ist >> tpoints >> nrots >> nr_col >> t_col >> sig_col >> data_file;

  //  ist >> wo >> sig_w;
 
  outf << "Using " << tpoints << " for fourier interpolation in frequency domain.\n";
  outf << "Nr of time steps: " << nrots << "\n";
  outf << "Nr of columns in input: " << nr_col << "\n";
  outf << "Time is in column: " << t_col << " (columns start with 1 !!!!!)\n";
  //  outf << "Signal is in column: " << sig_col << "(columns start with 1 !!!!!) \n";
  outf << "Data file is: " << data_file << "\n";
  //  outf << "Window omega is: " << wo << "\n";
  //  outf << "Sigma window is: " << sig_w << "\n";

  double* data     = new double[nrots];
  double* tvals    = new double[nrots];
  double* inp_line = new double[nr_col];

  double dt;
  ifstream inf(data_file);
  
  for(int x = 0; x < nrots; x++){
    for(int y = 0; y < nr_col; y++)
      inf >> inp_line[y];
    tvals[x] = inp_line[t_col-1];
    data[x]  = inp_line[sig_col-1];
  }

  dt = (tvals[nrots-1]-tvals[0])/(nrots-1.);
  clog << "dt: " << dt << "\n";

  outf << "Data read:\n";
  outf << "dt: " << dt << "\n";
  
  ofstream datf;
  sprintf(dumc,"%s.aet",argv[2]);
  outf << "-------------------------------------------------------------------------------\n";
  
  Complex*  fftw_vec_in  =  new Complex[tpoints];
  Complex*  fftw_vec_out =  new Complex[tpoints];

  fftwnd_plan plan_et, plan_bt;
  plan_et   = fftwnd_create_plan(1, &tpoints, FFTW_BACKWARD,   FFTW_ESTIMATE);
  plan_bt   = fftwnd_create_plan(1, &tpoints, FFTW_FORWARD ,  FFTW_ESTIMATE);

  for(int x = 0; x < tpoints; x++) fftw_vec_in[x] = 0.;
  for(int x = 0; x < nrots; x++) fftw_vec_in[x] = data[x];
  
  fftwnd(plan_et, 1,(FFTW_COMPLEX*) fftw_vec_in,1,0,(FFTW_COMPLEX*) fftw_vec_out,1,0);
  for(int x = 0; x < tpoints; x++) fftw_vec_out[x] *= 1./sqrt((double) tpoints);
  
  sprintf(dumc,"%s.ew",argv[2]);
  datf.open(dumc);
  outf << "Printing |E(w)|^2 to " << dumc << "\n";
  for(int x = 0; x < tpoints/2; x++)
    datf << 2.*M_PI/(dt*tpoints) *x << " " << norm(fftw_vec_out[x]) << " " 
      //	 <<  norm(fftw_vec_out[x]*exp(-pow(wo-x*2.*M_PI/(dt*tpoints),2)/(2.*sig_w*sig_w)))  
	 <<  "\n";
  datf.close();

//   for(int x = 0; x < tpoints; x++){
//     int X = x; 
//     if(x > tpoints/2) X -= tpoints;
//     fftw_vec_out[x] *= exp(-pow(wo-X*2.*M_PI/(dt*tpoints),2)/(2.*sig_w*sig_w));
//   }
  
//   fftwnd(plan_bt, 1,(FFTW_COMPLEX*) fftw_vec_out,1,0,(FFTW_COMPLEX*) fftw_vec_in,1,0);
  
//   sprintf(dumc,"%s.etw",argv[2]);
//   datf.open(dumc);
//   for(int x = 0; x < nrots; x++){
//     datf << x*dt << " " << norm(fftw_vec_in[x]) << "\n";
//   }
  

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



