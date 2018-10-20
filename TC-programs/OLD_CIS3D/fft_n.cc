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
  if(argc != 7 && argc != 8){
    cerr << "need index, #points, #timesteps, #columns, time and signal column, (input file)\n";
    exit(1);
  }

  //INPUT DATA
  int    index = atoi(argv[1]);      // if index = 0: input from terminal; if index = 1: input from file
  int    tpoints = atoi(argv[2]);          //nr of points to use for fourier transfromation should be minimal of nrots/2.
  int    nrots = atoi(argv[3]);            //nrots in input file
  int    nr_col = atoi(argv[4]);           //nr of columns in input file
  int    t_col = atoi(argv[5]);            //colum with t
  int    sig_col = atoi(argv[6]);          //colum with signal
  
  double* data     = new double[nrots];
  double* tvals    = new double[nrots];
  double* inp_line = new double[nr_col];

  double dt;

  if(index == 1){
   if(argc == 7){
    cerr << "need input file\n";
    exit(1);
   }
  ifstream inf;
  inf.open(argv[7]);
  for(int x = 0; x < nrots; x++){
    for(int y = 0; y < nr_col; y++)
      inf >> inp_line[y];
    tvals[x] = inp_line[t_col-1];
    data[x]  = inp_line[sig_col-1];
  }
  inf.close();
  }else{
  for(int x = 0; x < nrots; x++){
   for(int y = 0; y < nr_col; y++)
    cin >> inp_line[y];
    tvals[x] = inp_line[t_col-1];
    data[x] = inp_line[sig_col-1];
   }
  }
  dt = (tvals[nrots-1]-tvals[0])/(nrots-1.);

  Complex*  fftw_vec_in  =  new Complex[tpoints];
  Complex*  fftw_vec_out =  new Complex[tpoints];

  fftwnd_plan plan_et, plan_bt;
  plan_et   = fftwnd_create_plan(1, &tpoints, FFTW_BACKWARD,   FFTW_ESTIMATE);
  plan_bt   = fftwnd_create_plan(1, &tpoints, FFTW_FORWARD ,  FFTW_ESTIMATE);

  for(int x = 0; x < tpoints; x++) fftw_vec_in[x] = 0.;
  for(int x = 0; x < nrots; x++) fftw_vec_in[x] = data[x];
  
  fftwnd(plan_et, 1,(FFTW_COMPLEX*) fftw_vec_in,1,0,(FFTW_COMPLEX*) fftw_vec_out,1,0);
  for(int x = 0; x < tpoints; x++) fftw_vec_out[x] *= 1./sqrt((double) tpoints);
  
  for(int x = 0; x < tpoints/2; x++)
    cout << 2.*M_PI/(dt*tpoints) *x << " " << norm(fftw_vec_out[x]) << " " 
      //	 <<  norm(fftw_vec_out[x]*exp(-pow(wo-x*2.*M_PI/(dt*tpoints),2)/(2.*sig_w*sig_w)))  
	 <<  "\n";

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



