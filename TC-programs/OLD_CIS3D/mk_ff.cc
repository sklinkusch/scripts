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
  if(argc != 4 && argc != 5){
    cerr << "need index, input-file, signal column, (output-prefix)!\n";
    exit(1);
  }

int nrots;
double dt;
int* pfac = new int[3];

ifstream ief(argv[2]);
ief.read((char *) &nrots, sizeof(int));
ief.read((char *) &dt, sizeof(double));
ief.read((char *) pfac, 3*sizeof(int));

double* field_h = new double[3*nrots];
double* shape = new double[nrots];

ief.read((char *) field_h, 3*nrots*sizeof(double));
ief.read((char *) shape, nrots*sizeof(double));


  //INPUT DATA FOR FOURIER TRANSFORMATION
  int    index = atoi(argv[1]);            // if index = 0: output to terminal; if index = 1: output to file
  int    nr_col = 3;                       //nr of columns in input file
  int    sig_col = atoi(argv[3]);          //column with signal

  int    tpoints = 10*nrots;
  
  double* data     = new double[nrots];
  double* tvals    = new double[nrots];
  double* inp_line = new double[nr_col];

  for(int x = 0; x < nrots; x++){
   for(int y = 0; y < nr_col; y++)
    inp_line[y] = field_h[y*nrots+x];
    tvals[x] = x*dt;
    data[x] = inp_line[sig_col-1];
   }

  

  Complex*  fftw_vec_in  =  new Complex[tpoints];
  Complex*  fftw_vec_out =  new Complex[tpoints];

  fftwnd_plan plan_et, plan_bt;
  plan_et   = fftwnd_create_plan(1, &tpoints, FFTW_BACKWARD,   FFTW_ESTIMATE);
  plan_bt   = fftwnd_create_plan(1, &tpoints, FFTW_FORWARD ,  FFTW_ESTIMATE);

  for(int x = 0; x < tpoints; x++) fftw_vec_in[x] = 0.;
  for(int x = 0; x < nrots; x++) fftw_vec_in[x] = data[x];
  
  fftwnd(plan_et, 1,(FFTW_COMPLEX*) fftw_vec_in,1,0,(FFTW_COMPLEX*) fftw_vec_out,1,0);
  for(int x = 0; x < tpoints; x++) fftw_vec_out[x] *= 1./sqrt((double) tpoints);
  
if(index == 1){
 if(argc == 4){
  cerr << "Need output-prefix!\n";
  exit(1);
  }
  ofstream outf;
  char dumc[2048];
  sprintf(dumc, "%s.ew", argv[4]);
  outf.open(dumc); 
  for(int x = 0; x < nrots; x++){
   outf << tvals[x] << " " << data[x] << "\n";
  }
  outf << "#------------------------------------------------------------------------------------#\n";
  for(int x = 0; x < tpoints/2; x++){
    outf << 2.*M_PI/(dt*tpoints) *x << " " << norm(fftw_vec_out[x]) << "\n"; 
  }
  outf.close();
  }else{
  for(int x = 0; x < tpoints/2; x++){
    cout << 2.*M_PI/(dt*tpoints) *x << " " << norm(fftw_vec_out[x]) << "\n";
  }

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



