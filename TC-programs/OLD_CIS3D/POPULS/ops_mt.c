#include <pthread.h>
#include <math.h>
#include <fftw.h>
#include <stdio.h>

//C -> Very ugly code

//DEFINITION OF COMPLEX NUMBERS
typedef struct {
  double re, im;
} mycom;

//DEFINITION OF DATA ARRAY NEEDED FOR THE TRANSFORMATION (small space)
typedef struct {
  double* m;
  double* vi;
  double* vo;
  int*    nrop;
} transdata;


//DEFINITION OF DATA ARRAY NEEDED FOR THE TRANSFORMATION (large space)
typedef struct {
  double* m;
  double* vi;
  double* vo;
  long long int*    nrop;
} ltransdata;


//THREAD FUNCTION FOR FORWARD TRANSFORMATION (small space)
void transforeT(void* td){
  transdata* Td = (transdata *)td;
  double* Mat       = Td->m;
  double* Vec_in    = Td->vi;
  double* Vec_out   = Td->vo;
  int*       Nrop   = Td->nrop;
  int x, y;
  
  for(x = 0; x < 2* *Nrop; x++){
    Vec_out[x] = 0.;
  }
  
 
  for(x = 0; x < *Nrop; x++){
    for(y = 0; y < *Nrop; y++){
      Vec_out[2*x]    +=     Mat[(x* (*Nrop))+y]*Vec_in[2*y];
      Vec_out[2*x+1]  +=     Mat[(x* (*Nrop))+y]*Vec_in[2*y+1];
    }
  }
  
}

//THREAD FUNCTION FOR FORWARD TRANSFORMATION (large space)
void ltransforeT(void* td){
  ltransdata* Td = (ltransdata *)td;
  double* Mat       = Td->m;
  double* Vec_in    = Td->vi;
  double* Vec_out   = Td->vo;
  long long int*       Nrop   = Td->nrop;
  long long int x, y;
  
  for(x = 0; x < 2* *Nrop; x++){
    Vec_out[x] = 0.;
  }
  
 
  for(x = 0; x < *Nrop; x++){
    for(y = 0; y < *Nrop; y++){
      Vec_out[2*x]    +=     Mat[(x* (*Nrop))+y]*Vec_in[2*y];
      Vec_out[2*x+1]  +=     Mat[(x* (*Nrop))+y]*Vec_in[2*y+1];
    }
  }
  
}


//THREAD FUNCTION FOR BACKWARD TRANSFORMATION (small space)
void transbackT(void* td){
  transdata* Td = (transdata *)td;
  double* Mat       = Td->m;
  double* Vec_in    = Td->vi;
  double* Vec_out   = Td->vo;
  int*       Nrop   = Td->nrop;
  int x, y;
  
  for(x = 0; x < 2* *Nrop; x++){
    Vec_out[x] = 0.;
  }
  
  for(x = 0; x < *Nrop; x++){
    for(y = 0; y < *Nrop; y++){
      Vec_out[2*y]    +=     Mat[x* *Nrop+y]*Vec_in[2*x];
      Vec_out[2*y+1]  +=     Mat[x* *Nrop+y]*Vec_in[2*x+1];
    }
  }
  
}

//THREAD FUNCTION FOR BACKWARD TRANSFORMATION (large space)
void ltransbackT(void* td){
  ltransdata* Td = (ltransdata *)td;
  double* Mat       = Td->m;
  double* Vec_in    = Td->vi;
  double* Vec_out   = Td->vo;
  long long int*       Nrop   = Td->nrop;
  long long int x, y;
  
  for(x = 0; x < 2* *Nrop; x++){
    Vec_out[x] = 0.;
  }
  
  for(x = 0; x < *Nrop; x++){
    for(y = 0; y < *Nrop; y++){
      Vec_out[2*y]    +=     Mat[x* *Nrop+y]*Vec_in[2*x];
      Vec_out[2*y+1]  +=     Mat[x* *Nrop+y]*Vec_in[2*x+1];
    }
  }
  
}

//MAKE TRANSFORMATION MULTI THREADED (small space) DIR=1 is forward !!
void mk_thread_trans(int dir, double* mat1, double* mat2, mycom* vec_in1, mycom* vec_in2,
		     mycom* vec_out1, mycom* vec_out2, int* nrop1, int* nrop2){
  pthread_t thread1, thread2;  
  int  iret1, iret2;
  transdata td1, td2;

  td1.m    = mat1;
  td1.vi   = (double* ) vec_in1;
  td1.vo   = (double* ) vec_out1;
  td1.nrop = nrop1;  

  td2.m    = mat2;
  td2.vi   = (double* ) vec_in2;
  td2.vo   = (double* ) vec_out2;
  td2.nrop = nrop2;  

  if(dir == 1){
    iret1 = pthread_create( &thread1, NULL,  (void* ) &transforeT , (void*) &td1);
    iret2 = pthread_create( &thread2, NULL,  (void* ) &transforeT , (void*) &td2);
  }else{
    iret1 = pthread_create( &thread1, NULL,  (void* ) &transbackT , (void*) &td1);
    iret2 = pthread_create( &thread2, NULL,  (void* ) &transbackT , (void*) &td2);
  }
  
  
  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL); 
  
}

//MAKE TRANSFORMATION MULTI THREADED (large space) DIR=1 is forward !!
void lmk_thread_trans(int dir, double* mat1, double* mat2, mycom* vec_in1, mycom* vec_in2,
		      mycom* vec_out1, mycom* vec_out2, long long int* nrop1, long long int* nrop2){
  pthread_t thread1, thread2;  
  int  iret1, iret2;
  ltransdata td1, td2;

  td1.m    = mat1;
  td1.vi   = (double* ) vec_in1;
  td1.vo   = (double* ) vec_out1;
  td1.nrop = nrop1;  

  td2.m    = mat2;
  td2.vi   = (double* ) vec_in2;
  td2.vo   = (double* ) vec_out2;
  td2.nrop = nrop2;  

  if(dir == 1){
    iret1 = pthread_create( &thread1, NULL,  (void* ) &ltransforeT , (void*) &td1);
    iret2 = pthread_create( &thread2, NULL,  (void* ) &ltransforeT , (void*) &td2);
  }else{
    iret1 = pthread_create( &thread1, NULL,  (void* ) &ltransbackT , (void*) &td1);
    iret2 = pthread_create( &thread2, NULL,  (void* ) &ltransbackT , (void*) &td2);
  }
  
  
  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL); 
  
}
