#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

//Functions
int calc_diffs_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int nroao, int nroe);
void show_diff_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int p, int q, int nroao, int nroe);
void show_diffs_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int p, int q, int r, int s, int nroao, int nroe);
void produce_pops(int* pops, int i, int j, int a , int b, int nroao, int nroe);

int calc_diffs_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int nroao, int nroe){
  //first sort the values
  int* I = new int[nroao];
  int* J = new int[nroao];
  for(int x = 0; x < (nroe/2); x++){
    if(x != Ii && x != i){
      I[x] = 1;
    }else{
      I[x] = 0;
    }
    if(x != Ji && x != j){
      J[x] = 1;
    }else{
      J[x] = 0;
    }
  }
  for(int x = (nroe/2); x < nroao; x++){
    if(x == If || x == a){
      I[x] = 1;
    }else{
      I[x] = 0;
    }
    if(x == Jf || x == b){
      J[x] = 1;
    }else{
      J[x] = 0;
    }
  }
  int diffcount = 0;
  for(int x = 0; x < nroao; x++){
    if(fabs(I[x]-J[x]) == 1){
      diffcount++;
    }else{
      continue;
    }
  }
    return(diffcount);
}
  
void show_diff_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int p, int q, int nroao, int nroe){
   int* I = new int[nroao];
   int* J = new int[nroao];
   for(int x = 0; x < (nroe/2); x++){
     if(x != Ii && x != i){
       I[x] = 1;
     }else{
       I[x] = 0;
     }
     if(x != Ji && x != j){
       J[x] = 1;
     }else{
       J[x] = 0;
     }
   }
   for(int x = (nroe/2); x < nroao; x++){
     if(x == If || x == a){
       I[x] = 1;
     }else{
       I[x] = 0;
     }
     if(x == Jf || x == b){
       J[x] = 1;
     }else{
       J[x] = 0;
     }
   }
   int diffcount = 0;
   for(int x = 0; x < nroao; x++){
     if(fabs(I[x]-J[x]) == 1 && diffcount == 0){
       p = x;
       diffcount++;
     }else if(fabs(I[x]-J[x]) == 1 && diffcount == 1){
       q = x;
       diffcount++;
       break;
     }else{
       continue;
     }
   }
}

void show_diffs_dd(int Ii, int Ji, int i, int j, int If, int Jf, int a, int b, int p, int q, int r, int s, int nroao, int nroe){
  int* I = new int[nroao];
  int* J = new int[nroao];
  for(int x = 0; x < (nroe/2); x++){
    if(x != Ii && x != i){
      I[x] = 1;
    }else{
      I[x] = 0;
    }
    if(x != Ji && x != j){
      J[x] = 1;
    }else{
      J[x] = 0;
    }
  }
  for(int x = (nroe/2); x < nroao; x++){
    if(x == If || x == a){
      I[x] = 1;
    }else{
      I[x] = 0;
    }
    if(x == Jf || x == b){
      J[x] = 1;
    }else{
      J[x] = 0;
    }
  }
  int diffcount = 0;
  for(int x = 0; x < nroao; x++){
    if(fabs(I[x]-J[x]) == 1 && diffcount == 0){
      p = x;
      diffcount++;
    }else if(fabs(I[x]-J[x]) == 1 && diffcount == 1){
      q = x;
      diffcount++;
    }else if(fabs(I[x]-J[x]) == 1 && diffcount == 2){
      r = x;
      diffcount++;
    }else if(fabs(I[x]-J[x]) == 1 && diffcount == 3){
      s = x;
      diffcount++;
      break;
    }else{
      continue;
    }
  }
}

void produce_pops(int* pops, int i, int j, int a , int b, int nroao, int nroe){
  for(int x = 0; x < (nroe/2); x++){
    if(x == i || x == j){
      pops[x] = 0;
    }else{
      pops[x] = 1;
    }
  }
  for(int x = (nroe/2); x < nroao; x++){
    if(x == a || x == b){
      pops[x] = 1;
    }else{
      pops[x] = 0;
    }
  }
}
