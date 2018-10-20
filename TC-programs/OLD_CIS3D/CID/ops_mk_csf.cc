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

using namespace std;

void mk_csf_TI(int nroe, int* dets, double* facs, int* nr);                                 //Build HF GS
void mk_csf_TII(int nroe, int* dets, double* facs, int* nr, int a, int r);                  //Build single ex singlet CSFs 2. (a,r)
void mk_csf_TIII(int nroe, int* dets, double* facs, int* nr, int a, int r);                 //Built doubly CSFs 3. (aa,rr)
void mk_csf_TIV(int nroe, int* dets, double* facs, int* nr, int a, int r, int s);           //Built doubly CSFs 4. (aa,rs)
void mk_csf_TVa(int nroe, int* dets, double* facs, int* nr, int a, int b, int r);           //Built doubly CSFs 5. (ab,rr)
void mk_csf_TV(int nroe, int* dets, double* facs, int* nr, int a, int b, int r);            //Built doubly CSFs 5. (ab,rr)
void mk_csf_TVI(int nroe, int* dets, double* facs, int* nr, int a, int b, int r, int s);    //Built doubly CSFs 6. A(ab,rs)
void mk_csf_TVII(int nroe, int* dets, double* facs, int* nr, int a, int b, int r, int s);   //Built doubly CSFs 7. B(ab,rs)

void calc_det(int i, int nroe, int llim, int ulim, int* dets, double* facs, int* nr);
extern void print_dets_ket(int nroe, int* det2);


/********************************************************************************
 *                                                                              *
 * Build HF ground state                                                        *
 *                                                                              *
 ********************************************************************************/

void mk_csf_TI(int nroe, int* dets, double* facs, int* nr){
  *nr = 1;
  facs[0] = 1.;
  for(int x = 0; x < nroe; x++) dets[x] = x;
}

/********************************************************************************
 *                                                                              *
 * Build singly excited singlet CFS                                             *
 *                                                                              *
 ********************************************************************************/

void mk_csf_TII(int nroe, int* dets, double* facs, int* nr, int a, int r){
  *nr = 2;
  facs[0] = 1./sqrt(2.);
  facs[1] = 1./sqrt(2.);
  for(int x = 0; x < nroe; x++) dets[x]      = x;
  for(int x = 0; x < nroe; x++) dets[x+nroe] = x;

  dets[2*a       ] = 2*r;                      //psi_a^r alpha
  dets[2*a+1+nroe] = 2*r+1;                    //psi_a^r beta
}


/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TIII  (aa,rr)                               *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TIII(int nroe, int* dets, double* facs, int* nr, int a, int r){
  *nr = 1;
  facs[0] = 1.;
  for(int x = 0; x < nroe; x++) dets[x]      = x;
  dets[2*a     ] = 2*r;                      //psi_aa^rr alpha
  dets[2*a+1   ] = 2*r+1;                    //psi_aa^rr beta
}
  
/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TIV (aa,rs)                                 *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TIV(int nroe, int* dets, double* facs, int* nr, int a, int r, int s){
  *nr = 2;
  facs[0] = 1./sqrt(2.);
  facs[1] = 1./sqrt(2.);
  for(int x = 0; x < nroe; x++) dets[x]        = x;
  for(int x = 0; x < nroe; x++) dets[x+nroe]   = x;
  dets[2*a       ] = 2*r;                    // alpha a -> r
  dets[2*a+1     ] = 2*s+1;                  // beta  a -> s
  dets[2*a +nroe ] = 2*s;                    // alpha a -> s
  dets[2*a+1+nroe] = 2*r+1;                  // beta  a -> r
}

/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TV (ab,rr)                                  *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TVa(int nroe, int* dets, double* facs, int* nr, int a, int b, int r){
  *nr = 2;
  facs[0] = 1./sqrt(2.);
  facs[1] = 1./sqrt(2.);
  for(int x = 0; x < nroe; x++) dets[x]        = x;
  for(int x = 0; x < nroe; x++) dets[x+nroe]   = x;
  dets[2*a       ] = 2*r;                    // alpha: a -> r
  dets[2*b+1     ] = 2*r+1;                  // beta: b -> r
  dets[2*a+1+nroe] = 2*r+1;                  // beta: a -> r
  dets[2*b+nroe  ] = 2*r;                    // alpha: b -> r
}

/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TVb (ab,rr)                                 *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TV(int nroe, int* dets, double* facs, int* nr, int a, int b, int r){
  *nr = 2;
  facs[0] = 1./sqrt(2.);
  facs[1] = 1./sqrt(2.);
  for(int x = 0; x < nroe; x++) dets[x]        = x;
  for(int x = 0; x < nroe; x++) dets[x+nroe]   = x;
  dets[2*a+1     ] = 2*r+1;                  // beta: a -> r   
  dets[2*b       ] = 2*r;                    // alpha: b -> r
  dets[2*a+nroe  ] = 2*r;                    // alpha: a -> r
  dets[2*b+1+nroe] = 2*r+1;                  // beta: b -> r
}

/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TVI  A(ab,rs)                               *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TVI(int nroe, int* dets, double* facs, int* nr, int a, int b, int r, int s){
  static double wurz = sqrt(12.);
  *nr = 6;
  facs[0] =  2./wurz;
  facs[1] =  2./wurz;
  facs[2] = -1./wurz;
  facs[3] =  1./wurz;
  facs[4] =  1./wurz;
  facs[5] = -1./wurz;
  for(int y = 0; y < *nr; y++){
    for(int x = 0; x < nroe; x++) dets[x+nroe*y]   = x;
    dets[2*a         ] = 2*r;                   // alpha: a -> r    
    dets[2*b         ] = 2*s;                   // alpha: b -> s

    dets[2*a+1+nroe  ] = 2*r+1;                 // beta: a -> r
    dets[2*b+1+nroe  ] = 2*s+1;                 // beta: b -> s

    dets[2*a+1+nroe*2] = 2*s+1;                 // beta  a -> s 
    dets[2*b+nroe*2  ] = 2*r;                   // alpha: b -> r

    dets[2*a+1+nroe*3] = 2*r+1;                 // beta: a -> r    
    dets[2*b+nroe*3  ] = 2*s;                   // alpha: b -> s

    dets[2*a+nroe*4  ] = 2*r;                   // alpha: a -> r    
    dets[2*b+1+nroe*4] = 2*s+1;                 // beta: b -> s

    dets[2*a+nroe*5  ] = 2*s;                   // alpha: a -> s    
    dets[2*b+1+nroe*5] = 2*r+1;                 // beta: b -> r
  }
}

/********************************************************************************
 *                                                                              *
 * Build doubly excited singlet CSF TVII B(ab,rs)                               *
 *                                                                              *
 ********************************************************************************/
void mk_csf_TVII(int nroe, int* dets, double* facs, int* nr, int a, int b, int r, int s){
  static double halb = 1./2.;
  *nr = 4;
  facs[0] =  halb;
  facs[1] =  halb;
  facs[2] =  halb;
  facs[3] =  halb;
  for(int y = 0; y < *nr; y++){
  for(int x = 0; x < nroe; x++) dets[x+nroe*y]   = x;
    dets[2*a+1       ] = 2*s+1;                 // beta  a -> s 
    dets[2*b         ] = 2*r;                   // alpha: b -> r

    dets[2*a+1+nroe  ] = 2*r+1;                 // beta: a -> r    
    dets[2*b+nroe    ] = 2*s;                   // alpha: b -> s

    dets[2*a+nroe*2  ] = 2*r;                   // alpha: a -> r    
    dets[2*b+1+nroe*2] = 2*s+1;                 // beta: b -> s

    dets[2*a+nroe*3  ] = 2*s;                   // alpha: a -> s    
    dets[2*b+1+nroe*3] = 2*r+1;                 // beta: b -> r

  }
}

/********************************************************************************
 *                                                                              *
 *                                *
 *                                                                              *
 ********************************************************************************/

void calc_det(int i, int nroe, int llim, int ulim, int* dets, double* facs, int* nr){
  //  int typ = 0;
  static int omo      = nroe/2 - llim;
  static int umo      = ulim - nroe/2+1;
  
  
  // Loops to determine cisd_size
  static int sum1 = 0;
  static int sum2 = 0;

  if(sum1 == 0){
    for (int virt = 1; virt <= umo-1; virt++){
      sum1 += virt; 
    }
    for (int occ = 1; occ <= omo-1; occ++){
      sum2 += occ; 
    }
  }
  
  static int g1 = 0;
  static int g2 = 0;  //omo*umo;
  static int g3 = 1*omo*umo;
  static int g4 = 1*omo*umo + omo*sum1;
  static int g5 = 1*omo*umo + omo*sum1 + sum2*umo;
  static int g6 = 1*omo*umo + omo*sum1 + sum2*umo + sum1*sum2;
  static int g7 = 1*omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;

  //  static int cisd_size = 1 + 2*omo*umo + omo*sum1 + sum2*umo + 2*sum1*sum2;

  //=AB,RS====================================
  if(i > g6 && i <= g7){
    //    clog << i << " -> B-Psi_ab^rs "; // 7. ^1B(ab,rs)
    //     typ = cisd_size - sum1*sum2; 

    int newi = i - g6 - 1;

    int l = newi%(sum1);
    int sub_rs = ulim - nroe/2;   
    int count_rs = 0;

    int k = newi/(sum1);
    int sub_ab = nroe/2-1-llim;
    int count_ab = 0;

    while(k >= sub_ab){
      k = k - sub_ab;
      sub_ab = sub_ab - 1;
      count_ab++;
    } 
 
    while(l >= sub_rs){
      l = l - sub_rs;
      sub_rs = sub_rs - 1;
      count_rs++;
    }
    
    int a = count_ab + llim;
    int b = k + a + 1;
    int r = count_rs + nroe/2;
    int s = l + r + 1; 

    mk_csf_TVII(nroe, dets, facs, nr, a, b, r, s);
  }

  //=AB,RS====================================
  if(i > g5 && i <= g6){
    //    clog << i << " -> A-Psi_ab^rs "; // 6. ^1A(ab,rs)
    int newi = i - g5 - 1;

    int l = newi%(sum1);
    int sub_rs = ulim - nroe/2;   
    int count_rs = 0;

    int k = newi/(sum1);
    int sub_ab = nroe/2-1-llim;
    int count_ab = 0;

    while(k >= sub_ab){
      k = k - sub_ab;
      sub_ab = sub_ab - 1;
      count_ab++;
    } 
 
    while(l >= sub_rs){
      l = l - sub_rs;
      sub_rs = sub_rs - 1;
      count_rs++;
    }
    
    int a = count_ab + llim;
    int b = k + a + 1;

    int r = count_rs + nroe/2;
    int s = l + r + 1; 

    mk_csf_TVI(nroe, dets, facs, nr, a, b, r, s);
  }
  //=AB,RR====================================
  if(i > g4 && i <= g5){
    //    clog << i << " -> Psi_ab^rr "; // 5. ^1(ab,rr)
    int newi = i - g4 - 1;
    int r = newi%(ulim-nroe/2+1)+nroe/2;       //Rest von (det_nr/umo)+nroe/2
    int k = newi/(ulim-nroe/2+1);              // det_nr/umo
    int sub = nroe/2-1-llim;                   // omo-1
    int count = 0;
    while(k >= sub){
      k = k - sub;
      sub = sub - 1;
      count++;
    }
    int a = count + llim;
    int b = k + a + 1;

    mk_csf_TV(nroe, dets, facs, nr, a, b, r);
  }
  //=AA,RS====================================
  if(i > g3 && i <= g4){
    //    clog << i << " -> Psi_aa^rs "; // 4. ^1(aa,rs)
    int newi = i - g3 - 1;    
    int a = newi%(nroe/2-llim)+llim;  
    int l = newi/(nroe/2 - llim); 
    int sub = ulim - nroe/2;   
    int count = 0;
    while(l >= sub){
      l = l - sub;
      sub = sub - 1;
      count++;
    }
    int r = count + nroe/2;
    int s = l + r + 1;    

    mk_csf_TIV(nroe, dets, facs, nr, a, r, s);   
  }
  //=AA,RR====================================
  if(i > g2 && i <=  g3){
    //    clog << i << " -> Psi_aa^rr "; // 3. ^1(aa,rr)
    int newi = i - g2 - 1;
    int r = newi%(ulim-nroe/2+1)+nroe/2;
    int a = newi/(ulim-nroe/2+1)+llim; 

    mk_csf_TIII(nroe, dets, facs, nr, a, r);
  }

  //=A,R====================================
  if(i > g1 && i <= g2){
    //    clog << i << " -> Psi_a^r "; // 2. ^1(a,r)
    int newi = i - g1 - 1;
    int r = newi%(ulim-nroe/2+1)+nroe/2;
    int a = newi/(ulim-nroe/2+1)+llim; 

    mk_csf_TII(nroe, dets, facs, nr, a, r);
  }

  //=0====================================
  if(i <= g1){
    //    clog << i << " -> Psi_0 "; // 1. Psi_0   

    mk_csf_TI(nroe, dets, facs, nr);
    //    typ = 0;
  }
} 


