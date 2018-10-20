#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>

/******************************************************************************************************************************************************
* GMS_INP                                                                                                                                             *
* Routines to construct the $CONTRL group and the $MASS, $FORCE, $CIS, and $CC groups                                                                 *
*                                                                                                                         Stefan Klinkusch, 2008-2009 *
******************************************************************************************************************************************************/

using namespace std;

//Functions
void contrl_group(int nroa, int nroe, int mult_gs, int charge, ofstream *outf);
void cis_props(ofstream **outf);
void cc_props(ofstream **outf);
void mp2_props(ofstream **outf);
void thermo_props(ofstream **outf);
void iso_props(int nroa, ofstream **outf);
void stat_point(ofstream **outf);
void trudge(ofstream **outf);

//Extern Functions
extern void dft_insert(ofstream **outf);
extern void dft_props(ofstream **outf);
extern void dft_exch(ofstream ***outf);
extern void dft_corr(ofstream ***outf);
extern void dft_comb(ofstream ***outf);
extern void dft_hybr(ofstream ***outf);
extern void dft_meta(ofstream ***outf);
extern void tddft_props(ofstream **outf);

//write $CONTRL group

void contrl_group(int nroa, int nroe, int mult_gs, int charge, ofstream *outf){
 int method_nra = 0;
 int method_nrb = 0;
 int iso_index = 0;
 if(mult_gs == 1){
  while(method_nra < 1 || method_nra > 7){
   cout << "Insert Number of Method to be used: \n";
   cout << "(1) RHF\n";
   cout << "(2) GVB\n";
   cout << "(3) CIS\n";
   cout << "(4) CCSD\n";
   cout << "(5) DFT-Methods\n";
   cout << "(6) TD-DFT\n";
   cout << "(7) MP2\n";
   cin >> method_nra;
  }
  *outf << " $CONTRL ";
  switch(method_nra)
  {
   case 1: *outf << "SCFTYP=RHF ";
           break;
   case 2: *outf << "SCFTYP=GVB ";
           break;
   case 3: *outf << "SCFTYP=RHF CITYP=CIS ";
           break;
   case 4: *outf << "SCFTYP=RHF CCTYP=CCSD ";
           break;
   case 5: *outf << "SCFTYP=RHF DFTTYP=";
           dft_insert(&outf);
           break;
   case 6: *outf << "SCFTYP=RHF TDDFT=EXCITE DFTTYP=";
           dft_insert(&outf);
           break;
   case 7: *outf << "SCFTYP=RHF MPLEVL=2 ";
           break;
  }
 }else{
  while(method_nrb < 1 || method_nrb > 6){
   cout << "Insert Number of Method to be used: \n";
   cout << "(1) UHF\n";
   cout << "(2) ROHF\n";
   cout << "(3) GVB\n";
   cout << "(4) DFT-Methods\n";
   cout << "(5) TD-DFT\n";
   cout << "(6) MP2\n";
   cin >> method_nrb;
  }
  *outf << " $CONTRL ";
  switch(method_nrb)
  {
   case 1: *outf << "SCFTYP=UHF ";
           break;
   case 2: *outf << "SCFTYP=ROHF ";
           break;
   case 3: *outf << "SCFTYP=GVB ";
           break;
   case 4: *outf << "SCFTYP=UHF DFTTYP=";
           dft_insert(&outf);
           break;
   case 5: *outf << "SCFTYP=UHF TDDFT=EXCITE DFTTYP=";
           dft_insert(&outf);
           break;
   case 6: *outf << "SCFTYP=UHF MPLEVL=2 ";
           break;
  }
 }
 int calc_nr = 0;
 while(calc_nr < 1 || calc_nr > 6){
  cout << "Insert Number of Calculation to be used: \n";
  cout << "(1) Single-Point Energy\n";
  cout << "(2) Single-Point Energy + Gradient\n";
  cout << "(3) Frequency Analysis\n";
  cout << "(4) Geometry Optimization\n";
  cout << "(5) non-gradient Total Energy minimization\n";
  cout << "(6) Search for local saddlepoint\n";
  cin >> calc_nr;
 }
 switch(calc_nr)
 {
  case 1: *outf << "RUNTYP=ENERGY ";
          break;
  case 2: *outf << "RUNTYP=GRADIENT ";
          break;
  case 3: *outf << "RUNTYP=HESSIAN ";
          break;
  case 4: *outf << "RUNTYP=OPTIMIZE ";
          break;
  case 5: *outf << "RUNTYP=TRUDGE ";
          break;
  case 6: *outf << "RUNTYP=SADPOINT ";
          break;
 }
 *outf << "UNITS=ANGS \n";
 int gaussian = 0;
 while(gaussian < 1 || gaussian > 2){
  cout << "Want to punch a GAUSSIAN Job file?\n";
  cout << "(1) Yes\n";
  cout << "(2) No\n";
  cin >> gaussian;
 }
 if(gaussian == 1) *outf << "   FRIEND=GAUSSIAN\n";
 if(charge!= 0) *outf << "   ICHARG=" << charge << " MULT=" << mult_gs << " \n";
 *outf << "   COORD=UNIQUE MOLPLT=.T. PLTORB=.T. $END\n";
 switch(method_nra){
  case 3: cis_props(&outf);
          break;
  case 4: cc_props(&outf);
          break;
  case 5: dft_props(&outf);
          break;
  case 6: dft_props(&outf);
          tddft_props(&outf);
          break;
  case 7: mp2_props(&outf);
          break;
 }
 switch(method_nrb){
  case 4: dft_props(&outf);
          break;
  case 5: dft_props(&outf);
          tddft_props(&outf);
          break;
  case 6: mp2_props(&outf);
          break;
 }
 switch(calc_nr){
  case 3: thermo_props(&outf);
          while(iso_index < 1 || iso_index > 2){
           cout << "Other isotopes included?\n";
           cout << "(1) Yes\n";
           cout << "(2) No\n";
           cin >> iso_index;
          }
          if(iso_index == 1) iso_props(nroa, &outf);
          break;
  case 4: stat_point(&outf);
          break;
  case 5: trudge(&outf);
          break;
  case 6: stat_point(&outf);
          break;
 }
} 

//routine to set options for a CIS calculation
//construction of the $CIS group

void cis_props(ofstream **outf){
 int nros = 0;
 while(nros < 1){
  cout << "Insert nr of excited states to be calculated:";
  cin >> nros;
 }
 int state = 0;
 while(state < 1){
  cout << "Insert excited state to be used for geometry optimization and properties:";
  cin >> state;
  if(state < 1 || state > nros){
   state = 0;
  }
 }
 int mult_es = 0;
 while(mult_es != 1 && mult_es !=3){
  cout << "Insert spin multiplicity of excited states:";
  cin >> mult_es;
 }
 **outf << " $CIS NSTATE=" << nros << " IROOT=" << state << " MULT=" << mult_es << " $END\n";
}

//routine to set options for a CCSD calculation
//construction of the $CC group

void cc_props(ofstream **outf){
 int nfc = 0;
 cout << "Insert nr of frozen core orbitals (if negative: default):";
 cin >> nfc;
 unsigned int nfv = 0;
 cout << "Insert nr of frozen virtual orbitals:";
 cin >> nfv;
 **outf << " $CCINP ";
 if(nfc >= 0){
  **outf << "NCORE=" << nfc << " ";
 }
 **outf << "NFZV=" << nfv << " $END\n";
}

//routine to construct the $MP2 group
//settings for a MP2 calculation

void mp2_props(ofstream **outf){
 int nfc = 0;
 cout << "Insert nr of frozen core orbitals (if negative: default):";
 cin >> nfc;
 if(nfc >= 0){
  **outf << " $MP2 NACORE=" << nfc << " NBCORE=" << nfc << " $END\n";
 }
}

//routine to construct the $FORCE group
//thermochemistry for a finite temperature

void thermo_props(ofstream **outf){
 int thermochem = -1;
 while(thermochem < 0 || thermochem > 1){
  cout << "Want the thermochemistry for finite temperatures?\n";
  cout << "(1) Yes\n";
  cout << "(0) No\n";
  cin >> thermochem;
 }

 if(thermochem == 1){
  int nrotmp = 0;
  while(nrotmp < 1 || nrotmp > 10){
   cout << "Insert No of finite temperatures:\n";
   cin >> nrotmp;
  }
 
  **outf << " $FORCE ";
  double* temp = new double[nrotmp];
  for(int xtemp = 0; xtemp < nrotmp; xtemp++){
   temp[xtemp] = 0.;
  }
  **outf << "Temperatures must be positive doubles( > 0., in K, for 0. K, type 0.001)!\n";
  for(int xtemp = 0; xtemp < nrotmp; xtemp++){
   while(temp[xtemp] <= 0.){
    cout << "Insert Temperature No. " << xtemp + 1 << ": ";
    cin >> temp[xtemp];
   }
   **outf << "TEMP(" << xtemp + 1 << ")=" << temp[xtemp] << " ";
  }
  **outf << "$END\n";
 }
}

//routine to construct the $MASS group
//change isotopes for specified atoms

void iso_props(int nroa, ofstream **outf){
 int nrodi = 0;
 while(nrodi < 1 || nrodi > nroa){
  cout << "Insert No of atoms with different isotopes as usual:";
  cin >> nrodi;
 }
 double* atindex = new double[nrodi];
 **outf << " $MASS ";
 double* isotop = new double[nrodi];
 for(int x = 0; x < nrodi; x++){
  cout << "Insert atom No. with different isotop:";
  cin >> atindex[x];
  cout << "Insert atom mass of atom " << atindex[x] << ": "; 
  cin >> isotop[x];
  **outf << "AMASS(" << atindex[x] << ")=" << isotop[x] << " ";
 }
 **outf << "$END\n";
}

//routine to construct the $STATPT group
//settings for geometry optimizations and saddle point searches

void stat_point(ofstream **outf){
 **outf << " $STATPT ";
 int optmet = 0;
 while(optmet < 1 || optmet > 5){
  cout << "Choose optimization method:\n";
  cout << "(1) Quadratic Approximation (QA, default)\n";
  cout << "(2) Straight Newton-Raphson iterate (NR)\n";
  cout << "(3) Rational Function Optimization (RFO)\n";
  cout << "(4) Schlegel's quasi-NR optimizer\n";
  cout << "(5) Constrained Optimization (for saddle points, initial geometry must be a minimum!)\n";
  cin >> optmet;
 }
 switch(optmet){
  case 1: break;
  case 2: **outf << "METHOD=NR ";
          break;
  case 3: **outf << "METHOD=RFO ";
          break;
  case 4: **outf << "METHOD=SCHLEGEL ";
          break;
  case 5: **outf << "METHOD=CONOPT ";
          break;
 }
 double opttol = 0.0001;
 do{
  cout << "Insert Value for Gradient Convergence Tolerance (default: 0.0001):";
  cin >> opttol;
 }while(opttol <= 0);
 if(opttol != 0.0001) **outf << "OPTTOL=" << opttol << " ";
 int nstep = 20;
 do{
  cout << "Insert maximal Number of optimization steps (default:20):";
  cin >> nstep;
 }while(nstep < 1);
 if(nstep != 20) **outf << "NSTEP=" << nstep << " ";
 **outf << "$END\n";
}

void trudge(ofstream **outf){
 **outf << " $TRUDGE ";
 int optmiz = 0;
 while(optmiz < 1 || optmiz > 2){
  cout << "What do you want to be optimized:\n";
  cout << "(1) Geometry (default)\n";
  cout << "(2) Basis Set\n";
  cin >> optmiz;
 }
 switch(optmiz){
  case 1: **outf << "OPTMIZ=GEOMETRY ";
          break;
  case 2: **outf << "OPTMIZ=BASIS ";
          break;
 }
 int npar = 0;
 while(npar < 1){
  cout << "Number of parameters to be optimized:";
  cin >> npar;
 }
 **outf << "NPAR=" << npar << " ";
 for(int x = 0; x < npar; x++){
  int index_atom;
  int index_par;
  cout << "Number of atom:";
  cin >> index_atom;
  while(index_par < 1 || index_par > 3){
   cout << "What parameter:\n";
   cout << "(1) Bond Length\n";
   cout << "(2) Bond Angle\n";
   cout << "(3) Dihedral Angle\n";
   cin >> index_par;
  }
  double param_value;
  cout << "Insert value for parameter " << 10*index_atom+index_par << ":";
  cin >> param_value;
  **outf << "IEX(" << x+1 << ")=" << 10*index_atom+index_par << " ";
  **outf << "P(" << x+1 << ")=" << param_value << " ";
 }
 **outf << "$END\n";
}


 
