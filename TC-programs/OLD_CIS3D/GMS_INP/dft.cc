#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********************************************************************************************************************************
* GMS_INP                                                                                                                       *
* routines necessary for density functional calculations                                                                        *
*                                                                                                   Stefan Klinkusch, 2008-2009 *
********************************************************************************************************************************/

using namespace std;

//Functions
void dft_insert(ofstream **outf);
void dft_props(ofstream **outf);
void tddft_props(ofstream **outf);
void dft_exch(ofstream ***outf);
void dft_corr(ofstream ***outf);
void dft_comb(ofstream ***outf);
void dft_hybr(ofstream ***outf);
void dft_meta(ofstream ***outf);

// routine to insert method in $CONTRL group

void dft_insert(ofstream **outf){
 int dft_group = 0;
 while(dft_group < 1 || dft_group > 5){
  cout << "Insert Number of DFT Method Group to be used:\n";
  cout << "(1) Pure exchange (Slater, Becke, Gill, OPTX, PW91X, PBEX)\n";
  cout << "(2) Pure correlation (VWN(5), VWN1, PZ81, P86, LYP, PW91C, PBEC, OP)\n";
  cout << "(3) Combinations (SVWN, BLYP, BOP, BP86, GVWN, GPW91, PBEVWN, PBE0P, OLYP, PW91, PBE)\n";
  cout << "(4) Hybrid GGA (BHHLYP, B3LYP, B3LYP1, PBE0, X3LYP)\n";
  cout << "(5) Meta GGA (VS98, M05, M05-2X, M06, M06-L, M06-2X, M06-HF)\n";
  cin >> dft_group;
 }
 switch(dft_group)
 {
  case 1: dft_exch(&outf);
          break;
  case 2: dft_corr(&outf);
          break;
  case 3: dft_comb(&outf);
          break;
  case 4: dft_hybr(&outf);
          break;
  case 5: dft_meta(&outf);
          break;
 }
}

//routine to choose one of the pure exchange functionals

void dft_exch(ofstream ***outf){
 int nroexm = 0;
 while(nroexm < 1 || nroexm > 6){
  cout << "Insert No of functional to be used:\n";
  cout << "(1) Slater\n";
  cout << "(2) Becke\n";
  cout << "(3) Gill\n";
  cout << "(4) OPTX (Handy-Cohen)\n";
  cout << "(5) PW91 (Perdew-Wang 1991)\n";
  cout << "(6) PBEX (Perdew-Burke-Ernzerhof)\n";
  cin >> nroexm;
 }
 switch(nroexm)
 {
  case 1: ***outf << "SLATER ";
          break;
  case 2: ***outf << "BECKE ";
          break;
  case 3: ***outf << "GILL ";
          break;
  case 4: ***outf << "OPTX ";
          break;
  case 5: ***outf << "PW91 ";
          break;
  case 6: ***outf << "PBEX ";
          break;
 }
}

//routine to choose one of the pure correlation functionals

void dft_corr(ofstream ***outf){
 int nrocrm = 0;
 while(nrocrm < 1 || nrocrm > 8){
  cout << "Insert No of functional to be used:\n";
  cout << "(1) VWN (Vosko-Wilke-Nusair formula 5)\n";
  cout << "(2) VWN1 (Vosko-Wilke-Nusair formula 1)\n";
  cout << "(3) PZ81 (Perdew-Zener 1981)\n";
  cout << "(4) P86 (Perdew 1986)\n";
  cout << "(5) LYP (Lee-Yang-Parr)\n";
  cout << "(6) PW91C (Perdew-Wang 1991)\n";
  cout << "(7) PBEC (Perdew-Burke-Ernzerhof)\n";
  cout << "(8) OP (One-parameter progressive)\n";
  cin >> nrocrm;
 }
 switch(nrocrm)
 {
  case 1: ***outf << "VWN ";
          break;
  case 2: ***outf << "VWN1 ";
          break;
  case 3: ***outf << "PZ81 ";
          break;
  case 4: ***outf << "P86 ";
          break;
  case 5: ***outf << "LYP ";
          break;
  case 6: ***outf << "PW91C ";
          break;
  case 7: ***outf << "PBEC ";
          break;
  case 8: ***outf << "OP ";
          break;
 }
}

//routine to choose one of the combined exchange-correlation functionals

void dft_comb(ofstream ***outf){
 int nrocm = 0;
 while(nrocm < 1 || nrocm > 11){
  cout << "Insert No of functional to be used:\n";
  cout << "(1) SVWN (Slater exchange + Vosko-Wilke-Nusair correlation)\n";
  cout << "(2) BLYP (Becke exchange + Lee-Yang-Parr correlation)\n";
  cout << "(3) BOP (Becke exchange + One-parameter correlation)\n";
  cout << "(4) BP86 (Becke exchange + Perdew 1986 correlation)\n";
  cout << "(5) GVWN (Gill exchange + Vosko-Wilke-Nusair correlation)\n";
  cout << "(6) GPW91 (Gill exchange + Perdew-Wang 1991 correlation)\n";
  cout << "(7) PBEVWN (Perdew-Burke-Ernzerhof exchange + Vosko-Wilke-Nusair correlation)\n";
  cout << "(8) PBEOP (Perdew-Burke-Ernzerhof exchange + One-point correlation)\n";
  cout << "(9) OLYP (Handy-Cohen OPTX exchange + Lee-Yang-Parr correlation)\n";
  cout << "(10) PW91 (Perdew-Wang 1991 exchange + Perdew-Wang 1991 correlation)\n";
  cout << "(11) PBE (Perdew-Burke-Ernzerhof exchange + Perdew-Burke-Ernzerhof correlation)\n";
  cin >> nrocm;
 }
 switch(nrocm)
 {
  case 1: ***outf << "SVWN ";
          break;
  case 2: ***outf << "BLYP ";
          break;
  case 3: ***outf << "BOP ";
          break;
  case 4: ***outf << "BP86 ";
          break;
  case 5: ***outf << "GVWN ";
          break;
  case 6: ***outf << "GPW91 ";
          break;
  case 7: ***outf << "PBEVWN ";
          break;
  case 8: ***outf << "PBEOP ";
          break;
  case 9: ***outf << "OLYP ";
          break;
  case 10: ***outf << "PW91 ";
           break;
  case 11: ***outf << "PBE ";
           break;
 }
}

//routine to choose one of the hybrid functionals (including HF-like exchange)

void dft_hybr(ofstream ***outf){
 int nrohyb = 0;
 while(nrohyb < 1 || nrohyb > 5){
  cout << "Insert No of functional to be used:\n";
  cout << "(1) BHHLYP (Becke/HF exchange + Lee-Yang-Parr correlation)\n";
  cout << "(2) B3LYP (Becke/Slater/HF exchange + Lee-Yang-Parr/Vosko-Wilke-Nusair correlation)\n";
  cout << "(3) B3LYP1 (Becke/Slater/HF exchange + Lee-Yang-Parr/Vosko-Wilke-Nusair (form.1) correlation)\n";
  cout << "(4) PBE0 (Perdew-Burke-Ernzerhof/HF exchange + Perdew-Burke-Ernzerhof correlation)\n";
  cout << "(5) X3LYP (HF/Slater/Becke 1988/Perdew-Wang 1991 exchange + Lee-Yang-Parr/Vosko-Wilke-Nusair (form. 1) correlation)\n";
  cin >> nrohyb;
 }
 switch(nrohyb)
 {
  case 1: ***outf << "BHHLYP ";
          break;
  case 2: ***outf << "B3LYP ";
          break;
  case 3: ***outf << "B3LYP1 ";
          break;
  case 4: ***outf << "PBE0 ";
          break;
  case 5: ***outf << "X3LYP ";
          break;
 }
}

//routine to choose one of the meta functionals

void dft_meta(ofstream ***outf){
 int nromet = 0;
 while(nromet < 1 || nromet > 7){
  cout << "Insert No of functional to be used:\n";
  cout << "(1) VS98 (Voorhis-Scuseria 1998 exchange-correlation)\n";
  cout << "(2) M05 (Minnesota 2005 exchange-correlation, 28% HF exchange)\n";
  cout << "(3) M05-2X (Minnesota 2005 exchange-correlation, 56% HF exchange)\n";
  cout << "(4) M06 (Minnesota 2006 exchange-correlation, 27% HF exchange)\n";
  cout << "(5) M06-L (Minnesota 2006 exchange-correlation, local)\n";
  cout << "(6) M06-2X (Minnesota 2006 exchange-correlation, 54% HF exchange)\n";
  cout << "(7) M06-HF (HF exchange + Minnesota correlation)\n";
  cin >> nromet;
 }
 switch(nromet)
 {
  case 1: ***outf << "VS98 ";
          break;
  case 2: ***outf << "M05 ";
          break;
  case 3: ***outf << "M05-2X ";
          break;
  case 4: ***outf << "M06 ";
          break;
  case 5: ***outf << "M06-L ";
          break;
  case 6: ***outf << "M06-2X ";
          break;
  case 7: ***outf << "M06-HF ";
          break;
 }
}

//routine to construct the $DFT group
//including long-range correction with HF exchange

void dft_props(ofstream **outf){
 int lc_index = 0;
 double mu = -0.33;
 while(lc_index < 1 || lc_index > 2){
  cout << "Include long-range correction for large inter-electronic distances (works only with BLYP, BOP, and BVWN)\n";
  cout << "(1) Yes\n";
  cout << "(2) No\n";
  cin >> lc_index;
 }
 if(lc_index == 1){
  while(mu < 0 || mu > 1){
   cout << "Insert mu_parameter for HF part of exchange (default=0.33):";
   cin >> mu;
  }
  **outf << " $DFT LCFLAG=.T. EMU=" << mu << " $END\n";
 }
}

//routine to construct the $TDDFT group
//settings for a TDDFT calculation

void tddft_props(ofstream **outf){
 int nros = 0;
 while(nros < 1){
  cout << "Insert nr of excited states to be calculated:";
  cin >> nros;
 }
 int state = 0;
 while(state < 1){
  cout << "Insert excited state to bu used for geometry optimization and properties:";
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
 **outf << " $TDDFT NSTATE=" << nros << " IROOT=" << state << " MULT=" << mult_es << " $END\n";
}

