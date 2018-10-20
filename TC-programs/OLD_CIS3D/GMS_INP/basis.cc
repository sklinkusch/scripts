#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*******************************************************************************************************************************************
* GMS_INP                                                                                                                                  *
* Routines to construct the $BASIS group                                                                                                   *
*                                                                                                              Stefan Klinkusch, 2008-2009 *
*******************************************************************************************************************************************/


using namespace std;

//Functions
void basis_group(ofstream *outf);
void pople_base(ofstream **outf);
void dunning_base(ofstream **outf);
void huzinaga_base(ofstream **outf);
void semiempirical_base(ofstream **outf);
void mzv_base(ofstream **outf);
void jensen_base(ofstream **outf);

//routine to choose one group of basis sets

void basis_group(ofstream *outf){
 int nrobs = 0;
 while(nrobs < 1 || nrobs > 8){
  cout << "Insert Nr of Basis-Set Group:\n";
  cout << "(1) Pople's Basis sets (e.g. STO, 3-21G, 6-31G*)\n";
  cout << "(2) Dunning's Basis sets (e.g. cc-pVDZ)\n";
  cout << "(3) Huzinaga's Basis sets (MINI, MIDI)\n";
  cout << "(4) semiempirical Basis sets (PM3, AM1, MNDO)\n";
  cout << "(5) Multiple Zeta Valence Basis sets (DZV, TZV)\n";
  cout << "(6) Dunning/Hay Basis set (DH)\n";
  cout << "(7) McLean/Chandler Basis set (MC)\n";
  cout << "(8) Jensen Polarization Consistent Basis sets (pc-VnZ, aug-pc-VnZ)\n";
  cin >> nrobs;
 }
 *outf << " $BASIS ";
 switch(nrobs)
 {
  case 1: pople_base(&outf);
          break;
  case 2: dunning_base(&outf);
          break;
  case 3: huzinaga_base(&outf);
          break;
  case 4: semiempirical_base(&outf);
          break;
  case 5: mzv_base(&outf);
          break;
  case 6: *outf << "GBASIS=DH ";
          break;
  case 7: *outf << "GBASIS=MC ";
          break;
  case 8: jensen_base(&outf);
          break;
 }

 *outf << "$END\n";
}

//routine to choose a Pople basis set

void pople_base(ofstream **outf){
 int nropb = 0;
 while(nropb < 1 || nropb > 4){
  cout << "Insert Nr of Pople-Basis-Set:\n";
  cout << "(1) STO-NG\n";                     //STO-2G, STO-3G, also STO-3G* etc.
  cout << "(2) N-21G\n";                      //e.g. 3-21G
  cout << "(3) N-31G\n";                      //e.g. 6-31G
  cout << "(4) N-311G\n";                     //e.g. 6-311G
  cin >> nropb;
 }
 int ngauss = 0;
 while(ngauss <= 0){
  cout << "No of Gaussian functions:";        //decision, whether 4-31G, 5-31G, 6-31G, ...
  cin >> ngauss;
 }
 switch(nropb)
 {
  case 1: **outf << "GBASIS=STO NGAUSS=" << ngauss << " ";
          break;
  case 2: **outf << "GBASIS=N21 NGAUSS=" << ngauss << " ";
          break;
  case 3: **outf << "GBASIS=N31 NGAUSS=" << ngauss << " ";
          break;
  case 4: **outf << "GBASIS=N311 NGAUSS=" << ngauss << " ";
          break;
 }

 if(nropb > 1){
  int pol_ha = -1;
  while(pol_ha !=0 && pol_ha !=1){
   cout << "Add polarization functions to heavy atoms (0 no, 1 yes):\n";       //makes 6-31G* from 6-31G
   cin >> pol_ha;
  }
  if(pol_ha == 1) **outf << "NDFUNC=1 ";

  if(pol_ha == 1){
   int pol_la = -1;
   while(pol_la != 0 && pol_la !=1){
    cout << "Add polarization functions to H/He (0 no, 1 yes):\n";              //makes 6-31G** from 6-31G*
    cin >> pol_la;
   }
   if(pol_la == 1) **outf << "NPFUNC=1 ";
  }

  int dif_ha = -1;
  while(dif_ha != 0 && dif_ha != 1){
   cout << "Add diffuse functions to heavy atoms (0 no, 1 yes):\n";           //makes 6-31+G from 6-31G
   cin >> dif_ha;
  }
  if(dif_ha == 1) **outf << "DIFFSP=.T. ";

  if(dif_ha == 1){
   int dif_la = -1;
   while(dif_la != 0 && dif_la != 1){
    cout << "Add diffuse functions to H/He (0 no, 1 yes):\n";                  //makes 6-31++G from 6-31+G
    cin >> dif_la;
   }
   if(dif_la == 1) **outf << "DIFFS=.T. ";
  }
 }
}

//routine to choose a Dunning basis set

void dunning_base(ofstream **outf){
 int nrodb = 0;
 while(nrodb < 1 || nrodb > 6){
  cout << "Insert No of Dunning Basis Set:\n";
  cout << "(1) cc-pVDZ\n";
  cout << "(2) cc-pVTZ\n";
  cout << "(3) cc-pVQZ\n";
  cout << "(4) aug-cc-pVDZ\n";
  cout << "(5) aug-cc-pVTZ\n";
  cout << "(6) aug-cc-pVQZ\n";
  cin >> nrodb;
 }

 switch(nrodb)
 {
  case 1: **outf << "GBASIS=CCD ";
          break;
  case 2: **outf << "GBASIS=CCT ";
          break;
  case 3: **outf << "GBASIS=CCQ ";
          break;
  case 4: **outf << "GBASIS=ACCD ";
          break;
  case 5: **outf << "GBASIS=ACCT ";
          break;
  case 6: **outf << "GBASIS=ACCQ ";
          break;
 }
}

//routine to choose a Huzinaga basis set

void huzinaga_base(ofstream **outf){
 int nrohb = 0;
 while(nrohb < 1 || nrohb > 2){
  cout << "Insert No of Huzinaga Basis Set:\n";
  cout << "(1) MINI\n";
  cout << "(2) MIDI\n";
  cin >> nrohb;
 }

 switch(nrohb)
 {
  case 1: **outf << "GBASIS=MINI ";
          break;
  case 2: **outf << "GBASIS=MIDI ";
          break;
 }
}

//routine to choose a semiempirical model hamiltonian

void semiempirical_base(ofstream **outf){
 int nrosb = 0;
 while(nrosb < 1 || nrosb > 3){
  cout << "Insert No of semiempirical Basis set:\n";
  cout << "(1) MNDO\n";
  cout << "(2) AM1\n";
  cout << "(3) PM3\n";
  cin >> nrosb;
 }

 switch(nrosb)
 {
  case 1: **outf << "GBASIS=MNDO ";
          break;
  case 2: **outf << "GBASIS=AM1 ";
          break;
  case 3: **outf << "GBASIS=PM3 ";
          break;
 }
}

//routine to choose a multiple-zeta valence basis set

void mzv_base(ofstream **outf){
 int nromzv = 0;
 while(nromzv < 1 || nromzv > 2){
  cout << "Insert No of Multiple Zeta Valence Basis set:\n";
  cout << "(1) DZV\n";
  cout << "(2) TZV\n";
  cin >> nromzv;
 }
 switch(nromzv)
 {
  case 1: **outf << "GBASIS=DZV ";
          break;
  case 2: **outf << "GBASIS=TZV ";
          break;
 }
}

//routine to choose a Jensen basis set

void jensen_base(ofstream **outf){
 int nrojb = 0;
 while(nrojb < 1 || nrojb > 10){
  cout << "Insert No of Jensen Basis set:\n";
  cout << "(1) PC-0 unpolarized\n";
  cout << "(2) PC-1\n";
  cout << "(3) PC-2\n";
  cout << "(4) PC-3\n";
  cout << "(5) PC-4\n";
  cout << "(6) APC-0\n";
  cout << "(7) APC-1\n";
  cout << "(8) APC-2\n";
  cout << "(9) APC-3\n";
  cout << "(10) APC-4\n";
  cin >> nrojb;
 }
 switch(nrojb)
 {
  case 1: **outf << "GBASIS=PC0 ";
          break;
  case 2: **outf << "GBASIS=PC1 ";
          break;
  case 3: **outf << "GBASIS=PC2 ";
          break;
  case 4: **outf << "GBASIS=PC3 ";
          break;
  case 5: **outf << "GBASIS=PC4 ";
          break;
  case 6: **outf << "GBASIS=APC0 ";
          break;
  case 7: **outf << "GBASIS=APC1 ";
          break;
  case 8: **outf << "GBASIS=APC2 ";
          break;
  case 9: **outf << "GBASIS=APC3 ";
          break;
  case 10: **outf << "GBASIS=APC4 ";
           break;
 }
}

