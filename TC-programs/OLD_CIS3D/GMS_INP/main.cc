#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

/***********************************************************************************************************************************************************************
* GMS_INP                                                                                                                                                              *
* Program to generate GAMESS(US) input files from xyz files for usual quantum chemical calculations                                                                    *
*                                                                                                                                          Stefan Klinkusch, 2008-2009 *
***********************************************************************************************************************************************************************/

using namespace std;

//Extern Functions

extern double calc_oz(char ata, char atb);                        //Bestimmung der Ordnungszahl bei gegebenem Symbol, bestehend aus den Buchstaben ata und atb
extern void contrl_group(int nroa, int nroe, int mult_gs, int charge, ofstream *outf);  //Erstellung der $CONTRL- und weiterer Gruppen ($CIS, $DFT, $TDDFT, $MASS, $FORCE)
extern void basis_group(ofstream *outf);     //Erstellung der $BASIS-Gruppe

int main(int argc, char* argv[]){
 if(argc != 3){
  cerr << "Need xyz-file, output-prefix!\n";
  exit(1);
 }

 int nroa;          //nr of atoms
 ifstream inf;      //filestream to read in nroa und coordinates
 ofstream outf;     //filestream to write GAMESS input file
 inf.open(argv[1]); //open xyz-file
 char inpfile[256]; //name of GAMESS input file
 char titel[256];   //title of GAMESS calculation
 inf >> nroa;       //read in nroa
 double ordnr;      //charge number of nucleus
 int nroe = 0;      //nr of electrons
 int charge = 0;    //charge of molecule
 double xcoord, ycoord, zcoord; //x-,y- and z-coordinate for an atom

 sprintf(inpfile, "%s.inp", argv[2]); //generation of file name
 outf.open(inpfile);                  //open GAMESS input file
 cout << "Insert title:";   
 cin >> titel;                        //read in title from terminal
 char* atomtyp = new char[2];         //symbol for an atom
 for(int i = 0; i < nroa; i++){       //routine to compute the nr of electrons for a molecule (chargeless)
  inf >> atomtyp >> xcoord >> ycoord >> zcoord; 
  ordnr = calc_oz(atomtyp[0], atomtyp[1]);
  nroe += (int) ordnr;
}
 cout << "Insert charge of system: ";
 cin >> charge;                      //read in total charge of system (cations positive, anions negative)
 nroe = nroe - charge;               //calculation of final nr of electrons
 
 int mult_gs = 0;           //multiplicity of ground state
 
 while(mult_gs == 0){
 cout << "Insert multiplicity:";
 cin >> mult_gs;           //read in multiplicity
 if((nroe%2 == 1 && mult_gs%2 == 1) || (nroe%2 == 0 && mult_gs%2 == 0)){    //check multiplicity
  cerr << "Wrong multiplicity!\n";
  mult_gs = 0;
  }
 }

//$CONTRL group, ($CIS, $DFT, $TDDFT, $FORCE, $MASS)
 
 contrl_group(nroa, nroe, mult_gs, charge, &outf);        //construct $CONTRL group
 
//$SYSTEM group

 int mbytes = -1;
 while(mbytes < 1){
 cout << "Insert memory (in MB) to be used:";
 cin >> mbytes;
 }
 double mwordsd = 1024.*1024.* (double) mbytes/(64000000.);  //routine to calculate MWORDS from MB
 int mwords = (int) mwordsd;
 outf << " $SYSTEM MWORDS=" << mwords << " $END\n";
 
//$BASIS group

 basis_group(&outf);

//$DATA group

 outf << " $DATA\n";
 outf << titel << "\n";
 outf << "C1\n";
 inf.close();
 inf.open(argv[1]);
 int testx;
 inf >> testx;
for(int i = 0; i < nroa; i++){
 inf >> atomtyp >> xcoord >> ycoord >> zcoord;
 outf << atomtyp << ", \t" << calc_oz(atomtyp[0],atomtyp[1]) << ", \t" << xcoord << ", \t" << ycoord << ", \t" << zcoord << "\n"; //write coordinates
}
 outf << "$END\n";
 inf.close();
 outf.close();
}

