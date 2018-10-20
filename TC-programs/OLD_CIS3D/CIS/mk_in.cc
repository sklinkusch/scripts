/********************************************************************************
 * Extracts one and two electron integrals from  Gamess F08 and F10 files .     *
 * Tested only for 12 bytes per integral 15000 integrals per record.            *
 *                                                                              *
 *                                                                              *
 *                                         Tillmann Klamroth 2004               *
 ********************************************************************************/

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
#include <stdlib.h>

using namespace std;

//VARIABLES
int           nroao;         //Nr of basis functions
int           nroa;          //Nr of atoms 
int           bpi;           //bytes per integral
int           ipr;           //integrals per record
int           nrorec;        //nrofrecords
long long int nrofint;       //nr of nonzero to electron integrals
long long int sortcount[4];  //boundaries for different permutation patterns;


double*              intval;        //two electron integrals
unsigned short int*  intnums;       //two electron indecies
double*              coord;         //atomic ccordinates
double*              charges;       //atomic charges  (all set to the same value so far)
double*              mass;          //atomic masses
double*              hmat;          //one electron hamiltonian
double*              kmat;          //kinetic energy
double*              smat;          //overlapp
double*              Dx;            //dipole x
double*              Dy;            //dipole y
double*              Dz;            //dipole z 
char                 f10file[256];  //Gamess f10 file
char                 f08file[256];  //Gamess f08 file


//EXTERN FUNCTIONS
extern "C" void dafrd_(double* V, long long int* LEN, long long int* RECN, char* FNAME);

//FUNCTIONS
int rem_com(char* filename, char* streamstring, int string_length);
void read_gammes(double* data, int recnr, int nromo, int type, char* filename);
void read_integrals_12(void);
void read_integrals_16(void);
void resort_integrals(ofstream *outf);


int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need input-file output-prefix\n";
    exit(1);
  }
  
  char dumc[2048];
  sprintf(dumc,"%s.log",argv[2]);
  ofstream outf(dumc);
  
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Reading input form " << argv[1] << "\n";

  int buff_length  =  65536;
  char*  file_buff =  new char[buff_length];
  rem_com(argv[1], file_buff,  buff_length);
  istringstream ist(file_buff);
  
  double allcharge;   //value of all atom core charges !!!!
  double allmass;     // value of all atom core masses !!!!
  
  int extra_charge = 0; 

  ist >> f08file >> f10file >> nroao >>  nroa >> allcharge >> allmass >> bpi >> ipr >> nrorec >> nrofint;

  outf << "-------------------------------------------------------------------------------\n";
  outf << "F08-file : " << f08file << "\n";
  outf << "F10-file : " << f10file << "\n";
  outf << "Nr of basis functions: " << nroao << "\n";
  outf << "Nr of atoms: " << nroa << "\n";
  outf << "Charge of atom cores: " << allcharge << "\n";
  outf << "Mass of atom cores: " << allmass << "\n";
  if(allcharge < 0.){
    allcharge = -allcharge;
    outf << "Default charge and mass set to " << allcharge << " " <<  allmass << "\n";
    outf << "Extra charges and masses read in at end of input\n";
    extra_charge = 1;
  }
  outf << "...............................................................................\n";
  outf << "2el Data:\n";
  outf << "Bytes per integral: " << bpi << "\n";
  outf << "Integrals per record: " << ipr << "\n";
  outf << "Nr of records: " << nrorec << "\n";
  outf << "Nr of 2el-Integrals: " << nrofint << "\n";
  outf << "-------------------------------------------------------------------------------\n";

  if(bpi != 12 && bpi != 16){
    outf << bpi << " Bytes per Integral not implemented yet\n";
    outf.flush();
    exit(2);
  }

  outf << "Allocating Memory\n";
  
  int aointl = nroao*nroao*6+nroa*5;
  outf << "Need " <<  aointl*8 << " Bytes for one el ints, coords, and charges.\n";
  outf << "Need " <<  nrofint*8 << " bytes for two el indices and the same for two el values.\n";
  
  double* dumd = new double[aointl+1024];
  
  int incre = 0;
  coord   = &(dumd[incre]); incre += 3*nroa;
  charges = &(dumd[incre]); incre +=   nroa;    
  mass    = &(dumd[incre]); incre +=   nroa;    
  hmat    = &(dumd[incre]); incre += nroao*nroao;       
  kmat    = &(dumd[incre]); incre += nroao*nroao;       
  smat    = &(dumd[incre]); incre += nroao*nroao;       
  Dx      = &(dumd[incre]); incre += nroao*nroao;       
  Dy      = &(dumd[incre]); incre += nroao*nroao;       
  Dz      = &(dumd[incre]); incre += nroao*nroao;       
  
  intval = new double[nrofint];
  intnums = new unsigned short[nrofint*4+4];  //4 numbers extra for big litlle endian quatsch

  outf << "Done\n";
  outf << "-------------------------------------------------------------------------------\n";
  
  outf.flush();

  outf << "Adding zeros to F10 file\n";
  ofstream outb(f10file,ios::app);
  char dc = 0;
  for(int x = 0; x < 4090*8; x++)
    outb << dc; 
  outb.close();
  
  outf << "Reading one electron matrices\n";

  read_gammes(coord , 1, 3*nroa, 0,f10file); 
  
  read_gammes(hmat , 11, nroao, 1,f10file);
  read_gammes(smat , 12, nroao, 1,f10file);
  read_gammes(kmat , 13, nroao, 1,f10file);
  read_gammes(Dx   , 95, nroao, 1,f10file);
  read_gammes(Dy   , 96, nroao, 1,f10file);
  read_gammes(Dz   , 97, nroao, 1,f10file);

  outf << "Done\n";
  outf << "-------------------------------------------------------------------------------\n";
  outf << "Atomic coordinates:\n";
  
  int count = 0;
  char  dumchar[256];
  for(int x = 0; x < nroa; x++){
    outf << x << "\t";
    sprintf(dumchar,"%+.5f", coord[count++]);    outf  << dumchar << "\t";
    sprintf(dumchar,"%+.5f", coord[count++]);    outf  << dumchar << "\t";
    sprintf(dumchar,"%+.5f", coord[count++]);    outf  << dumchar << "\n";
    charges[x] = allcharge;
    mass[x] = allmass;
  }

  outf << "All charges set to " << allcharge << "\n";
  outf << "All masses set  to " << allmass   << "\n";
  outf << "-------------------------------------------------------------------------------\n";
  if(extra_charge == 1){
    int nroex, curr_ex;
    ist >> nroex;
    outf << "Reading " << nroex << " non default charges and masses from input\n";
    for(int x = 0; x < nroex; x++){
      ist >> curr_ex;
      if(curr_ex >= nroa){
	outf << "Center nr out of range!!\n";
	exit(5);
      }
      ist >> charges[curr_ex] >> mass[curr_ex];
      outf << "Charge/Mass for center nr: " << curr_ex  << " " <<  charges[curr_ex] << " " << mass[curr_ex] << "\n";
    }
    outf << "-------------------------------------------------------------------------------\n";
  }
  outf.flush();
    
  
  outf << "Reading 2el integrals\n";
  
  if(bpi == 12) 
    read_integrals_12();

  if(bpi == 16) 
    read_integrals_16();

  
  outf << "Done \n";
  outf.flush();
  outf << "-------------------------------------------------------------------------------\n";

  
  resort_integrals(&outf);  

  sprintf(dumc,"%s.sys",argv[2]);
  ofstream datf(dumc);

  outf << "Writing binary data to " <<  dumc << "\n";
  
  //SYSTEM DATA
  datf.write((char *) &nroao , sizeof(int));
  datf.write((char *) &nroa  , sizeof(int));
  datf.write((char *) &nrofint,  sizeof(long long int));
  datf.write((char *) coord  , sizeof(double)*3*nroa);
  datf.write((char *) charges, sizeof(double)*nroa);
  datf.write((char *) mass,    sizeof(double)*nroa);
  
  //ONEL EL INTEGRAL DATA
  datf.write((char * ) hmat  , sizeof(double)*nroao*nroao);
  datf.write((char * ) kmat  , sizeof(double)*nroao*nroao);
  datf.write((char * ) smat  , sizeof(double)*nroao*nroao);
  datf.write((char * ) Dx    , sizeof(double)*nroao*nroao);
  datf.write((char * ) Dy    , sizeof(double)*nroao*nroao);
  datf.write((char * ) Dz    , sizeof(double)*nroao*nroao);
  
  //TWO EL INTEGRAL DATA
  datf.write((char *) sortcount, sizeof(long long int)*4);
  datf.write((char *) intval,    sizeof(double)*nrofint);
  datf.write((char *) intnums,   sizeof(unsigned short)*nrofint*4);
  
  datf.close();
  outf << "...............................................................................\n";
}


void read_integrals_12(void){
           char     dumchar[256];
  unsigned char     num_buf[15000*4];   //(!!!)
  double  val_buf[15000];

  long long int count = 0;
  
  ifstream inf(f08file);
  
  for(int x = 0; x < nrorec; x++){
    inf.read(dumchar, 12);   //12  bytes junk at start ??
    
    inf.read((char *) num_buf,ipr*4);               //read charakters
    inf.read((char *) val_buf,ipr*sizeof(double));  //read doubles

    int buffcount = 0;
    
    int offset = ipr;

    if(count + ipr >  nrofint) offset = nrofint - count;

    for(long long int y = count; y < count+offset; y++)
      intval[y] = val_buf[buffcount++];
    
    buffcount = 0;
    //Bloody high low word kacke !!!
    for(long long int y = count; y < count+offset; y+=2){
      intnums[y*4+4] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+5] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+6] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+7] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+0] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+1] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+2] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+3] = (unsigned short ) num_buf[buffcount++]-1;
    }
    
    count += offset;
    inf.read( dumchar , 4); //4 bytes junk at end ??
  }  
  
}

void read_integrals_16(void){
           char     dumchar[256];
  unsigned short     num_buf[15000*4];   //(!!!)
  double  val_buf[15000];

  long long int count = 0;
  
  ifstream inf(f08file);
  
  for(int x = 0; x < nrorec; x++){
    inf.read(dumchar, 12);   //??12 bytes junk at start 
    
    inf.read((char *) num_buf,ipr*8);               //read charakters
    inf.read((char *) val_buf,ipr*sizeof(double));  //read doubles

    int buffcount = 0;
    
    int offset = ipr;

    if(count + ipr >  nrofint) offset = nrofint - count;

    for(long long int y = count; y < count+offset; y++)
      intval[y] = val_buf[buffcount++];
    
    buffcount = 0;
    //Bloody high low word kacke !!!
    for(long long int y = count; y < count+offset; y++){
      intnums[y*4+2] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+3] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+0] = (unsigned short ) num_buf[buffcount++]-1;
      intnums[y*4+1] = (unsigned short ) num_buf[buffcount++]-1;
    }
    
    count += offset;
    inf.read( dumchar , 4); //4?? bytes junk at end 
  }  
  
}


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



void read_gammes(double* data, int recnr, int nromo, int type, char* filename){

  long long int len;
  //SET CORRECT LENGTH
  if(type == 0) len = nromo;
  if(type == 1) len = (nromo*nromo+nromo)/2;
  if(type == 2) len = nromo*nromo;

  //len = 2;

  char dumn[513];
  int fnamelen = strlen(filename);
  sprintf(dumn,"%s",filename);
  for(int x = fnamelen; x < 512; x++) dumn[x] = ' ';
  dumn[512] = 0;

  long long int Recnr = recnr;

  dafrd_(data, &len, &Recnr, dumn);


  //RESORT MATRIX
  if(type == 1){
    int dcount = len-1;
    for(int x = nromo-1; x >= 0; x--){
      for(int y = x; y >= 0; y--){
        data[x*nromo+y] = data[dcount--];
      }
    }
    for(int x = 0; x < nromo; x++){
      for(int y = x; y < nromo; y++) data[x*nromo+y] = data[y*nromo+x];
    }
  }

}

inline void swap_ints(long long int a, long long int b, double* intvals, unsigned short* intnums){
  //temp = a
  unsigned short ta = intnums[a*4+0];
  unsigned short tb = intnums[a*4+1];
  unsigned short tc = intnums[a*4+2];
  unsigned short td = intnums[a*4+3];
  double tval = intvals[a];
  
  //a = b
  intnums[a*4+0] =   intnums[b*4+0];     
  intnums[a*4+1] =   intnums[b*4+1];  
  intnums[a*4+2] =   intnums[b*4+2];  
  intnums[a*4+3] =   intnums[b*4+3];  
  intvals[a]     =   intvals[b]    ;    
  
  //b = temp
  intnums[b*4+0] =    ta; 
  intnums[b*4+1] =    tb; 
  intnums[b*4+2] =    tc; 
  intnums[b*4+3] =    td; 
  intvals[b]     =    tval;
  
}


void resort_integrals(ofstream *outf){
    //RESORT INTEGRALS

  //STEP 1: BRING to basis types
  for(long long int x = 0; x < nrofint; x++){
    unsigned short a,b,c,d;
    
    //INPUT AUCH  CHEMIKER ??
    a      = intnums[x*4+0];  //a
    b      = intnums[x*4+1];  //b
    c      = intnums[x*4+2];  //c
    d      = intnums[x*4+3];  //d
    
    
    //TYP IIb: reorder to IIa
    if(a==b && b==d && c!=d){
      intnums[x*4+0] = a;  //a
      intnums[x*4+1] = a;  //b
      intnums[x*4+2] = a;  //c
      intnums[x*4+3] = c;  //d
    }

    //TYP IIc: reorder to IIa
    if(a!=b && a==c && c==d){
      intnums[x*4+0] = a;  //a
      intnums[x*4+1] = a;  //b
      intnums[x*4+2] = a;  //c
      intnums[x*4+3] = b;  //d
    }

    //TYP IId: reorder to IIa
    if(a!=b && b==c && c==d){
      intnums[x*4+0] = b;  //a
      intnums[x*4+1] = b;  //b
      intnums[x*4+2] = b;  //c
      intnums[x*4+3] = a;  //d
    }
    
    //TYP IIIc: reorder to IIIb
    if(a==d && b==c && a!=b){
      intnums[x*4+0] = a;  //a
      intnums[x*4+1] = b;  //b
      intnums[x*4+2] = d;  //c
      intnums[x*4+3] = c;  //d      
    }
    

    //TYP IVc:  reorder IVb
    if(b==c && b!=a && b!=d && a!=d){
      intnums[x*4+0] = b;  //a
      intnums[x*4+1] = a;  //b
      intnums[x*4+2] = c;  //c
      intnums[x*4+3] = d;  //d      
    }

     //TYP IVd:  reorder to IVa
    if(c==d && c!=a && c!=b && a!=b){
      intnums[x*4+0] = c;  //a
      intnums[x*4+1] = c;  //b
      intnums[x*4+2] = a;  //c
      intnums[x*4+3] = b;  //d      
    }

    //TYP IVe:  reorder IVb
    if(a==d && a!=b && a!=c && b!=c){
      intnums[x*4+0] = a;  //a
      intnums[x*4+1] = b;  //b
      intnums[x*4+2] = d;  //c
      intnums[x*4+3] = c;  //d      
    }

    //TYP IVf: perm_all,  reorder IVb
    if(b==d && b!=a && b!=c && a!=c){
      intnums[x*4+0] = b;  //a
      intnums[x*4+1] = a;  //b
      intnums[x*4+2] = d;  //c
      intnums[x*4+3] = c;  //d      
    }    
  }

  //STEP2: BRING PERM1 to front
  long long int sorted = 0;
  for(long long int x = 0; x < nrofint; x++){
    unsigned short a,b,c,d;
    
    //INPUT AUCH  CHEMIKER ??
    a      = intnums[x*4+0];  //a
    b      = intnums[x*4+1];  //b
    c      = intnums[x*4+2];  //c
    d      = intnums[x*4+3];  //d
    if(a==b && b==c && c==d){
      swap_ints(x, sorted,  intval,  intnums);
      sorted++;
    }
  }
  *outf << sorted << " integrals after search for perm_1\n";
  sortcount[0] = sorted;
  //STEP3: BRING PERM1_5 thereafter
  
  for(long long int x = sorted; x < nrofint; x++){
    unsigned short a,b,c,d;
    
    //INPUT AUCH  CHEMIKER ??
    a      = intnums[x*4+0];  //a
    b      = intnums[x*4+1];  //b
    c      = intnums[x*4+2];  //c
    d      = intnums[x*4+3];  //d
    if(a==b && c==d){
      swap_ints(x, sorted,  intval,  intnums);
      sorted++;
    }
  }
  *outf << sorted << " integrals after search for perm_15\n";
  sortcount[1] = sorted;


  //STEP4: BRING PERM1234
  for(int x = sorted; x < nrofint; x++){
    int a,b,c,d;
    
    //INPUT AUCH  CHEMIKER ??
    a      = intnums[x*4+0];  //a
    b      = intnums[x*4+1];  //b
    c      = intnums[x*4+2];  //c
    d      = intnums[x*4+3];  //d
    if(a==c && b==d){
      swap_ints(x, sorted,  intval,  intnums);
      sorted++;
    }
  }
  *outf << sorted << " integrals after search for perm_1234\n";
  sortcount[2] = sorted;

  //STEP5: BRING PERM1256
  for(int x = sorted; x < nrofint; x++){
    int a,b,c,d;
    
    //INPUT AUCH  CHEMIKER ??
    a      = intnums[x*4+0];  //a
    b      = intnums[x*4+1];  //b
    c      = intnums[x*4+2];  //c
    d      = intnums[x*4+3];  //d
    if(a==b){
      swap_ints(x, sorted,  intval,  intnums);
      sorted++;
    }
  }
  *outf << sorted << " integrals after search for perm_1256\n";
  sortcount[3] = sorted;

  outf->flush();
  *outf << "-------------------------------------------------------------------------------\n";
  *outf << "Undoing GAMESS - scaling\n";
  
  //PERM_1
  for(long long int x = 0; x < sortcount[0]; x++)
    intval[x] *= 8.;

  //PERM_15
  for(long long int x = sortcount[0]; x < sortcount[1]; x++)
    intval[x] *= 4.;

  //PERM_1234
  for(long long int x = sortcount[1]; x < sortcount[2]; x++)
    intval[x] *= 2.;

  //PERM_1256
  for(long long int x = sortcount[2]; x < sortcount[3]; x++)
    intval[x] *= 2.;

  *outf << "-------------------------------------------------------------------------------\n";
  outf->flush();
}
