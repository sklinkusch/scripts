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

void read_integrals_12(char* f08file, int nrofint, int nrorec, int ipr, unsigned short* intnums, double* intval);

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Usage: ./readf08 f08file nrofint nrorec ipr outfile\n";
	exit(0);
    }
    char f08file[512];
    sprintf(f08file, "%s", argv[1]);
    int nrofint = atoi(argv[2]);
    int nrorec = atoi(argv[3]);
    int ipr = atoi(argv[4]);
    char outfile[512];
    sprintf(outfile, "%s", argv[5]);
    unsigned short* intnums = new unsigned short[nrofint*4];
    double* intval = new double[nrofint];

    read_integrals_12(f08file, nrofint, nrorec, ipr, intnums, intval);

    ofstream outf;
    outf.open(outfile);
    char dumv[2048];
    for(int x = 0; x < nrofint; x++){
	sprintf(dumv, "%6d %2d %2d %2d %2d % 10.6e", x, intnums[x*4+0], intnums[x*4+1], intnums[x*4+2], intnums[x*4+3], intval[x]);
	outf << dumv << "\n";
	outf.flush();
    }
    outf.close();
}


void read_integrals_12(char* f08file, int nrofint, int nrorec, int ipr, unsigned short* intnums, double* intval){
 char     dumchar[256];
 unsigned char     num_buf[15000*4];   //(!!!)
 double  val_buf[15000];

 long long int count = 0;

 ifstream inf(f08file);

 for(int x = 0; x < nrorec; x++){
    inf.read(dumchar, 16);   //16  bytes junk at start ??
    inf.read((char *) num_buf,ipr*4);               //read charakters
    inf.read((char *) val_buf,ipr*sizeof(double));  //read doubles

    int buffcount = 0;

    int offset = ipr;

    if(count + ipr >  nrofint) offset = nrofint - count;

    for(long long int y = count; y < count+offset; y++)
     intval[y] = val_buf[buffcount++];

     buffcount = 0;
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
      inf.read( dumchar , 8); //4 bytes junk at end ??
     }
}

