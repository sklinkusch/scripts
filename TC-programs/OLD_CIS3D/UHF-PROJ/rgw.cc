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

void read_line(ifstream* inf, char* line_buf, int max_count){
  int count = 0;
  char c;
  int end = 0;
  while(end == 0){
    inf->get(c);
    line_buf[count++] = c;
    if(c == '\n') end = 1;
    if(count == max_count){
      cerr << "Buffer overflow while parsing\n"; exit(1);
    }
  }
  line_buf[count-1] = 0;
}

int main(int argc, char* argv[]){
  if(argc != 3){
    cerr << "Need nroao dat-file\n";
  }
  
  int nroao = atoi(argv[1]);
  double* MOens = new double[nroao];
  double* MOs   = new double[nroao*nroao];
  for(int x = 0; x < nroao; x++) MOens[x] = 0.;
  
  ifstream inf(argv[2]);
  char numchar[16];
  numchar[15] = 0;
  
  clog << "Reading from " << argv[2] << " Nroao is " << nroao << "\n";
  char line_buf[64000];
  line_buf[1] = 0;
  
  while(strcmp(" $VEC",line_buf)!=0)
    read_line(&inf, line_buf,64000);
  
  int MOcount = 0;
  int AOcount = 0;
  while(strcmp(" $END",line_buf)!=0){
    read_line(&inf, line_buf,64000);
    if(strcmp(" $END",line_buf)!=0){
    int line_count = 0;
    while(line_count < 5){
      for(int x = 0; x < 15; x++)
	numchar[x] = line_buf[5+x+15*line_count];
      line_count++;
      MOs[MOcount*nroao+AOcount] =  strtod(numchar, NULL);
      AOcount++;
      if(AOcount == nroao){
	line_count = 5; 
	AOcount = 0;
	MOcount++;
      }
    }
    }
  }
  clog << "MOcount is " <<  MOcount << "\n";
  sprintf(line_buf,"%s.hwf",argv[2]);
  ofstream outf(line_buf);
  outf.write((char *) &nroao, sizeof(int));
  outf.write((char *) MOens,  sizeof(double)*nroao);
  outf.write((char *) MOs,    sizeof(double)*nroao*nroao);

  outf.close();

}
