# include<iostream>
# include<fstream>
# include<math.h>
# include<stdio.h>
# include<stdlib.h>
# include<string.h>
# include<string>
# include <limits.h>
# include <time.h>
# include <unistd.h>
# include <ctype.h>
# include <sys/time.h> 
# include <sys/resource.h>
# include <sstream>

using namespace std;

//Funktionen
int  rem_com(char* filename, char* streamstring, int string_length);
void day_x(int freq, int nrod, int* day, int* month, int schaltjahr, int &count);
void day_y(int freq, int nrod, int* day, int* month, int schaltjahr, int &count, int* cday, int* cmonth, int* cfreq, int cmode);
int schaltjahr(int year);
void subdate(int nrod, int nrosd, int* day, int* month, int year, int* dayin, int* dayout, int* monthin, int* monthout, int* yearin, int* yearout, string* subtype, char typus);
bool print(string::size_type n, string s);
void cleardate(int nrod, int* day, int* month, int &count);
void zerodate(int &count);

int main(int argc, char* argv[]){
 if(argc != 2){
  cerr << "Need input-file\n";
  exit(1);
 }

int buff_length  =  65536;
char*  file_buff =  new char[buff_length];
rem_com(argv[1], file_buff,  buff_length);
istringstream ist(file_buff);


 const int nrod = 55;  // Zahl der maximalen Abholtage (55 > Anzahl der Wochen, sollte also ausreichen)
 const int nrowd = 2;  // Zahl der Weihnachtsbaumabholtage
 int alba_mode, org_mode, bsr_mode, paper_mode, glass_mode, bio_mode, tree_mode; // Modus der Abholung (0: keine Abholung, 1: Abholung, 2: eine Änderung, 3: 2 Änderungen, ...)
 int alba_freq, org_freq, bsr_freq, paper_freq, glass_freq, bio_freq; // Frequenz der Abholung (in Tagen)
 int albaccount = 0; 
 int orgccount = 0;
 int bsrccount = 0;
 int paperccount = 0;
 int glassccount = 0;
 int bioccount = 0;
 int* albaday    = new int[nrod];  int* albamonth  = new int[nrod]; // Termine zur Abholung der gelben Tonnen/gelben Säcke
 int* orgday     = new int[nrod];  int* orgmonth   = new int[nrod]; // Termine zur Abholung der orangen Wertstofftonne
 int* bsrday     = new int[nrod];  int* bsrmonth   = new int[nrod]; // Termine zur Abholung der Restmülltonne
 int* paperday   = new int[nrod];  int* papermonth = new int[nrod]; // Termine zur Abholung der Altpapiertonne
 int* glassday   = new int[nrod];  int* glassmonth = new int[nrod]; // Termine zur Abholung der Altglastonne
 int* bioday     = new int[nrod];  int* biomonth   = new int[nrod]; // Termine zur Abholung der Biomülltonne
 int* treeday    = new int[nrowd];  int* treemonth  = new int[nrowd]; // Termine zur Abholung der Weihnachtsbäume
 int  year; // Kalenderjahr
 int* albacday = new int[nrod]; int* albacmonth = new int[nrod]; int* albacfreq = new int[nrod]; // Änderungstermine
 int* orgcday = new int[nrod]; int* orgcmonth = new int[nrod]; int* orgcfreq = new int[nrod];
 int* bsrcday = new int[nrod]; int* bsrcmonth = new int[nrod]; int* bsrcfreq = new int[nrod];
 int* papercday = new int[nrod]; int* papercmonth = new int[nrod]; int* papercfreq = new int[nrod];
 int* glasscday = new int[nrod]; int* glasscmonth = new int[nrod]; int* glasscfreq = new int[nrod];
 int* biocday = new int[nrod]; int* biocmonth = new int[nrod]; int* biocfreq = new int[nrod];

 ofstream outf; // Ausgabe aller Termine nach termine.dat
 outf.open("termine.dat");

// Einlesen von Modus, erstem Abholtermin, Abholfrequenz und Änderungsterminen
ist >> alba_mode; // gelbe Säcke, gelbe Tonne
if(alba_mode == 1){
 ist >> alba_freq;
 ist >> albaday[0] >> albamonth[0];
}else if(alba_mode > 1){
  ist >> alba_freq >> albaday[0] >> albamonth[0];
 for(int i = 0; i < (alba_mode - 1); i++){
  ist >> albacday[i] >> albacmonth[i] >> albacfreq[i];
  albaccount++;
 }
}else{
 alba_freq = 0;
 albaday[0] = 0;
 albamonth[0] = 0;
}
ist >> org_mode; // orange Tonne
if(org_mode == 1){
 ist >> org_freq;
 ist >> orgday[0] >> orgmonth[0];
}else if(org_mode > 1){
  ist >> org_freq >> orgday[0] >> orgmonth[0];
 for(int i = 0; i < (org_mode - 1); i++){
  ist >> orgcday[i] >> orgcmonth[i] >> orgcfreq[i];
  orgccount++;
 }
}else{
 org_freq = 0;
 orgday[0] = 0;
 orgmonth[0] = 0;
}
ist >> bsr_mode; // Restmülltonne
if(bsr_mode == 1){
 ist >> bsr_freq;
 ist >> bsrday[0] >> bsrmonth[0];
}else if(bsr_mode > 1){
  ist >> bsr_freq >> bsrday[0] >> bsrmonth[0]; 
 for(int i = 0; i < (bsr_mode - 1); i++){
  ist >> bsrcday[i] >> bsrcmonth[i] >> bsrcfreq[i];
  bsrccount++;
 }
}else{
 bsr_freq = 0;
 bsrday[0] = 0;
 bsrmonth[0] = 0;
}
ist >> paper_mode; // Papiertonne
if(paper_mode == 1){
 ist >> paper_freq;
 ist >> paperday[0] >> papermonth[0];
}else if(paper_mode > 1){
 ist >> paper_freq >> paperday[0] >> papermonth[0];
 for(int i = 0; i < (paper_mode - 1); i++){
  ist >> papercday[i] >> papercmonth[i] >> papercfreq[i];
  paperccount++;
 }
}else{
 paper_freq = 0;
 paperday[0] = 0;
 papermonth[0] = 0;
}
ist >> glass_mode; // Glastonne
if(glass_mode == 1){
 ist >> glass_freq;
 ist >> glassday[0] >> glassmonth[0];
}else if(glass_mode > 1){
 ist >> glass_freq >> glassday[0] >> glassmonth[0];
 for(int i = 0; i < (glass_mode - 1); i++){
  ist >> glasscday[i] >> glasscmonth[i] >> glasscfreq[i];
  glassccount++;
 }
}else{
 glass_freq = 0;
 glassday[0] = 0;
 glassmonth[0] = 0;
}
ist >> bio_mode; // Biomülltonne
if(bio_mode == 1){
 ist >> bio_freq;
 ist >> bioday[0] >> biomonth[0];
}else if(bio_mode > 1){
  ist >> bio_freq >> bioday[0] >> biomonth[0];
 for(int i = 0; i < (bio_mode - 1); i++){
  ist >> biocday[i] >> biocmonth[i] >> biocfreq[i];
  bioccount++;
 }
}else{
 bio_freq = 0;
 bioday[0] = 0;
 biomonth[0] = 0;
}
ist >> tree_mode; // Abholung von Weihnachtsbäumen
if(tree_mode == 1){
 ist >> treeday[0] >> treemonth[0] >> treeday[1] >> treemonth[1];
}else{
 treeday[0] = 0;
 treemonth[0] = 0;
 treeday[1] = 0;
 treemonth[1] = 0;
}

ist >> year; // Kalenderjahr

int q = schaltjahr(year); // Feststellen, ob das Jahr ein Schaltjahr ist

int albacount; // Zahl der Abholtermine
int orgcount;
int bsrcount;
int papercount;
int glasscount;
int biocount;
 
// Berechnen der Termine (auch über das Jahr hinaus, kein Halt am 31.12., auch ein 13. Monat mit unendlich vielen Tagen ist möglich)
if(alba_mode > 1){
 day_y(alba_freq, nrod, albaday, albamonth, q, albacount, albacday, albacmonth, albacfreq, alba_mode - 1);
}else{
 day_x(alba_freq, nrod, albaday, albamonth, q, albacount);
}
if(org_mode > 1){
 day_y(org_freq, nrod, orgday, orgmonth, q, orgcount, orgcday, orgcmonth, orgcfreq, org_mode - 1);
}else{
 day_x(org_freq, nrod, orgday, orgmonth, q, orgcount);
}
if(bsr_mode > 1){
 day_y(bsr_freq, nrod, bsrday, bsrmonth, q, bsrcount, bsrcday, bsrcmonth, bsrcfreq, bsr_mode - 1);
}else{
 day_x(bsr_freq, nrod, bsrday, bsrmonth, q, bsrcount);
}
if(paper_mode > 1){
 day_y(paper_freq, nrod, paperday, papermonth, q, papercount, papercday, papercmonth, papercfreq, paper_mode - 1);
}else{
 day_x(paper_freq, nrod, paperday, papermonth, q, papercount);
}
if(glass_mode > 1){
 day_y(glass_freq, nrod, glassday, glassmonth, q, glasscount, glasscday, glasscmonth, glasscfreq, glass_mode - 1);
}else{
 day_x(glass_freq, nrod, glassday, glassmonth, q, glasscount);
}
if(bio_mode > 1){
 day_y(bio_freq, nrod, bioday, biomonth, q, biocount, biocday, biocmonth, biocfreq, bio_mode - 1);
}else{
 day_x(bio_freq, nrod, bioday, biomonth, q, biocount);
}

// Einlesen der Ersatztermine (durch Feiertagsverschiebung o.ä.)
 ifstream subf;
 subf.open("ersatz.dat");
 int nrosd;                                                          // nr of substitute dates
 subf >> nrosd;
 int*    dayin    = new int[nrosd];
 int*    monthin  = new int[nrosd];
 int*    yearin   = new int[nrosd];
 int*    dayout   = new int[nrosd];
 int*    monthout = new int[nrosd];
 int*    yearout  = new int[nrosd];
 string subtype[nrosd];
 for(int x = 0; x < nrosd; x++){
  subf >> dayin[x] >> monthin[x] >> yearin[x] >> dayout[x] >> monthout[x] >> yearout[x] >> subtype[x];
 }
 subf.close();

// Änderung der durch Ersatztermine betroffenen Abholtermine
char typus = 'Y'; // gelbe Tonnen, gelbe Säcke
subdate(nrod, nrosd, albaday, albamonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus);
typus = 'O'; // orange Wertstofftonne
subdate(nrod, nrosd, orgday, orgmonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus);
typus = 'M'; // Restmülltonne
subdate(nrod, nrosd, bsrday, bsrmonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus);
typus = 'P'; // Altpapiertonne
subdate(nrod, nrosd, paperday, papermonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus);
typus = 'G'; // Altglastonne
subdate(nrod, nrosd, glassday, glassmonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus);
typus = 'B'; // Biotonne
subdate(nrod, nrosd, bioday, biomonth, year, dayin, dayout, monthin, monthout, yearin, yearout, subtype, typus); 

// Streichung der Termine, die über das Jahr hinaus gehen; Herabsetzung der Terminanzahl
cleardate(nrod, albaday, albamonth, albacount);
cleardate(nrod, orgday, orgmonth, orgcount);
cleardate(nrod, bsrday, bsrmonth, bsrcount);
cleardate(nrod, paperday, papermonth, papercount);
cleardate(nrod, glassday, glassmonth, glasscount);
cleardate(nrod, bioday, biomonth, biocount);

zerodate(albacount);
zerodate(orgcount);
zerodate(bsrcount);
zerodate(papercount);
zerodate(glasscount);
zerodate(biocount);

// Ausgabe der Termine nach termine.dat
 outf << "ALBA-Termine: " << albacount << "\n";
 for(int i = 0; i < albacount; i++){
   outf << albaday[i] << "." << albamonth[i] << "." << year << "\n";
 }
 outf << "Orange-Termine: " << orgcount << "\n";
 for(int i = 0; i < orgcount; i++){
  outf << orgday[i] << "." << orgmonth[i] << "." << year << "\n";
 }
 outf << "BSR-Termine: " << bsrcount << "\n";
 for(int i = 0; i < bsrcount; i++){
   outf << bsrday[i] << "." << bsrmonth[i] << "." << year << "\n";
 }
 outf << "Papier-Termine: " << papercount << "\n";
 for(int i = 0; i < papercount; i++){
  outf << paperday[i] << "." << papermonth[i] << "." << year << "\n";
 }
 outf << "Glas-Termine: " << glasscount << "\n";
 for(int i = 0; i < glasscount; i++){
  outf << glassday[i] << "." << glassmonth[i] << "." << year << "\n";
 }
 outf << "Biomüll-Termine: " << biocount << "\n";
 for(int i = 0; i < biocount; i++){
  outf << bioday[i] << "." << biomonth[i] << "." << year << "\n";
 }
 outf << "Weihnachtsbaumabholung: " << nrowd << "\n";
 for(int i = 0; i < nrowd; i++){
  outf << treeday[i] << "." << treemonth[i] << "." << year << "\n";
 }
 outf.close();

// Ausgabe der Termine in die einzelnen Dateien für die weitere Verwendung im Kalender:
 ofstream albf;
 albf.open("alba.dat");
 albf << year << "\n";
 albf << albacount << "\n";
 for(int i = 0; i < albacount; i++){
  albf << albaday[i] << " " << albamonth[i] << "\n";
 }
 albf.close();

 ofstream orgf;
 orgf.open("orange.dat");
 orgf << year << "\n";
 orgf << orgcount << "\n";
 for(int i = 0; i < orgcount; i++){
  orgf << orgday[i] << " " << orgmonth[i] << "\n";
 }
 orgf.close();

 ofstream bsrf;
 bsrf.open("bsr.dat");
 bsrf << year << "\n";
 bsrf << bsrcount << "\n";
 for(int i = 0; i < bsrcount; i++){
  bsrf << bsrday[i] << " " << bsrmonth[i] << "\n";
 }
 bsrf.close();

 ofstream pref;
 pref.open("pappe.dat");
 pref << year << "\n";
 pref << papercount << "\n";
 for(int i = 0; i < papercount; i++){
  pref << paperday[i] << " " << papermonth[i] << "\n";
 }
 pref.close();

 ofstream glaf;
 glaf.open("glass.dat");
 glaf << year << "\n";
 glaf << glasscount << "\n";
 for(int i = 0; i < glasscount; i++){
  glaf << glassday[i] << " " << glassmonth[i] << "\n";
 }
 glaf.close();

 ofstream biof;
 biof.open("bio.dat");
 biof << year << "\n";
 biof << biocount << "\n";
 for(int i = 0; i < biocount; i++){
  biof << bioday[i] << " " << biomonth[i] << "\n";
 }
 biof.close();

 ofstream tref;
 tref.open("baum.dat");
 tref << year << "\n";
 tref << nrowd << "\n";
 for(int i = 0; i < nrowd; i++){
  tref << treeday[i] << " " << treemonth[i] << "\n";
 }
 tref.close();

}

// Funktionen
// Einlesen aus kommentierten Eingabedateien
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

// Berechnen der Abholtermine ohne Änderung
void day_x(int freq, int nrod, int* day, int* month, int schaltjahr, int &count){
if(freq != 0){
 count = 1;
 for(int i = 1; i < nrod; i++){
 if(day[i-1] > (31 - freq) && (month[i-1] == 1 || month[i-1] == 3 || month[i-1] == 5 || month[i-1] == 7 || month[i-1] == 8 || month[i-1] == 10 || month[i-1] == 12)){
  day[i] = day[i-1] - (31 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (30 - freq) && (month[i-1] == 4 || month[i-1] == 6 || month[i-1] == 9 || month[i-1] == 11 )){
  day[i] = day[i-1] - (30 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (29 - freq) && month[i-1] == 2 && schaltjahr == 1){
  day[i] = day[i-1] - (29 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (28 - freq) && month[i-1] == 2 && schaltjahr == 0){
  day[i] = day[i-1] - (28-freq);
  month[i] = month[i-1] + 1;
  count++;
 }else{
  day[i] = day[i-1] + freq;
  month[i] = month[i-1];
  count++;
 }
}
}else{
count = 0;
}
}

// Berechnen der Abholtermine mit Änderungen
void day_y(int freq, int nrod, int* day, int* month, int schaltjahr, int &count, int* cday, int* cmonth, int* cfreq, int cmode){
 int i = 1;
 int j = 0;
 if(freq != 0){
  count = 1;
 }else{
  count = 0;
 }
 start:
if(j <= cmode){
if(freq != 0){
 while(i < nrod){
 if(day[i-1] > (31 - freq) && (month[i-1] == 1 || month[i-1] == 3 || month[i-1] == 5 || month[i-1] == 7 || month[i-1] == 8 || month[i-1] == 10 || month[i-1] == 12)){
  day[i] = day[i-1] - (31 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (30 - freq) && (month[i-1] == 4 || month[i-1] == 6 || month[i-1] == 9 || month[i-1] == 11 )){
  day[i] = day[i-1] - (30 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (29 - freq) && month[i-1] == 2 && schaltjahr == 1){
  day[i] = day[i-1] - (29 - freq);
  month[i] = month[i-1] + 1;
  count++;
 }else if(day[i-1] > (28 - freq) && month[i-1] == 2 && schaltjahr == 0){
  day[i] = day[i-1] - (28-freq);
  month[i] = month[i-1] + 1;
  count++;
 }else{
  day[i] = day[i-1] + freq;
  month[i] = month[i-1];
  count++;
 }
 if(j < cmode){
  if(day[i] >= cday[j] && month[i] >= cmonth[j]){
   day[i] = cday[j];
   month[i] = cmonth[j];
   freq = cfreq[j];
   i++;
   j++;
   goto start;
  }else{
   i++;
  }
 }else{
  i++;
 }
}
}else{
 j++;
 goto start;
}
}
}

// Berechnen, ob ein Jahr ein Schaltjahr ist
int schaltjahr(int year){
 int z;
 if(year%400 == 0){
  z = 1;
 }else if(year%100 == 0){
  z = 0;
 }else if(year%4 == 0){
  z = 1;
 }else{
  z = 0;
 }
 return(z);
}

// Berechnen der Ersatztermine
void subdate(int nrod, int nrosd, int* day, int* month, int year, int* dayin, int* dayout, int* monthin, int* monthout, int* yearin, int* yearout, string* subtype, char typus){
 bool c;
 for(int i = 0; i < nrod; i++){
  for(int j = 0; j < nrosd; j++){
   string str = subtype[j];
   if(day[i] == dayin[j] && month[i] == monthin[j] && year == yearin[j]){
    string::size_type n = str.find(typus, 0);
    c = print(n, str);
    if(c){
     day[i] = dayout[j];
     month[i] = monthout[j];
     break;
    }
   }
  }
 }
}

// Vergleich, ob ein Zeichen (character) in einer Zeichenfolge (string) vorkommt
bool print(string::size_type n, string s){
 bool x;
 if(n == string::npos){
  x = false;
 }else{
  x = true;
 }
 return(x);
}

// Streichen der unmöglichen Termine
void cleardate(int nrod, int* day, int* month, int &count){
 int c = 0;
 if(count < nrod){
  c = 1;
 }
 for(int i = 0; i < nrod; i++){
  if(month[i] > 12 || month[i] < 1 || day[i] < 1 || day[i] > 31){
   day[i] = 0;
   month[i] = 0;
   if(c == 0){
    count--;
   }
  }
 }
}

// Setzen der Abholterminanzahl auf null, wenn negativ
void zerodate(int &count){
 if(count < 0){
  count = 0;
 }
}

