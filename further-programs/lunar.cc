# include <iostream>
# include <fstream>
# include <math.h>
# include <string.h>
# include <stdio.h>
# include <stdlib.h>

using namespace std;

// Programm zur Berechnung der Mondphasen für Daten nach 1999

//Funktionen
double lunarday(int day, int month, int year, int hour, int minute);
int leapyear(int year);

int main(int argc, char* argv[]){
 if(argc != 6){
  cerr << "Usage: lunar <day> <month> <year> <hour> <minute>\n";
  exit(1);
 }
 int day, month, year, hour, minute; // Datumsangaben in Tagen, Minuten, Jahren, Stunden, Minuten
 day    = atoi(argv[1]);
 month  = atoi(argv[2]);
 year   = atoi(argv[3]);
 hour   = atoi(argv[4]);
 minute = atoi(argv[5]);
 double lday = lunarday(day, month, year, hour, minute); // Datumsangabe in Dezimaltagen seit 01.01.1999
 cout << "Tage seit 1.1.1999 3:49: " << lday << "\n";
 const double startvollmond = 1.159027778;
 const double lunarmonth = 29.530588;
 double diffday = (lday - startvollmond);
 double moonmonths = diffday / lunarmonth;
 int monthmoon = (int) moonmonths;
 cout << "Zahl der Mondumläufe seit 1.1.1999: " << monthmoon << "\n";
 double moonphase = moonmonths - (double) monthmoon;
 cout << "Mondphase: " << moonphase << "\n";
 cout << "Mondphase = 0.00: Vollmond\n";
 cout << "Mondphase = 0.25: abnehmender Mond\n";
 cout << "Mondphase = 0.50: Neumond\n";
 cout << "Mondphase = 0.75: zunehmender Mond\n";
 double s;
 if(moonphase <= 0.5){
  s = cos(2*moonphase*M_PI);
 }else{
  s = -1.*cos(2*moonphase*M_PI);
 }
 cout << "Parameter s = " << s << "\n";
}

double lunarday(int day, int month, int year, int hour, int minute){
 double lday = 0;
 const int startyear = 1999;
 if(year >= startyear){
  int i = startyear;
  while(i < year){
   int j = leapyear(i);
   if(j == 1){
    lday += 366.;
   }else{
    lday += 365.;
   }
   i++;
  }
 }
 int j = leapyear(year);
 if(j == 1 && month > 2){
  lday += 1.;
 }else{
  lday += 0.;
 }
 switch(month)
 {
 case 1:	lday += 0.;
		break;
 case 2:        lday += 31.;
		break;
 case 3:	lday += 59.;
		break;
 case 4:	lday += 90.;
		break;
 case 5:	lday += 120.;
		break;
 case 6:	lday += 151.;
		break;
 case 7:	lday += 181.;
		break;
 case 8:	lday += 212.;
		break;
 case 9:	lday += 243.;
		break;
 case 10:	lday += 273.;
		break;
 case 11:	lday += 304.;
		break;
 case 12:	lday += 334.;
		break;
 }
 lday += ((double) day - 1.);
 lday += ((double) hour / 24.);
 lday += ((double) minute / (24.*60.));
 return(lday);
}

int leapyear(int year){
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
 
 
