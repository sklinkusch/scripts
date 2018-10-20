# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

using namespace std;
// Funktionen
char* minu(int minute);
void time(long double sun, int &hour, int &minute);
void mez(long double &rmez, long double rmoz, long double &smez, long double smoz, long double longitude, long double delta_long, long double zone);
void woz(long double &rwoz, long double &swoz, long double delta_t);
void moz(long double rwoz, long double swoz, long double &rmoz, long double &smoz, long double wozmoz);
void dst(long double sunrise_moz, long double &sunrise_mez, long double sunset_moz, long double &sunset_mez, long double longitude, long double delta_long, int dsc_an, int dsc_aus, int day, int month);
void check_date(int day, int month, int year, int schaltjahr);
int leapyear(int year);
long double yearday(int day, int month, int year);
int dsc_on(int year);
int dsc_off(int year);
int wotag(int d, int m, int j);
int gauss(double zahl);
long double normalize(long double x, long double pi2);

int main(int argc, char* argv[]){
 if(argc != 4){
  cout << "Gebrauch: ./sun <Tag> <Monat> <Jahr> für Jahre ab 2000\n";
  exit(1);
 }
 // Konstanten
 long double pi = (long double) 4.* atanl((long double) 1.);
 long double pi2 = (long double) 2.*pi;
 long double halfcircle = (long double) 180.;
 long double fullcircle = (long double) 360.;
 long double latitude = 52.51483;
 long double longitude = 13.17990;
 long double delta_long = 15.;
 long double B = pi*latitude/halfcircle; // geographische Breite im Bogenmass
 long double h_sun = ((long double) -50./(long double) 60.)*pi/halfcircle; // - 50 Bogenminuten zur Berechnung des Sonnenauf-/untergangs
 long double h_dawnc = ((long double) -6.)*pi/halfcircle; // - 6° zur Berechnung der bürgerlichen Dämmerung
 long double h_dawnn = ((long double) -12.)*pi/halfcircle; // - 12° zur Berechnung der nautischen Dämmerung
 long double h_dawna = ((long double) -18.)*pi/halfcircle; // - 18° zur Berechnung der astronomischen Dämmerung
 // Einlesen der Parameter
 int day = atoi(argv[1]);
 int month = atoi(argv[2]);
 int year = atoi(argv[3]);
 int schaltjahr = leapyear(year);
// Formeln für den Sonnenlauf
 check_date(day, month, year, schaltjahr);
 long double T = yearday(day, month, year); // Bestimmung der Anzahl der Tage seit dem 1.1.2000
 long double M = pi2*((long double) 0.993133 + ((long double) 99.997361*T));
 M = normalize(M, pi2);
 long double L = pi2*((long double) 0.7859453 + M/pi2 + ((long double) 6893.*sinl(M) + (long double) 72.*sinl((long double) 2.*M) + (long double) 6191.2*T)/(long double) 1296000); // zurückgelegter Weg der Sonne seit Frühlingsanfang
 L = normalize(L, pi2);
 long double e = pi2*((long double) 23.43929111 + ((long double) -46.8150*T - (long double) 0.00059*T*T + (long double) 0.001813*T*T*T)/(long double) 3600.)/fullcircle; // Neigung der Erdachse
 long double DK = asinl(sinl(e)*sinl(L)); // Deklination
 long double RA = atanl(tanl(L)*cosl(e)); // Rektaszension
// Korrekturen für die Rektaszension (arcustangens nicht eindeutig)
 if(RA < (long double) 0.){
  RA += pi;
 }
 if(L > pi){
  RA += pi;
 }
 long double RAh = (long double) 24.*RA/pi2;
// mittlere Rektaszension
 long double RAm = (long double) 18.71506921 + (long double) 2400.0513369*T + ((long double) 0.000025862 - (long double) 0.00000000172*T)*T*T;
int quotient = (int) (RAm/(long double) 24.);
 RAm = RAm - (long double) quotient * (long double) 24.;
// Zeitgleichung
long double woz_moz = 1.0027379*(RAm - RAh);
// Deklination der Sonne
// double d = 0.40954*sin(0.0172*(T-79.35));
// Zeitdifferenzen
 long double delta_trs = (long double) 12.*acosl((sinl(h_sun) - (sinl(B)*sinl(DK)))/(cosl(B)*cosl(DK)))/pi; // Zeitdifferenz für Sonnenauf-/untergang
 long double delta_ttc = (long double) 12.*acosl((sinl(h_dawnc) - (sinl(B)*sinl(DK)))/(cosl(B)*cosl(DK)))/pi; // Zeitdifferenz für bürgerliche Dämmerung
 long double delta_ttn = (long double) 12.*acosl((sinl(h_dawnn) - (sinl(B)*sinl(DK)))/(cosl(B)*cosl(DK)))/pi; // Zeitdifferenz für nautische Dämmerung
 long double delta_tta = (long double) 12.*acosl((sinl(h_dawna) - (sinl(B)*sinl(DK)))/(cosl(B)*cosl(DK)))/pi; // Zeitdifferenz für astronomische Dämmerung
 long double sunrise_woz, sunset_woz, dawnmc_woz, dawnec_woz, dawnmn_woz, dawnen_woz, dawnma_woz, dawnea_woz;
 woz(sunrise_woz, sunset_woz, delta_trs);
 woz(dawnmc_woz, dawnec_woz, delta_ttc);
 woz(dawnmn_woz, dawnen_woz, delta_ttn);
 woz(dawnma_woz, dawnea_woz, delta_tta);
 long double sunrise_moz, sunset_moz, dawnmc_moz, dawnec_moz, dawnmn_moz, dawnen_moz, dawnma_moz, dawnea_moz;
 moz(sunrise_woz, sunset_woz, sunrise_moz, sunset_moz, woz_moz);
 moz(dawnmc_woz, dawnec_woz, dawnmc_moz, dawnec_moz, woz_moz);
 moz(dawnmn_woz, dawnen_woz, dawnmn_moz, dawnen_moz, woz_moz);
 moz(dawnma_woz, dawnea_woz, dawnma_moz, dawnea_moz, woz_moz);
 int dsc_an = dsc_on(year);
 int dsc_aus = dsc_off(year);
 long double sunrise_mez, sunset_mez, dawnmc_mez, dawnec_mez, dawnmn_mez, dawnen_mez, dawnma_mez, dawnea_mez;
 dst(sunrise_moz, sunrise_mez, sunset_moz, sunset_mez, longitude, delta_long, dsc_an, dsc_aus, day, month);
 dst(dawnmc_moz, dawnmc_mez, dawnec_moz, dawnec_mez, longitude, delta_long, dsc_an, dsc_aus, day, month);
 dst(dawnmn_moz, dawnmn_mez, dawnen_moz, dawnen_mez, longitude, delta_long, dsc_an, dsc_aus, day, month);
 dst(dawnma_moz, dawnma_mez, dawnea_moz, dawnea_mez, longitude, delta_long, dsc_an, dsc_aus, day, month);
 if(dawnma_mez != dawnma_mez){
  dawnma_mez = 0.;
 }
 if(dawnea_mez != dawnea_mez){
  dawnea_mez = 0.;
 }
 int hour_rise, minute_rise, hour_set, minute_set, hour_dmc, minute_dmc, hour_dec, minute_dec, hour_dmn, minute_dmn, hour_den, minute_den, hour_dma, minute_dma, hour_dea, minute_dea;
 char* hr_rise = new char[3];
 char* hr_set = new char[3];
 char* hr_dmc = new char[3];
 char* hr_dec = new char[3];
 char* hr_dmn = new char[3];
 char* hr_den = new char[3];
 char* hr_dma = new char[3];
 char* hr_dea = new char[3];
 char* min_rise = new char[3];
 char* min_set = new char[3];
 char* min_dmc = new char[3];
 char* min_dec = new char[3];
 char* min_dmn = new char[3];
 char* min_den = new char[3];
 char* min_dma = new char[3];
 char* min_dea = new char[3];
 time(sunrise_mez, hour_rise, minute_rise);
 time(sunset_mez, hour_set, minute_set);
 time(dawnmc_mez, hour_dmc, minute_dmc);
 time(dawnec_mez, hour_dec, minute_dec);
 time(dawnmn_mez, hour_dmn, minute_dmn);
 time(dawnen_mez, hour_den, minute_den);
 time(dawnma_mez, hour_dma, minute_dma);
 time(dawnea_mez, hour_dea, minute_dea);
 if(hour_dea >= 24){
  hour_dea -= 24;
 }
 hr_rise = minu(hour_rise);
 hr_set = minu(hour_set);
 hr_dmc = minu(hour_dmc);
 hr_dec = minu(hour_dec);
 hr_dmn = minu(hour_dmn);
 hr_den = minu(hour_den);
 hr_dma = minu(hour_dma);
 hr_dea = minu(hour_dea);
 min_rise = minu(minute_rise);
 min_set = minu(minute_set);
 min_dmc = minu(minute_dmc);
 min_dec = minu(minute_dec);
 min_dmn = minu(minute_dmn);
 min_den = minu(minute_den);
 min_dma = minu(minute_dma);
 min_dea = minu(minute_dea);
 if(strcmp(hr_dma,hr_dea) == 0 && strcmp(min_dma,min_dea) == 0){
  sprintf(hr_dma, "--");
  sprintf(hr_dea, "--");
  sprintf(min_dma, "--");
  sprintf(min_dea, "--");
 }
 cout << "Daten für den " << day << "." << month << "." << year << ":\n"; 
 cout << "astronomische Morgendämmerung ab " << hr_dma << ":" << min_dma << " Uhr\n";
 cout << "nautische Morgendämmerung ab " << hr_dmn << ":" << min_dmn << " Uhr\n";
 cout << "bürgerliche Morgendämmerung ab " << hr_dmc << ":" << min_dmc << " Uhr\n";
 cout << "Sonnenaufgang um " << hr_rise << ":" << min_rise << " Uhr\n";
 cout << "Sonnenuntergang um " << hr_set << ":" << min_set << " Uhr\n";
 cout << "bürgerliche Abenddämmerung bis " << hr_dec << ":" << min_dec << " Uhr\n";
 cout << "nautische Abenddämmerung bis " << hr_den << ":" << min_den << " Uhr\n";
 cout << "astronomische Abenddämmerung bis " << hr_dea << ":" << min_dea << " Uhr\n";
}

char* minu(int minute){
 char* x = new char[3];
 if(minute < 10){
  sprintf(x, "0%d", minute);
 }else{
  sprintf(x, "%d", minute);
 }
 return(x);
}

void time(long double sun, int &hour, int &minute){
 hour = (int) sun;
 minute = (int) (60.*(sun - (long double) hour));
}

void mez(long double &rmez, long double rmoz, long double &smez, long double smoz, long double longitude, long double delta_long, long double zone){
rmez = rmoz - (longitude/delta_long) + zone;
smez = smoz - (longitude/delta_long) + zone;
}

void woz(long double &rwoz, long double &swoz, long double delta_t){
 rwoz = (long double) 12. - delta_t;
 swoz = (long double) 12. + delta_t;
}

void moz(long double rwoz, long double swoz, long double &rmoz, long double &smoz, long double wozmoz){
 rmoz = rwoz + wozmoz;
 smoz = swoz + wozmoz;
}

void dst(long double sunrise_moz, long double &sunrise_mez, long double sunset_moz, long double &sunset_mez, long double longitude, long double delta_long, int dsc_an, int dsc_aus, int day, int month){
 if(month < 3 || month > 10){
  mez(sunrise_mez, sunrise_moz, sunset_mez, sunset_moz, longitude, delta_long, (long double) 1.);
 }else if(month > 3 && month < 10){
  mez(sunrise_mez, sunrise_moz, sunset_mez, sunset_moz, longitude, delta_long, (long double) 2.);
 }else if((month == 3 && day < dsc_an) || (month == 10 && day >= dsc_aus)){
  mez(sunrise_mez, sunrise_moz, sunset_mez, sunset_moz, longitude, delta_long, (long double) 1.);
 }else if((month == 3 && day >= dsc_an) || (month == 10 && day < dsc_aus)){
  mez(sunrise_mez, sunrise_moz, sunset_mez, sunset_moz, longitude, delta_long, (long double) 2.);
 }
}

void check_date(int day, int month, int year, int schaltjahr){
 if(day < 1 || day > 31){
  cout << "Datum nicht existent!\n"; exit(2);
 }
 if(day > 30 && (month == 4 || month == 6 || month == 9 || month == 11)){
  cout << "Datum nicht existent!\n"; exit(2);
 }
 if(day > 29 && month == 2 && schaltjahr == 1){
  cout << "Datum nicht existent!\n"; exit(2);
 }
 if(day > 28 && month == 2 && schaltjahr == 0){
  cout << "Datum nicht existent!\n"; exit(2);
 }
 if(month < 1 || month > 12){
  cout << "Datum nicht existent!\n"; exit(2);
 }
 if(year < 2000){
  cout << "Datum nicht zulaessig!\n"; exit(2);
 }
}

int leapyear(int year){
 int x = 0;
 if(year%400 == 0){
  x = 1;
 }else if(year%100 == 0){
  x = 0;
 }else if(year%4 == 0){
  x = 1;
 }else{
  x = 0;
 }
 return(x);
}

long double yearday(int day, int month, int year){
 long double x = 0.;
 int z = 0;
 for(int aktyear = 2000; aktyear < year; aktyear++){
  z = leapyear(aktyear);
  if(z == 1){
   x += (long double) 366.;
  }else{
   x += (long double) 365.;
  }
 }
 int schaltjahr = leapyear(year);
 if(schaltjahr == 1){
  switch(month){
   case 1:     x += (long double) day; break;
   case 2:     x += (long double) day + (long double) 31.; break;
   case 3:     x += (long double) day + (long double) 60.; break;
   case 4:     x += (long double) day + (long double) 91.; break;
   case 5:     x += (long double) day + (long double) 121.; break;
   case 6:     x += (long double) day + (long double) 152.; break;
   case 7:     x += (long double) day + (long double) 182.; break;
   case 8:     x += (long double) day + (long double) 213.; break;
   case 9:     x += (long double) day + (long double) 244.; break;
   case 10:    x += (long double) day + (long double) 274.; break;
   case 11:    x += (long double) day + (long double) 305.; break;
   case 12:    x += (long double) day + (long double) 335.; break;
   default:    cout << "Falscher Monat angegeben!\n"; exit(2);
  }
 }else{
  switch(month){
   case 1:     x += (long double) day; break;
   case 2:     x += (long double) day + (long double) 31.; break;
   case 3:     x += (long double) day + (long double) 59.; break;
   case 4:     x += (long double) day + (long double) 90.; break;
   case 5:     x += (long double) day + (long double) 120.; break;
   case 6:     x += (long double) day + (long double) 151.; break;
   case 7:     x += (long double) day + (long double) 181.; break;
   case 8:     x += (long double) day + (long double) 212.; break;
   case 9:     x += (long double) day + (long double) 243.; break;
   case 10:    x += (long double) day + (long double) 273.; break;
   case 11:    x += (long double) day + (long double) 304.; break;
   case 12:    x += (long double) day + (long double) 334.; break;
   default:    cout << "Falscher Monat angegeben!\n"; exit(2);
  }
 }
 x -= (long double) 1.;
 x /= (long double) 36525.;
 return(x);
}

int dsc_on(int year){
 int result = 31 - wotag(31,3,year);
 return(result);
}

int dsc_off(int year){
 int result = 31 - wotag(31,10,year);
 return(result);
}

int wotag(int d, int m, int j){
 int result;
 int interim;
 int c = gauss((double) j/100.);
 int y = j - 100*c;
 int m_neu;
 if(m == 1){
  m_neu = 11;
  if(y == 0){
   y = 99;
   c = c - 1;
  }else{
   y = y - 1;
  }
 }else if(m == 2){
  m_neu = 12;
   if(y == 0){
    y = 99;
    c = c - 1;
   }else{
    y = y - 1;
   }
  }else{
   m_neu = m - 2;
  }
 interim = d + gauss(2.6*(double) m_neu - 0.2) + y + gauss((double) y/4.) + gauss((double) c/4.) - 2*c;
 result = interim%7;
 while(result < 0){
  result = result + 7;
 }
return(result);
}

int gauss(double zahl){
 int result;
 if(zahl > 0){
  result = (int) zahl;
 }else if(zahl < 0){
  result = (int) zahl - 1;
 }else{
  result = 0;
 }
 return(result);
}

long double normalize(long double x, long double pi2){
 while(x < (long double) 0.){
  x = x + pi2;
 }
 while(x > pi2){
  x = x - pi2;
 }
return(x);
}

