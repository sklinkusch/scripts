# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

using namespace std;
// Funktionen
char* minu(int minute);
void time(long double sun, int &hour, int &minute);
void dst(long double sunrise_moz, long double &sunrise_mez, long double sunset_moz, long double &sunset_mez, long double longitude, long double delta_long, int dsc_an, int dsc_aus, int day, int month);
void check_date(int day, int month, int year, int schaltjahr);
int leapyear(int year);
long double yearday(int day, int month, int year, long double hours);
int dsc_on(int year);
int dsc_off(int year);
int wotag(int d, int m, int j);
int gauss(double zahl);
long double normalize(long double x, long double pi2);
long double normalize_deg(long double x, long double fullcircle);
void calc_hours(int hora, int minuta, long double timedelay, long double &x);

int main(int argc, char* argv[]){
 if(argc != 7 && argc != 10){
  cout << "Gebrauch: ./sun3 1 <Tag> <Monat> <Jahr> <Stunde> <Minute> <Breite> <Länge> <Zeitverschiebung zu UTC> für Jahre ab 2000\n";
  cout << "oder: ./sun3 0 <Tag> <Monat> <Jahr> <Stunde> <Minute> für Berlin-Spandau\n";
  exit(1);
 }
 // Konstanten
 long double pi = (long double) 4.* atanl((long double) 1.);
 long double pi2 = (long double) 2.*pi;
 long double halfcircle = (long double) 180.;
 long double fullcircle = (long double) 360.;
 long double delta_long = 15.;
 // Einlesen der Parameter
 int mode = atoi(argv[1]);
 int day = atoi(argv[2]);
 int month = atoi(argv[3]);
 int year = atoi(argv[4]);
 int hora = atoi(argv[5]);
 int minuta = atoi(argv[6]);
 long double latitude, longitude;
 if(mode == 1){
  latitude = strtod(argv[7],NULL);
  longitude = strtod(argv[8],NULL);
 }else{
  latitude = 52.51483;
  longitude = 13.17990;
 }
 long double B = pi*latitude/halfcircle; // geographische Breite im Bogenmass
 int schaltjahr = leapyear(year);
// Formeln für den Sonnenlauf
 check_date(day, month, year, schaltjahr);
// Sommerzeit
 int dsc_an = dsc_on(year);
 int dsc_aus = dsc_off(year);
 long double timedelay = (long double) 0.;
 if(mode == 1){
  timedelay = strtod(argv[9],NULL);
 }else{
  if(month > 3 && month < 10){
   timedelay = (long double) 2.;
  }else if(month < 3 || month > 10){
   timedelay = (long double) 1.;
  }else if((month == 3 && day < dsc_an) || (month == 10 && day >= dsc_aus)){
   timedelay = (long double) 1.;
  }else if((month == 3 && day >= dsc_an) || (month == 10 && day < dsc_aus)){
   timedelay = (long double) 2.;
  }
 }
 long double hours = 0.;
 calc_hours(hora, minuta, timedelay, hours);
 long double T = yearday(day, month, year, hours); // Bestimmung der Anzahl der Tage seit dem 1.1.2000
 long double L = (long double) 280.460 + (long double) 0.9856474*T; // mittlere ekliptikale Länge
 long double Ln = normalize_deg(L, fullcircle);
 long double g = (long double) 357.528 + (long double) 0.9856003*T; // Anomalie
 long double gn = normalize_deg(g, fullcircle);
 long double e = pi2*((long double) 23.43929111 - 0.0000004*T)/fullcircle; // Neigung der Erdachse
 long double Lambda = Ln + 1.915 * sinl(pi*gn/halfcircle) + (0.020 * sinl(pi2*gn/halfcircle)); // ekliptikale Länge der Sonne
 long double DK = asinl(sinl(e)*sinl(pi*Lambda/halfcircle)); // Deklination
 long double RA = atanl(tanl(pi*Lambda/halfcircle)*cosl(e)); // Rektaszension
// Korrekturen für die Rektaszension (arcustangens nicht eindeutig)
 if(RA < (long double) 0.){
  RA += pi;
 }
 if(RA > pi){
  RA -= pi;
 }
// Horizontalkoordinaten der Sonne
long double T_midnight = yearday(day, month, year, 0);
long double t_midnight = ((T_midnight)/36525);
long double t_mstar = 6.697376 + 2400.05134*t_midnight + 1.002738*hours;
long double theta_g = delta_long*t_mstar;
long double theta = theta_g + longitude;
long double theta_rad = pi*theta/halfcircle;
long double tau = theta_rad - RA;
long double nenn = (cosl(tau)*sinl(B)-(tanl(DK)*cosl(B)));
long double azimut = atanl(sinl(tau)/nenn);
 if(nenn < 0){
  azimut += pi;
 }
 if(azimut > pi){
  azimut -= pi2;
 }
 if(azimut < (long double) -1.*pi){
  azimut += pi2;
 }
long double azimut_deg = halfcircle*azimut/pi;
long double height = asinl(cosl(DK)*cos(tau)*cosl(B) + (sinl(DK)*sinl(B)));
long double height_deg = halfcircle*height/pi;
// Refraktionskorrektur
long double R_korr = 1.02 / tanl(pi*(height_deg+(10.3/(height_deg+5.11)))/halfcircle);
long double h_korr = height_deg + R_korr/60.;
char* dies = new char[2];
char* mens = new char[2];
char* minute = new char[2];
dies = minu(day);
mens = minu(month);
minute = minu(minuta);
int lat_grad, long_grad, lat_min, long_min, lat_sec, long_sec;
lat_grad = (int) latitude;
lat_min = (int) (60.*(latitude - (long double) lat_grad));
lat_sec = (int) (3600.*(latitude - (long double) lat_grad - (long double) lat_min/60.));
long_grad = (int) longitude;
long_min = (int) (60.*(longitude - (long double) long_grad));
long_sec = (int) (3600.*(longitude - (long double) long_grad - (long double) long_min/60.));
// Ausgabe
 cout << "Daten für den " << dies << "." << mens << "." << year << " um " << hora << ":" << minute << " Ortszeit:\n"; 
 cout << "Geographische Breite: " << lat_grad << " Grad " << lat_min << " Minuten " << lat_sec << " Sekunden\n";
 cout << "Geographische Länge: " << long_grad << " Grad " << long_min << " Minuten " << long_sec << " Sekunden\n";
 cout << "Zeitverschiebung: " << timedelay << " Stunden\n";
 cout << "Azimut: " << azimut_deg << " Grad\n";
 cout << "Höhe (bei 1010 mbar und 10 Grad Celsius): " << h_korr << " Grad\n"; 
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

long double yearday(int day, int month, int year, long double hours){
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
 x += (hours - 12.)/24.;
//  x /= (long double) 36525.;
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

long double normalize_deg(long double x, long double fullcircle){
 while(x < (long double) 0.){
  x = x + fullcircle;
 }
 while(x > fullcircle){
  x = x - fullcircle;
 }
return(x);
}

void calc_hours(int hora, int minuta, long double timedelay, long double &x){
 x = (long double) hora + ((long double) minuta/(long double) 60.);
 x -= timedelay;
/* if(x < 0.){
  x += 24.;
  T -= 1;
 }else if(x > 24.){
  x -= 24.;
  T += 1;
 }
*/
}

