# include<iostream>
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<string.h>

using namespace std;

// Berechnung des Taupunkts

int main(int argc, char* argv[]){
 if(argc != 4){
  cerr << "Format: ./dewpoint <mode> <Temperatur in °C> <relative Luftfeuchte in %>\n";
  cerr << "<mode> = 1: ueber Wasser, <mode> = 0: ueber Eis\n";
  exit(1);
 }

int mode; // mode = 1: über Wasser, mode = 0: über Eis;
mode = atoi(argv[1]);
double theta = strtod(argv[2],NULL);
double varphi = strtod(argv[3],NULL);
double phi = varphi/100.;
double K_zwei, K_drei;
if(mode == 1){
K_zwei = 17.62; //
K_drei = 243.12; // °C
}else{
K_zwei = 22.46;
K_drei = 272.62; // °C
}
if(theta < -273.15){
 cerr << "Berechnung nicht moeglich, weil Temperatur unter -273.15°C\n";
 exit(2);
}
if(phi < 0. || phi > 1.){
 cerr << "Berechnung nicht moeglich, weil relative Luftfeuchte ausserhalb des zugelassenen Wertebereichs\n";
 exit(3);
}
double zahler = ((K_zwei * theta)/(K_drei + theta)) + log(phi);
double nenner = ((K_zwei * K_drei)/(K_drei + theta)) - log(phi);
double dewpoint = K_drei * zahler / nenner;
cout << "Taupunkt: " << dewpoint << " °C\n";
}
