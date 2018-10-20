#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>

using namespace std;

double croot(double x);

int main(int argc, char* argv[]){
 if(argc != 5){
  cerr << "usage: 'cardano A B C D' fuer Gleichungen Ax^3 + Bx^2 + Cx + D = 0\n";
  exit(0);
 }
// Definition der imaginaeren Einheit i
const complex double c_i = _Complex_I;
// Einlesen der Koeffizienten
double coeffA = strtod(argv[1],NULL); // Kubischer Term
double coeffB = strtod(argv[2],NULL); // Quadratischer Term
double coeffC = strtod(argv[3],NULL); // Linearer Term
double coeffD = strtod(argv[4],NULL); // Konstanter Term
// Division durch den kubischen Term
double a = coeffB / coeffA;
double b = coeffC / coeffA;
double c = coeffD / coeffA;
// Bildung der reduzierten Form
double p = b - (1./3.)*pow(a,2.);
double q = (2./27.)*pow(a,2.) - (1./3.)*a*b + c;
// Berechnung der Diskriminante
double D = pow(q/2.,2.)+pow(p/3.,3.);
// Fallunterscheidung
int dcase = 0;
if(D > 0){
 cout << "eine reelle Loesung und zwei komplexe Loesungen\n";
 dcase = 1;
}else if(D == 0){
 cout << "eine doppelte und eine einfache reelle Loesung oder eine dreifache reelle Loesung\n";
 dcase = 2;
}else if(D < 0){
 cout << "drei verschiedene reelle Loesungen\n";
 dcase = 3;
}else{
 exit(2);
}
complex double z_a, z_b, z_c;
// Fall 1 (D > 0)
if(dcase == 1){
double qh = (q/2.);
double sqd = sqrt(D);
double ut = (-qh + sqd);
double vt = (-qh - sqd);
double u = croot(ut);
double v = croot(vt);
// Bestimmung von u und v
 z_a =  u + v;
 z_b = -((u+v)/2)+((u-v)/2)*c_i*sqrt(3);
 z_c = -((u+v)/2)-((u-v)/2)*c_i*sqrt(3);
}
// Fall 2 (D = 0)
if(dcase == 2){
 z_a = 3*q/p;
 z_b = -3*q/(2*p);
 z_c = -3*q/(2*p);
}
// Fall 3 (D < 0)
if(dcase == 3){
 z_a = sqrt(-4.*p/3.)*cos((1./3.)*acos((-1.*q/2.)*sqrt(-27./pow(p,3))));
 z_b = -1*sqrt(-4.*p/3.)*cos((1./3.)*acos((-1.*q/2.)*sqrt(-27./pow(p,3)))+(M_PI/3.));
 z_c = -1*sqrt(-4.*p/3.)*cos((1./3.)*acos((-1.*q/2.)*sqrt(-27./pow(p,3)))-(M_PI/3.));
}
// andere Faelle (eigentlich nicht moeglich)
if(dcase < 1 || dcase > 3){
 exit(3);
}
//Ruecksubstitution
complex double x_a, x_b, x_c;
 x_a = z_a - (a/3.);
 x_b = z_b - (a/3.);
 x_c = z_c - (a/3.);
// Ausgabe
 cout << "Lösungen: \n";
 cout << "x_1 = " << creal(x_a) << " " << cimag(x_a) << "\n";
 cout << "x_2 = " << creal(x_b) << " " << cimag(x_b) << "\n";
 cout << "x_3 = " << creal(x_c) << " " << cimag(x_c) << "\n";
}

double croot(double x){
 double y;
 if (x < 0){
  y = -1*pow((-1.*x),(1./3.));
 }else{
  y = pow(x,(1./3.));
 }
 return(y);
}

