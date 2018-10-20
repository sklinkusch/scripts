#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Programm zur Erstellung von Kalendervorlagen

using namespace std;

// Funktionen:
int gauss(double zahl);                                                                                                               // Routine zur Bestimmung des Ostersonntagsdatums
int wotag(int day, int month, int year);                                                                                              // Routine zum Bestimmen des Wochentags
void ostermontag(int year, int *tag_om, int *mon_om, int *tag_os, int *mon_os);                                                       // Berechnen des Ostermontagsdatums aus Ostersonntagsdatum
void karfreitag(int tag_os, int mon_os, int *tag_kf, int *mon_kf);                                                                    // Berechnen des Karfreitagsdatums aus Ostersonntagsdatum
void himmelfahrt(int tag_os, int mon_os, int *tag_hf, int *mon_hf);                                                                   // Berechnen des Himmelfahrtsdatums aus Ostersonntagsdatum
void pfingsten(int tag_hf, int mon_hf, int *tag_pm, int *mon_pm);                                                                     // Berechnen des Pfingstmontagsdatums aus Himmelfahrtsdatum
void parsedate(int count, int i, int e, int* day, int* month, int aktmonth, int &z, ofstream& outf, char* inn, char* c, char* off);   // Prüfung auf Übereinstimmung und Ausgabe
void read_data(int count, int* day, int* month, ifstream& inf);                                                                       // Routine zum Einlesen der zu markierenden Daten
int leapyear(int year);                                                                                                               // Routine zum Feststellen eines Schaltjahres

int main(int argc, char* argv[]){
if (argc != 3){
 cerr << "Bitte Kalenderjahr und Ausgabedatei(tex) angeben!\n";
 exit(1);
}
int year = atoi(argv[1]);
ofstream outf;
outf.open(argv[2]);

ifstream tref; // Weihnachtsbaum-Abholung
ifstream pref; // Altpapiertonne
ifstream bsrf; // Mülltonne
ifstream albf; // Abholung der gelben Säcke
ifstream orgf; // orange Wertstofftonne
ifstream glaf; // Altglastonne
ifstream biof; // Bioguttonne
tref.open("baum.dat"); pref.open("pappe.dat"); bsrf.open("bsr.dat"); albf.open("alba.dat");
orgf.open("orange.dat"); glaf.open("glass.dat"); biof.open("bio.dat");
int baumyear, papyear, bsryear, albayear, orgyear, glasyear, bioyear; // Jahreszahlen für die einzelnen Leerungsdaten
tref >> baumyear;
pref >> papyear;
bsrf >> bsryear;
albf >> albayear;
orgf >> orgyear;
glaf >> glasyear;
biof >> bioyear;
if(baumyear != year || papyear != year || bsryear != year || albayear != year || orgyear != year || glasyear != year || bioyear != year){
 cerr << "Jahr der Abholungstermine und Kalenderjahr stimmen nicht überein!\n";
 exit(2);
}
int wbcount, papcount, bsrcount, albacount, orgcount, glascount, biocount; // Zahl der Leerungen in diesem Jahr
tref >> wbcount;
pref >> papcount;
bsrf >> bsrcount;
albf >> albacount;
orgf >> orgcount;
glaf >> glascount;
biof >> biocount;
int* baumday = new int[wbcount]; int* baummon = new int[wbcount];      // Daten der Weihnachtsbaumabholung
int* papday = new int[papcount]; int* papmon = new int[papcount];      // Leerung der Altpapiertonne
int* bsrday = new int[bsrcount]; int* bsrmon = new int[bsrcount];      // Leerung der Restmülltonne
int* albaday = new int[albacount]; int* albamon = new int[albacount];  // Leerung der gelben Wertstofftonne/Abholung des gelben Sackes
int* orgday = new int[orgcount]; int* orgmon = new int[orgcount];      // Leerung der orangen Wertstofftonne
int* glasday = new int[glascount]; int* glasmon = new int[glascount];  // Leerung der Altglastonne
int* bioday = new int[biocount]; int* biomon = new int[biocount];      // Leerung der Bioguttonne
read_data(wbcount, baumday, baummon, tref);
read_data(papcount, papday, papmon, pref);
read_data(bsrcount, bsrday, bsrmon, bsrf);
read_data(albacount, albaday, albamon, albf);
read_data(orgcount, orgday, orgmon, orgf);
read_data(glascount, glasday, glasmon, glaf);
read_data(biocount, bioday, biomon, biof);
int tag_os = 0; int mon_os = 0; //Ostersonntag
int tag_om = 0; int mon_om = 0; //Ostermontag
int tag_kf = 0; int mon_kf = 0; //Karfreitag
int tag_hf = 0; int mon_hf = 0; //Christi Himmelfahrt
int tag_pm = 0; int mon_pm = 0; //Pfingstmontag
ostermontag(year, &tag_om, &mon_om, &tag_os, &mon_os);
karfreitag(tag_os, mon_os, &tag_kf, &mon_kf);
himmelfahrt(tag_os, mon_os, &tag_hf, &mon_hf);
pfingsten(tag_hf, mon_hf, &tag_pm, &mon_pm);
outf << "\\documentclass[10pt,a4paper,landscape]{article}\n";                              // LaTeX-Präambel: Dokumentklasse, Format A4 quer, Schriftgröße 10 pt
outf << "\\usepackage[dvips,left=5.0cm,right=1.0cm,top=10.5cm,bottom=1.0cm]{geometry}\n";  // Festlegung der Seitenränder (Kalender wird nach rechts unten gedruckt), Ausgabetreiber dvips
outf << "\\usepackage[utf8]{inputenc}\n";                                                  // Festlegung der Zeichencodierung auf Unicode UTF-8
outf << "\\usepackage[ngerman]{babel}\n";                                                  // Festlegung der Silbentrennung auf neue deutsche Rechtschreibung
outf << "\\usepackage[usenames,dvipsnames]{color}\n";                                      // Einbeziehung von Farben (sowohl Farbpalette 'usenames' als auch 'dvipsnames')
outf << "\\usepackage{colortbl}\n";                                                        // Farben in Tabellen (muss zwingend NACH color angegeben werden)
outf << "\\pagestyle{empty}\n";                                                            // Seite ohne Kopfzeile, Fußzeile, Seitennummer
outf << "\\newcommand{\\bb}[1]{\\cellcolor{MidnightBlue}\\textcolor{white}{\\bf #1}}\n";   // \bb als Kurzbefehl für weiße Fettschrift auf blauem Hintergrund
outf << "\\newcommand{\\gb}[1]{\\cellcolor{SpringGreen}\\textcolor{black}{\\bf #1}}\n";    // \gb als Kurzbefehl für schwarze Fettschrift auf hellgrünem Hintergrund
outf << "\\newcommand{\\db}[1]{\\cellcolor{ForestGreen}\\textcolor{white}{\\bf #1}}\n";    // \db als Kurzbefehl für weiße Fettschrift auf dunkelgrünem Hintergrund
outf << "\\newcommand{\\yb}[1]{\\cellcolor{yellow}\\textcolor{black}{\\bf #1}}\n";         // \yb als Kurzbefehl für schwarze Fettschrift auf gelbem Hintergrund
outf << "\\newcommand{\\ob}[1]{\\cellcolor{Orange}\\textcolor{white}{\\bf #1}}\n";         // \ob als Kurzbefehl für weiße Fettschrift auf orangem Hintergrund
outf << "\\newcommand{\\iv}[1]{\\cellcolor{black}\\textcolor{white}{\\bf #1}}\n";          // \iv als Kurzbefehl für weiße Fettschrift auf schwarzem Hintergrund
outf << "\\newcommand{\\bv}[1]{\\cellcolor{Brown}\\textcolor{white}{\\bf #1}}\n";          // \bv als Kurzbefehl für weiße Fettschrift auf braunem Hintergrund
outf << "\\newcommand{\\rb}[1]{\\textbf{\\textcolor{red}{#1}}}\n";                         // \rb als Kurzbefehl für rote Fettschrift
outf << "\\newcommand{\\hv}[1]{\\textbf{\\textcolor{Gray}{#1}}}\n\n";                      // \hv als Kurzbefehl für graue Fettschrift
outf << "\\begin{document}\n";                                                             // Beginn des Dokuments
outf << "\\begin{tabular}{|ccccccc|ccccccc|ccccccc|ccccccc|}\n";                           // Beginn einer Tabelle: |: vertikale Linie, c: zentrierte Spalte
outf << "\\hline\n";                                                                       // horizontale Linie
outf << "\\multicolumn{28}{|c|}{\\textbf{Entsorgungskalender " << year << "}} \\\\\\hline\n";
outf << "\\multicolumn{7}{|c|}{\\bf Januar} & \\multicolumn{7}{|c|}{\\bf Februar} & \\multicolumn{7}{|c|}{\\bf M\"arz} & \\multicolumn{7}{|c|}{\\bf April} \\\\\n";
char* bb = new char[16]; char* gb = new char[16]; char* db = new char[16]; char* yb = new char[16]; char* ob = new char[16];
char* iv = new char[16]; char* bv = new char[16]; char* oy = new char[16]; // Definition der einzelnen Farbkategorien
sprintf(bb, "\\bb{"); sprintf(gb, "\\gb{"); sprintf(db, "\\db{"); sprintf(yb, "\\yb{"); sprintf(ob, "\\ob{");
sprintf(iv, "\\iv{"); sprintf(bv, "\\bv{"); sprintf(oy, "} & ");
int schaltjahr = leapyear(year);
int a = wotag(1,1,year);
int b = wotag(1,2,year);
int c = wotag(1,3,year);
int d = wotag(1,4,year);
int az = 31;
int bz;
if(schaltjahr == 1){
 bz = 29;
}else{
 bz = 28;
}
int cz = 31;
int dz = 30;
int ae = -a;
int be = -b;
int ce = -c;
int de = -d;
char ac[3];
char bc[3];
char cc[3];
char dc[3];
int z;

for(int q = 0; q < 6; q++){
for(int i = 0; i < 7; i++){
if((i+1+ae) > 0 && (i+1+ae) < 10){
sprintf(ac,"0%d",i+1+ae);
}else if((i+1+ae) >= 10 && (i+1+ae) <= az){
sprintf(ac,"%d",i+1+ae);
}
z = 0;
parsedate(wbcount, i, ae, baumday, baummon, 1, z, outf, gb, ac, oy);
parsedate(papcount, i, ae, papday, papmon, 1, z, outf, bb, ac, oy);
parsedate(bsrcount, i, ae, bsrday, bsrmon, 1, z, outf, iv, ac, oy);
parsedate(albacount, i, ae, albaday, albamon, 1, z, outf, yb, ac, oy);
parsedate(orgcount, i, ae, orgday, orgmon, 1, z, outf, ob, ac, oy);
parsedate(biocount, i, ae, bioday, biomon, 1, z, outf, bv, ac, oy);
parsedate(biocount, i, ae, glasday, glasmon, 1, z, outf, db, ac, oy);
if((i+1+ae) > 0 && (i+1+ae) <= az && i == 6 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) == 1 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && i == 5 && z == 0){
outf << "\\hv{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && z == 0){
outf << ac << " & ";
}else if(z == 0){
outf << "& ";
}
}
ae += 7;
for(int i = 0; i < 7; i++){
if((i+1+be) > 0 && (i+1+be) < 10){
sprintf(bc, "0%d", i+1+be);
}else if((i+1+be) >= 10 && (i+1+be) <= bz){
sprintf(bc, "%d", i+1+be);
}
z = 0;
parsedate(wbcount, i, be, baumday, baummon, 2, z, outf, gb, bc, oy);
parsedate(papcount, i, be, papday, papmon, 2, z, outf, bb, bc, oy);
parsedate(bsrcount, i, be, bsrday, bsrmon, 2, z, outf, iv, bc, oy);
parsedate(albacount, i, be, albaday, albamon, 2, z, outf, yb, bc, oy);
parsedate(orgcount, i, be, orgday, orgmon, 2, z, outf, ob, bc, oy);
parsedate(biocount, i, be, bioday, biomon, 2, z, outf, bv, bc, oy);
parsedate(biocount, i, be, glasday, glasmon, 2, z, outf, db, bc, oy);
if((i+1+be) > 0 && (i+1+be) <= bz && i == 6 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && i == 5 && z == 0){
outf << "\\hv{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && z == 0){
outf << bc << " & ";
}else if(z == 0){
outf << "& ";
}
}
be += 7;
for(int i = 0; i < 7; i++){
if((i+1+ce) > 0 && (i+1+ce) < 10){
sprintf(cc, "0%d", i+1+ce);
}else if((i+1+ce) >= 10 && (i+1+ce) <= cz){
sprintf(cc, "%d", i+1+ce);
}
z = 0;
parsedate(wbcount, i, ce, baumday, baummon, 3, z, outf, gb, cc, oy);
parsedate(papcount, i, ce, papday, papmon, 3, z, outf, bb, cc, oy);
parsedate(bsrcount, i, ce, bsrday, bsrmon, 3, z, outf, iv, cc, oy);
parsedate(albacount, i, ce, albaday, albamon, 3, z, outf, yb, cc, oy);
parsedate(orgcount, i, ce, orgday, orgmon, 3, z, outf, ob, cc, oy);
parsedate(biocount, i, ce, bioday, biomon, 3, z, outf, bv, cc, oy);
parsedate(biocount, i, ce, glasday, glasmon, 3, z, outf, db, cc, oy);
if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 6 && z == 0){
outf << "\\rb{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 5 && z == 0){
outf << "\\hv{" << cc << "} & ";
}else if(mon_om == 3 && (i+1+ce) == tag_om && i == 0 && z == 0){
outf << "\\rb{" << cc << "} & ";
}else if(mon_kf == 3 && (i+1+ce) == tag_kf && i == 4 && z == 0){
outf << "\\rb{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && z == 0){
outf << cc << " & ";
}else if(z == 0){
outf << "& ";
}
}
ce += 7;
for(int i = 0; i < 7; i++){
if((i+1+de) > 0 && (i+1+de) < 10){
sprintf(dc, "0%d", i+1+de);
}else if((i+1+de) >= 10 && (i+1+de) <= dz){
sprintf(dc, "%d", i+1+de);
}
z = 0;
parsedate(wbcount, i, de, baumday, baummon, 4, z, outf, gb, dc, oy);
parsedate(papcount, i, de, papday, papmon, 4, z, outf, bb, dc, oy);
parsedate(bsrcount, i, de, bsrday, bsrmon, 4, z, outf, iv, dc, oy);
parsedate(albacount, i, de, albaday, albamon, 4, z, outf, yb, dc, oy);
parsedate(orgcount, i, de, orgday, orgmon, 4, z, outf, ob, dc, oy);
parsedate(biocount, i, de, bioday, biomon, 4, z, outf, bv, dc, oy);
parsedate(biocount, i, de, glasday, glasmon, 4, z, outf, db, dc, oy);
if((i+1+de) > 0 && (i+1+de) <= dz && i == 6 && z == 0){
outf << "\\rb{" << dc << "} ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && i == 5 && z == 0){
outf << "\\hv{" << dc << "} & ";
}else if(mon_om == 4 && (i+1+de) == tag_om && i == 0 && z == 0){
outf << "\\rb{" << dc << "} & ";
}else if(mon_kf == 4 && (i+1+de) == tag_kf && i == 4 && z == 0){
outf << "\\rb{" << dc << "} & ";
}else if(mon_hf == 4 && (i+1+de) == tag_hf && i == 3 && z == 0){
outf << "\\rb{" << dc << "} & ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && z == 0){
outf << dc << " & ";
}else if(i == 6 && z == 0){
outf << " ";
}else if(z == 0){
outf << " & ";
}
}
de += 7;
outf << "\\\\\n";
}
outf << "\\hline\n";
outf << "\\multicolumn{7}{|c|}{\\bf Mai} & \\multicolumn{7}{|c|}{\\bf Juni} & \\multicolumn{7}{|c|}{\\bf Juli} & \\multicolumn{7}{|c|}{\\bf August} \\\\\n";
a = wotag(1,5,year);
b = wotag(1,6,year);
c = wotag(1,7,year);
d = wotag(1,8,year);
az = 31;
bz = 30;
cz = 31;
dz = 31;
ae = -a;
be = -b;
ce = -c;
de = -d;

for(int q = 0; q < 6; q++){
for(int i = 0; i < 7; i++){
if((i+1+ae) > 0 && (i+1+ae) < 10){
sprintf(ac,"0%d",i+1+ae);
}else if((i+1+ae) >= 10 && (i+1+ae) <= az){
sprintf(ac,"%d",i+1+ae);
}
z = 0;
parsedate(wbcount, i, ae, baumday, baummon, 5, z, outf, gb, ac, oy);
parsedate(papcount, i, ae, papday, papmon, 5, z, outf, bb, ac, oy);
parsedate(bsrcount, i, ae, bsrday, bsrmon, 5, z, outf, iv, ac, oy);
parsedate(albacount, i, ae, albaday, albamon, 5, z, outf, yb, ac, oy);
parsedate(orgcount, i, ae, orgday, orgmon, 5, z, outf, ob, ac, oy);
parsedate(biocount, i, ae, bioday, biomon, 5, z, outf, bv, ac, oy);
parsedate(biocount, i, ae, glasday, glasmon, 5, z, outf, db, ac, oy);
if((i+1+ae) > 0 && (i+1+ae) <= az && i == 6 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) == 1 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && i == 5 && z == 0){
outf << "\\hv{" << ac << "} & ";
}else if(mon_hf == 5 && (i+1+ae) == tag_hf && i == 3 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if(mon_pm == 5 && (i+1+ae) == tag_pm && i == 0 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && z == 0){
outf << ac << " & ";
}else if(z == 0){
outf << "& ";
}
}
ae += 7;
for(int i = 0; i < 7; i++){
if((i+1+be) > 0 && (i+1+be) < 10){
sprintf(bc, "0%d", i+1+be);
}else if((i+1+be) >= 10 && (i+1+be) <= bz){
sprintf(bc, "%d", i+1+be);
}
z = 0;
parsedate(wbcount, i, be, baumday, baummon, 6, z, outf, gb, bc, oy);
parsedate(papcount, i, be, papday, papmon, 6, z, outf, bb, bc, oy);
parsedate(bsrcount, i, be, bsrday, bsrmon, 6, z, outf, iv, bc, oy);
parsedate(albacount, i, be, albaday, albamon, 6, z, outf, yb, bc, oy);
parsedate(orgcount, i, be, orgday, orgmon, 6, z, outf, ob, bc, oy);
parsedate(biocount, i, be, bioday, biomon, 6, z, outf, bv, bc, oy);
parsedate(biocount, i, be, glasday, glasmon, 6, z, outf, db, bc, oy);
if((i+1+be) > 0 && (i+1+be) <= bz && i == 6 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && i == 5 && z == 0){
outf << "\\hv{" << bc << "} & ";
}else if(mon_hf == 6 && (i+1+be) == tag_hf && i == 3 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if(mon_pm == 6 && (i+1+be) == tag_pm && i == 0 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && z == 0){
outf << bc << " & ";
}else if(z == 0){
outf << "& ";
}
}
be += 7;
for(int i = 0; i < 7; i++){
if((i+1+ce) > 0 && (i+1+ce) < 10){
sprintf(cc, "0%d", i+1+ce);
}else if((i+1+ce) >= 10 && (i+1+ce) <= cz){
sprintf(cc, "%d", i+1+ce);
}
z = 0;
parsedate(wbcount, i, ce, baumday, baummon, 7, z, outf, gb, cc, oy);
parsedate(papcount, i, ce, papday, papmon, 7, z, outf, bb, cc, oy);
parsedate(bsrcount, i, ce, bsrday, bsrmon, 7, z, outf, iv, cc, oy);
parsedate(albacount, i, ce, albaday, albamon, 7, z, outf, yb, cc, oy);
parsedate(orgcount, i, ce, orgday, orgmon, 7, z, outf, ob, cc, oy);
parsedate(biocount, i, ce, bioday, biomon, 7, z, outf, bv, cc, oy);
parsedate(biocount, i, ce, glasday, glasmon, 7, z, outf, db, cc, oy);
if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 6 && z == 0){
outf << "\\rb{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 5 && z == 0){
outf << "\\hv{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && z == 0){
outf << cc << " & ";
}else if(z == 0){
outf << "& ";
}
}
ce += 7;
for(int i = 0; i < 7; i++){
if((i+1+de) > 0 && (i+1+de) < 10){
sprintf(dc, "0%d", i+1+de);
}else if((i+1+de) >= 10 && (i+1+de) <= dz){
sprintf(dc, "%d", i+1+de);
}
z = 0;
parsedate(wbcount, i, de, baumday, baummon, 8, z, outf, gb, dc, oy);
parsedate(papcount, i, de, papday, papmon, 8, z, outf, bb, dc, oy);
parsedate(bsrcount, i, de, bsrday, bsrmon, 8, z, outf, iv, dc, oy);
parsedate(albacount, i, de, albaday, albamon, 8, z, outf, yb, dc, oy);
parsedate(orgcount, i, de, orgday, orgmon, 8, z, outf, ob, dc, oy);
parsedate(biocount, i, de, bioday, biomon, 8, z, outf, bv, dc, oy);
parsedate(biocount, i, de, glasday, glasmon, 8, z, outf, db, dc, oy);
if((i+1+de) > 0 && (i+1+de) <= dz && i == 6 && z == 0){
outf << "\\rb{" << dc << "} ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && i == 5 && z == 0){
outf << "\\hv{" << dc << "} & ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && z == 0){
outf << dc << " & ";
}else if(i == 6 && z == 0){
outf << " ";
}else if(z == 0){
outf << "& ";
}
}
de += 7;
outf << "\\\\\n";
}
outf << "\\hline\n";
outf << "\\multicolumn{7}{|c|}{\\bf September} & \\multicolumn{7}{|c|}{\\bf Oktober} & \\multicolumn{7}{|c|}{\\bf November} & \\multicolumn{7}{|c|}{\\bf Dezember} \\\\\n";
a = wotag(1,9,year);
b = wotag(1,10,year);
c = wotag(1,11,year);
d = wotag(1,12,year);
az = 30;
bz = 31;
cz = 30;
dz = 31;
ae = -a;
be = -b;
ce = -c;
de = -d;

for(int q = 0; q < 6; q++){
for(int i = 0; i < 7; i++){
if((i+1+ae) > 0 && (i+1+ae) < 10){
sprintf(ac,"0%d",i+1+ae);
}else if((i+1+ae) >= 10 && (i+1+ae) <= az){
sprintf(ac,"%d",i+1+ae);
}
z = 0;
parsedate(wbcount, i, ae, baumday, baummon, 9, z, outf, gb, ac, oy);
parsedate(papcount, i, ae, papday, papmon, 9, z, outf, bb, ac, oy);
parsedate(bsrcount, i, ae, bsrday, bsrmon, 9, z, outf, iv, ac, oy);
parsedate(albacount, i, ae, albaday, albamon, 9, z, outf, yb, ac, oy);
parsedate(orgcount, i, ae, orgday, orgmon, 9, z, outf, ob, ac, oy);
parsedate(biocount, i, ae, bioday, biomon, 9, z, outf, bv, ac, oy);
parsedate(biocount, i, ae, glasday, glasmon, 9, z, outf, db, ac, oy);
if((i+1+ae) > 0 && (i+1+ae) <= az && i == 6 && z == 0){
outf << "\\rb{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && i == 5 && z == 0){
outf << "\\hv{" << ac << "} & ";
}else if((i+1+ae) > 0 && (i+1+ae) <= az && z == 0){
outf << ac << " & ";
}else if(z == 0){
outf << "& ";
}
}
ae += 7;
for(int i = 0; i < 7; i++){
if((i+1+be) > 0 && (i+1+be) < 10){
sprintf(bc, "0%d", i+1+be);
}else if((i+1+be) >= 10 && (i+1+be) <= bz){
sprintf(bc, "%d", i+1+be);
}
z = 0;
parsedate(wbcount, i, be, baumday, baummon, 10, z, outf, gb, bc, oy);
parsedate(papcount, i, be, papday, papmon, 10, z, outf, bb, bc, oy);
parsedate(bsrcount, i, be, bsrday, bsrmon, 10, z, outf, iv, bc, oy);
parsedate(albacount, i, be, albaday, albamon, 10, z, outf, yb, bc, oy);
parsedate(orgcount, i, be, orgday, orgmon, 10, z, outf, ob, bc, oy);
parsedate(biocount, i, be, bioday, biomon, 10, z, outf, bv, bc, oy);
parsedate(biocount, i, be, glasday, glasmon, 10, z, outf, db, bc, oy);
if((i+1+be) > 0 && (i+1+be) <= bz && i == 6 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if((i+1+be) == 3 && z == 0){
outf << "\\rb{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && i == 5 && z == 0){
outf << "\\hv{" << bc << "} & ";
}else if((i+1+be) > 0 && (i+1+be) <= bz && z == 0){
outf << bc << " & ";
}else if(z == 0){
outf << "& ";
}
}
be += 7;
for(int i = 0; i < 7; i++){
if((i+1+ce) > 0 && (i+1+ce) < 10){
sprintf(cc, "0%d", i+1+ce);
}else if((i+1+ce) >= 10 && (i+1+ce) <= cz){
sprintf(cc, "%d", i+1+ce);
}
z = 0;
parsedate(wbcount, i, ce, baumday, baummon, 11, z, outf, gb, cc, oy);
parsedate(papcount, i, ce, papday, papmon, 11, z, outf, bb, cc, oy);
parsedate(bsrcount, i, ce, bsrday, bsrmon, 11, z, outf, iv, cc, oy);
parsedate(albacount, i, ce, albaday, albamon, 11, z, outf, yb, cc, oy);
parsedate(orgcount, i, ce, orgday, orgmon, 11, z, outf, ob, cc, oy);
parsedate(biocount, i, ce, bioday, biomon, 11, z, outf, bv, cc, oy);
parsedate(biocount, i, ce, glasday, glasmon, 11, z, outf, db, cc, oy);
if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 6 && z == 0){
outf << "\\rb{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && i == 5 && z == 0){
outf << "\\hv{" << cc << "} & ";
}else if((i+1+ce) > 0 && (i+1+ce) <= cz && z == 0){
outf << cc << " & ";
}else if(z == 0){
outf << "& ";
}
}
ce += 7;
for(int i = 0; i < 7; i++){
if((i+1+de) > 0 && (i+1+de) < 10){
sprintf(dc, "0%d", i+1+de);
}else if((i+1+de) >= 10 && (i+1+de) <= dz){
sprintf(dc, "%d", i+1+de);
}
z = 0;
parsedate(wbcount, i, de, baumday, baummon, 12, z, outf, gb, dc, oy);
parsedate(papcount, i, de, papday, papmon, 12, z, outf, bb, dc, oy);
parsedate(bsrcount, i, de, bsrday, bsrmon, 12, z, outf, iv, dc, oy);
parsedate(albacount, i, de, albaday, albamon, 12, z, outf, yb, dc, oy);
parsedate(orgcount, i, de, orgday, orgmon, 12, z, outf, ob, dc, oy);
parsedate(biocount, i, de, bioday, biomon, 12, z, outf, bv, dc, oy);
parsedate(biocount, i, de, glasday, glasmon, 12, z, outf, db, dc, oy);
if((i+1+de) > 0 && (i+1+de) <= dz && i == 6 && z == 0){
outf << "\\rb{" << dc << "} ";
}else if((i+1+de) == 25 || (i+1+de) == 26 ){
outf << "\\rb{" << dc << "} & ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && i == 5 && z == 0){
outf << "\\hv{" << dc << "} & ";
}else if(((i+1+de) == 24 || (i+1+de) == 31) && z == 0){
outf << "\\hv{" << dc << "} & ";
}else if((i+1+de) > 0 && (i+1+de) <= dz && z == 0){
outf << dc << " & ";
}else if(i == 6 && z == 0){
outf << " ";
}else if(z == 0){
outf << "& ";
}
}
de += 7;
outf << "\\\\\n";
}
outf << "\\hline\n";
outf << "\\end{tabular}\n";
outf << "\\end{document}\n";
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

int wotag(int day, int month, int year){
 int tag;  int result;  int interim;
 int century = gauss((double) year/100.);
 int jahr = year - 100*century;
 int month_new;
 if(month == 1){
  month_new = 11;
  if(jahr == 0){
   jahr = 99;
   century = century - 1;
  }else{
   jahr = jahr - 1;
  }
 }else if(month == 2){
  month_new = 12;
   if(jahr == 0){
    jahr = 99;
    century = century - 1;
   }else{
    jahr = jahr - 1;
   }
  }else{
   month_new = month - 2;
  }
 interim = day + gauss(2.6*(double) month_new - 0.2) + jahr + gauss((double) jahr/4.) + gauss((double) century/4.) - 2*century;
 result = interim%7;
 while(result < 0){
  result = result + 7;
 }
switch(result){
 case 0: 
	tag = 6;
	break;
 case 1: 
	tag = 0;
	break;
 case 2: 
	tag = 1;
	break;
 case 3: 
	tag = 2;
	break;
 case 4: 
	tag = 3;
	break;
 case 5: 
	tag = 4;
	break;
 case 6: 
	tag = 5;
	break;
}
return(tag);
}

void ostermontag(int year, int *tag_om, int *mon_om, int *tag_os, int *mon_os){
int g    = year%19;
int c    = gauss((double) year/100.);
int a    = gauss((double) c/4.);
int d    = gauss((8.*(double) c + 13.)/25.);
int h    = (int)(c - a - d + 19*g + 15)%30;
int f    = gauss((double) h/20.);
int k    = gauss(29./((double) h + 1.));
int m    = gauss((21.-(double) g)/11.);
int n    = gauss((double) year/4.);
int i    = h - f*(1 - k*m);
int j    = (int)(year + n + i + 2 - c + a)%7;
int l    = i - j;

if(l<=3){
 *mon_os = 3;
 *tag_os = l + 28;
}else{
 *mon_os = 4;
 *tag_os = l - 3;
}

if(*tag_os == 31 && *mon_os == 3){
 *tag_om = 1;
 *mon_om = 4;
}else{
 *tag_om = *tag_os + 1;
 *mon_om = *mon_os;
}
}

void karfreitag(int tag_os, int mon_os, int *tag_kf, int *mon_kf){
if(tag_os < 3 && mon_os == 4){
 *tag_kf = tag_os + 29;
 *mon_kf = 3;
}else{
 *tag_kf = tag_os - 2;
 *mon_kf = mon_os;
}
}

void himmelfahrt(int tag_os, int mon_os, int *tag_hf, int *mon_hf){
if(tag_os < 23 && mon_os == 3){
 *tag_hf = tag_os + 8;
 *mon_hf = 4;
}else if(tag_os > 22 && mon_os == 3){
 *tag_hf = tag_os - 22;
 *mon_hf = 5;
}else if(tag_os > 21 && mon_os == 4){
 *tag_hf = tag_os - 22;
 *mon_hf = 6;
}else{
 *tag_hf = tag_os + 9;
 *mon_hf = 5;
}
}

void pfingsten(int tag_hf, int mon_hf, int *tag_pm, int *mon_pm){
int tag_ps = 0;
int mon_ps = 0;
if(mon_hf == 4){
 tag_ps = tag_hf - 20;
 mon_ps = 5;
}else if(mon_hf == 5 && tag_hf < 22){
 tag_ps = tag_hf + 10;
 mon_ps = 5;
}else if(mon_hf == 5 && tag_hf > 21){
 tag_ps = tag_hf - 21;
 mon_ps = 6;
}else if(mon_hf == 6){
 tag_ps = tag_hf + 10;
 mon_ps = 6;
}

if(mon_ps == 5 && tag_ps == 31){
 *tag_pm = 1;
 *mon_pm = 6;
}else{
 *tag_pm = tag_ps + 1;
 *mon_pm = mon_ps;
}
}

void parsedate(int count, int i, int e, int* day, int* month, int aktmonth, int &z, ofstream& outf, char* inn, char* c, char* off){
 for(int j = 0; j < count; j++){
  if((i+1+e) == day[j] && month[j] == aktmonth && z == 0){
   outf << inn << c << off;
   z++;
   break;
  }
 }
}

void read_data(int count, int* day, int* month, ifstream& inf){
 for(int i = 0; i < count; i++){
  inf >> day[i] >> month[i];
 }
}

int leapyear(int year){
 int q;
 if(year%400 == 0){
  q = 1;
 }else if(year%100 == 0){
  q = 0;
 }else if(year%4 == 0){
  q = 1;
 }else{
  q = 0;
 }
 return(q);
}

