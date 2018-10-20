#include <iostream> // Bibliothek zum Einlesen/Ausschreiben aus/in Konsole
#include <fstream>  // Bibliothek zum Einlesen aus Dateien/Ausschreiben in Dateien
#include <math.h>   // Bibliothek fuer mathematische Funktionen
#include <stdio.h>  // Standard-Input-Output-Bibliothek
#include <stdlib.h> // Standard-Bibliothek
#include <string.h> // Bibliothek fuer Zeichenketten

using namespace std; 

/******************************************************************************************
* Programm zum Berechnen des Restguthabens bei der PB-Sofortrente                         *
*                                                                 Stefan Klinkusch, 2013  *
******************************************************************************************/

int main(int argc, char* argv[]){
// Abfrage ob genug Parameter eingegeben wurden, sonst Abbruch und Fehlermeldung
 if(argc != 4){
  cerr << "Gebrauch: ./sofortrente <Zinssatz> <Auszahlbetrag> <Ausgabedatei>\n";
  exit(1);
 }
// Definition von Variablen
 double zinssatz;                                // jaehrl. Zinssatz
 double G_o;                                     // Startguthaben
 double rente;                                   // monatlicher Auszahlbetrag
 double G_start;                                 // Guthaben zu Beginn eines Monats
 int i;                                          // Zahl der Monate (Laufvariable)
 double G;                                       // Guthaben nach dem Monat
 double zins;                                    // errechnete Zinsen
// Zuweisung von Werten
 zinssatz = strtod(argv[1],NULL);                // Zinssatz (in %) einlesen
 zinssatz /= 100.;                               // Zinssatz ohne % 
 G_o  = 100000.;                                 // Startguthaben: 100.000 EUR
 rente = strtod(argv[2],NULL);                   // Rente einlesen
 ofstream outf;                                  // Ausgabedatei
 outf.open(argv[3]);                             // Ausgabedatei oeffnen
// Berechnung fuer Monat 0 (bei Abschluss der Rente)
 G_start = G_o;                                  // Kapital zu Beginn des Monats
 i = 0;                                          // Zahl der Monate (Laufvariable)
 outf << "# Monate      # Guthaben (EUR)\n";     // schreibe Titelzeile in Datei
 outf << i << "       " << G_start << "\n";      // schreibe Monat 0 in Datei
// Schleife ueber 360 Monate (30 Jahre), Berechnung und Ausschreiben in Datei
 for(i = 1; i <= 360;  i++){                    
  zins = G_start * zinssatz / 12.;               // Zinsen fuer den Monat
  G = G_start + zins - rente;                    // neues Guthaben
  outf << i << "       " << G << "\n";           // neues Guthaben in Datei
  G_start = G;                                   // altes Guthaben fuer Monat i+1
  outf.flush();                                  // Zwischenspeichern
 }
 outf.close();                                   // Schliessen der Ausgabedatei
}

