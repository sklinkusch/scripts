#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

int gauss(double zahl)
{
  int result;
  if (zahl > 0)
  {
    result = (int)zahl;
  }
  else if (zahl < 0)
  {
    result = (int)zahl - 1;
  }
  else
  {
    result = 0;
  }
  return (result);
}

char *wotag(int d, int m, int j)
{
  int result;
  int interim;
  char *tag = new char[3];
  int c = gauss((double)j / 100.);
  int y = j - 100 * c;
  int m_neu;
  if (m == 1)
  {
    m_neu = 11;
    if (y == 0)
    {
      y = 99;
      c = c - 1;
    }
    else
    {
      y = y - 1;
    }
  }
  else if (m == 2)
  {
    m_neu = 12;
    if (y == 0)
    {
      y = 99;
      c = c - 1;
    }
    else
    {
      y = y - 1;
    }
  }
  else
  {
    m_neu = m - 2;
  }
  interim = d + gauss(2.6 * (double)m_neu - 0.2) + y + gauss((double)y / 4.) + gauss((double)c / 4.) - 2 * c;
  result = interim % 7;
  while (result < 0)
  {
    result = result + 7;
  }
  if (result == 0)
  {
    sprintf(tag, "So");
  }
  else if (result == 1)
  {
    sprintf(tag, "Mo");
  }
  else if (result == 2)
  {
    sprintf(tag, "Di");
  }
  else if (result == 3)
  {
    sprintf(tag, "Mi");
  }
  else if (result == 4)
  {
    sprintf(tag, "Do");
  }
  else if (result == 5)
  {
    sprintf(tag, "Fr");
  }
  else if (result == 6)
  {
    sprintf(tag, "Sa");
  }
  return (tag);
}

int wochentag(int d, int m, int j)
{
  int result;
  int interim;
  int c = gauss((double)j / 100.);
  int y = j - 100 * c;
  int m_neu;
  if (m == 1)
  {
    m_neu = 11;
    if (y == 0)
    {
      y = 99;
      c = c - 1;
    }
    else
    {
      y = y - 1;
    }
  }
  else
  {
    m_neu = m - 2;
  }
  interim = d + gauss(2.6 * (double)m_neu - 0.2) + y + gauss((double)y / 4.) + gauss((double)c / 4.) - 2 * c;
  result = interim % 7;
  while (result < 0)
  {
    result = result + 7;
  }
  return (result);
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cerr << "Bitte Jahr eingeben!"
         << "\n";
    exit(1);
  }

  int schalt = 0;
  char *monat = new char[10];
  int monnr;
  int tag;
  int tag_rosen;
  int tag_fast;
  int tag_asch;
  int tag_kafr;
  int tag_osmo;
  int tag_hifa;
  int tag_mutt;
  int tag_pfiso;
  int tag_pfimo;
  int tag_trini;
  int tag_fron;
  int tag_volk;
  int tag_buss;
  int tag_tot;
  int tag_one;
  int tag_two;
  int tag_three;
  int tag_four;
  int monat_rosen;
  int monat_fast;
  int monat_asch;
  int monat_kafr;
  int monat_osmo;
  int monat_hifa;
  int monat_mutt;
  int monat_pfiso;
  int monat_pfimo;
  int monat_trini;
  int monat_fron;
  int monat_one;

  int jahr = atoi(argv[1]);
  if (jahr % 4 == 0 && jahr % 100 != 0)
  {
    schalt = 1;
  }
  else if (jahr % 400 == 0)
  {
    schalt = 1;
  }
  else
  {
    schalt = 0;
  }

  int g = jahr % 19;
  // cout << "g = " << g << "\n";
  int c = gauss((double)jahr / 100.);
  // cout << "c = " << c << "\n";
  int a = gauss((double)c / 4.);
  // cout << "a = " << a << "\n";
  int d = gauss((8. * (double)c + 13.) / 25.);
  // cout << "d = " << d << "\n";
  int h = (int)(c - a - d + 19 * g + 15) % 30;
  // cout << "h = " << h << "\n";
  int f = gauss((double)h / 20.);
  // cout << "f = " << f << "\n";
  int k = gauss(29. / ((double)h + 1.));
  // cout << "k = " << k << "\n";
  int m = gauss((21. - (double)g) / 11.);
  // cout << "m = " << m << "\n";
  int n = gauss((double)jahr / 4.);
  // cout << "n = " << n << "\n";
  int i = h - f * (1 - k * m);
  // cout << "i = " << i << "\n";
  int j = (int)(jahr + n + i + 2 - c + a) % 7;
  // cout << "j = " << j << "\n";
  int l = i - j;
  // cout << "l = " << l << "\n";

  if (l <= 3)
  {
    sprintf(monat, "Maerz");
    monnr = 3;
    tag = l + 28;
  }
  else
  {
    sprintf(monat, "April");
    monnr = 4;
    tag = l - 3;
  }

  if (l > 20)
  {
    tag_rosen = tag - 17;
    monat_rosen = 3;
  }
  else if (l < 21 && l > 3 && schalt == 1)
  {
    tag_rosen = 29 + tag - 17;
    monat_rosen = 2;
  }
  else if (l < 21 && l > 3 && schalt == 0)
  {
    tag_rosen = 28 + tag - 17;
    monat_rosen = 2;
  }
  else if (l < 4 && schalt == 1)
  {
    tag_rosen = tag - 19;
    monat_rosen = 2;
  }
  else if (l < 4 && schalt == 0)
  {
    tag_rosen = tag - 20;
    monat_rosen = 2;
  }
  else
  {
    tag_rosen = 0;
    monat_rosen = 0;
  }

  if (tag_rosen == 28 && schalt == 0)
  {
    tag_fast = 1;
    monat_fast = 3;
  }
  else if (tag_rosen == 29 && schalt == 1)
  {
    tag_fast = 1;
    monat_fast = 3;
  }
  else
  {
    tag_fast = tag_rosen + 1;
    monat_fast = monat_rosen;
  }

  if (tag_fast == 28 && schalt == 0)
  {
    tag_asch = 1;
    monat_asch = 3;
  }
  else if (tag_fast == 29 && schalt == 1)
  {
    tag_asch = 1;
    monat_asch = 3;
  }
  else
  {
    tag_asch = tag_fast + 1;
    monat_asch = monat_fast;
  }

  if (tag < 3 && monnr == 4)
  {
    tag_kafr = tag + 29;
    monat_kafr = 3;
  }
  else
  {
    tag_kafr = tag - 2;
    monat_kafr = monnr;
  }

  if (tag == 31 && monnr == 3)
  {
    tag_osmo = 1;
    monat_osmo = 4;
  }
  else
  {
    tag_osmo = tag + 1;
    monat_osmo = monnr;
  }

  if (tag < 23 && monnr == 3)
  {
    tag_hifa = tag + 8;
    monat_hifa = 4;
  }
  else if (tag > 22 && monnr == 3)
  {
    tag_hifa = tag - 22;
    monat_hifa = 5;
  }
  else if (tag > 21 && monnr == 4)
  {
    tag_hifa = tag - 22;
    monat_hifa = 6;
  }
  else
  {
    tag_hifa = tag + 9;
    monat_hifa = 5;
  }

  int mut_dex = wochentag(1, 5, jahr);
  tag_mutt = 1 - mut_dex + 14;
  monat_mutt = 5;

  if (monat_hifa == 4)
  {
    tag_pfiso = tag_hifa - 20;
    monat_pfiso = 5;
  }
  else if (monat_hifa == 5 && tag_hifa < 22)
  {
    tag_pfiso = tag_hifa + 10;
    monat_pfiso = 5;
  }
  else if (monat_hifa == 5 && tag_hifa > 21)
  {
    tag_pfiso = tag_hifa - 21;
    monat_pfiso = 6;
  }
  else if (monat_hifa == 6)
  {
    tag_pfiso = tag_hifa + 10;
    monat_pfiso = 6;
  }

  if (monat_pfiso == 5 && tag_pfiso == 31)
  {
    tag_pfimo = 1;
    monat_pfimo = 6;
  }
  else
  {
    tag_pfimo = tag_pfiso + 1;
    monat_pfimo = monat_pfiso;
  }

  if (monat_pfimo == 5 && tag_pfimo > 25)
  {
    tag_trini = tag_pfimo - 25;
    monat_trini = 6;
  }
  else
  {
    tag_trini = tag_pfimo + 6;
    monat_trini = monat_pfimo;
  }

  if (monat_trini == 5 && tag_trini > 27)
  {
    tag_fron = tag_trini - 27;
    monat_fron = 6;
  }
  else
  {
    tag_fron = tag_trini + 4;
    monat_fron = monat_trini;
  }
  tag_volk = 19 - wochentag(19, 11, jahr);
  tag_buss = tag_volk + 3;
  tag_tot = 26 - wochentag(26, 11, jahr);
  int day = 3 - wochentag(3, 12, jahr);
  if (day < 1)
  {
    tag_one = 30 + day;
    monat_one = 11;
  }
  else
  {
    tag_one = day;
    monat_one = 12;
  }
  tag_two = 10 - wochentag(10, 12, jahr);
  tag_three = 17 - wochentag(17, 12, jahr);
  tag_four = 24 - wochentag(24, 12, jahr);

  char *wotag_neujahr = new char[3];
  char *wotag_epiph = new char[3];
  char *wotag_frau = new char[3];
  char *wotag_mai = new char[3];
  char *wotag_maria = new char[3];
  char *wotag_einheit = new char[3];
  char *wotag_reform = new char[3];
  char *wotag_allheil = new char[3];
  char *wotag_allseel = new char[3];
  char *wotag_heilig = new char[3];
  char *wotag_1weihn = new char[3];
  char *wotag_2weihn = new char[3];
  char *wotag_silv = new char[3];

  wotag_neujahr = wotag(1, 1, jahr);
  wotag_epiph = wotag(6, 1, jahr);
  wotag_frau = wotag(8, 3, jahr);
  wotag_mai = wotag(1, 5, jahr);
  wotag_maria = wotag(15, 8, jahr);
  wotag_einheit = wotag(3, 10, jahr);
  wotag_reform = wotag(31, 10, jahr);
  wotag_allheil = wotag(1, 11, jahr);
  wotag_allseel = wotag(2, 11, jahr);
  wotag_heilig = wotag(24, 12, jahr);
  wotag_1weihn = wotag(25, 12, jahr);
  wotag_2weihn = wotag(26, 12, jahr);
  wotag_silv = wotag(31, 12, jahr);

  cout << "Neujahr am 1.1." << jahr << " (" << wotag_neujahr << ")\n";
  cout << "Heilige Drei Koenige (Epiphanias) am 6.1." << jahr << " (" << wotag_epiph << ")\n";
  if (monat_rosen == 2 || tag_rosen < 6)
  {
    cout << "Rosenmontag am " << tag_rosen << "." << monat_rosen << "." << jahr << " (Mo)\n";
    cout << "Fastnacht am " << tag_fast << "." << monat_fast << "." << jahr << " (Di)\n";
    cout << "Aschermittwoch am " << tag_asch << "." << monat_asch << "." << jahr << " (Mi)\n";
    cout << "Internationaler Frauentag am 8.3." << jahr << " (" << wotag_frau << ")\n";
  }
  else if (tag_rosen == 6)
  {
    cout << "Rosenmontag am 06.03." << jahr << " (Mo)\n";
    cout << "Fastnacht am 07.03." << jahr << " (Di)\n";
    cout << "Internationaler Frauentag und Aschermittwoch am 8.3." << jahr << " (Mi)\n";
  }
  else if (tag_rosen == 7)
  {
    cout << "Rosenmontag am 7.3." << jahr << " (Mo)\n";
    cout << "Internationaler Frauentag und Fastnacht am 8.3." << jahr << " (Di)\n";
    cout << "Aschermittwoch am 9.3." << jahr << " (Mi)\n";
  }
  else
  {
    cout << "Internationaler Frauentag und Rosenmontag am 8.3." << jahr << " (Mo)\n";
    cout << "Fastnacht am 9.3." << jahr << " (Di)\n";
    cout << "Aschermittwoch am 10.3." << jahr << " (Mi)\n";
  }
  cout << "Karfreitag am " << tag_kafr << "." << monat_kafr << "." << jahr << " (Fr)\n";
  cout << "Ostersonntag am " << tag << "." << monnr << "." << jahr << " (So)\n";
  cout << "Ostermontag am " << tag_osmo << "." << monat_osmo << "." << jahr << " (Mo)\n";
  if (tag_hifa == 30 && monat_hifa == 4)
  {
    cout << "Christi Himmelfahrt am 30.4." << jahr << " (Do)\n";
  }
  if (tag_hifa == 1 && monat_hifa == 5)
  {
    cout << "Christi Himmelfahrt und Maifeiertag am 1.5." << jahr << " (Do)\n";
  }
  else
  {
    cout << "Maifeiertag am 1.5." << jahr << " (" << wotag_mai << ")\n";
  }
  if (monat_hifa == 5 && tag_hifa > 1 && tag_hifa < tag_mutt)
  {
    cout << "Christi Himmelfahrt am " << tag_hifa << ".5." << jahr << " (Do)\n";
  }
  cout << "Muttertag am " << tag_mutt << "." << monat_mutt << "." << jahr << " (So)\n";
  if ((monat_hifa == 5 && tag_hifa > tag_mutt) || (monat_hifa == 6))
  {
    cout << "Christi Himmelfahrt am " << tag_hifa << "." << monat_hifa << "." << jahr << " (Do)\n";
  }
  cout << "Pfingstsonntag am " << tag_pfiso << "." << monat_pfiso << "." << jahr << " (So)\n";
  cout << "Pfingstmontag am " << tag_pfimo << "." << monat_pfimo << "." << jahr << " (Mo)\n";
  cout << "Trinitatis am " << tag_trini << "." << monat_trini << "." << jahr << " (So)\n";
  cout << "Fronleichnam am " << tag_fron << "." << monat_fron << "." << jahr << " (Do)\n";
  cout << "Mariae Himmelfahrt am 15.8." << jahr << " (" << wotag_maria << ")\n";
  cout << "Tag der Deutschen Einheit am 3.10." << jahr << " (" << wotag_einheit << ")\n";
  cout << "Reformationstag am 31.10." << jahr << " (" << wotag_reform << ")\n";
  cout << "Allerheiligen am 1.11." << jahr << " (" << wotag_allheil << ")\n";
  cout << "Allerseelen am 2.11." << jahr << " (" << wotag_allseel << ")\n";
  cout << "Volkstrauertag am " << tag_volk << ".11." << jahr << " (So)\n";
  cout << "Buß- und Bettag am " << tag_buss << ".11." << jahr << " (Mi)\n";
  cout << "Totensonntag am " << tag_tot << ".11." << jahr << " (So)\n";
  cout << "1. Advent am " << tag_one << "." << monat_one << "." << jahr << " (So)\n";
  cout << "2. Advent am " << tag_two << ".12." << jahr << " (So)\n";
  cout << "3. Advent am " << tag_three << ".12." << jahr << " (So)\n";
  cout << "4. Advent am " << tag_four << ".12." << jahr << " (So)\n";
  cout << "Heiligabend am 24.12." << jahr << " (" << wotag_heilig << ")\n";
  cout << "1. Weihnachtstag am 25.12." << jahr << " (" << wotag_1weihn << ")\n";
  cout << "2. Weihnachtstag am 26.12." << jahr << " (" << wotag_2weihn << ")\n";
  cout << "Silvester am 31.12." << jahr << " (" << wotag_silv << ")\n";
}
