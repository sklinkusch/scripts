# include<iostream>
# include<fstream> 
# include<stdio.h>
# include<stdlib.h>
# include<string.h>

/*****************************************************************************************************************
* GMS_INP                                                                                                        *
* Routine zur Bestimmung von Ordnungszahlen                                                                      *
*                                                                                    Stefan Klinkusch, 2008-2009 *
*****************************************************************************************************************/

using namespace std;

double calc_oz(char symbol1, char symbol2){

double oz;
int sym[2];
sym[0] = (int) symbol1;
sym[1] = (int) symbol2;
if(sym[0] == 'A' || sym[0] == 'a') {
    if(sym[1] == 'c' || sym[1] == 'C') oz = 89.0;        //Ac, actinium, 89.0
    else if(sym[1] == 'g' || sym[1] == 'G') oz = 47.0;   //Ag, silver, 47.0
    else if(sym[1] == 'l' || sym[1] == 'L') oz = 13.0;   //Al, aluminium, 13.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 95.0;   //Am, americium, 95.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 18.0;   //Ar, argon, 18.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 33.0;   //As, arsenic, 33.0
    else if(sym[1] == 't' || sym[1] == 'T') oz = 85.0;   //At, astatine, 85.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 79.0;   //Au, gold, 79.0
    else exit(1);
}else if(sym[0] == 'B' || sym[0] == 'b'){
    if(sym[1] == '\0') oz = 5.0;                         //B, boron, 5.0
    else if(sym[1] == 'a' || sym[1] == 'A') oz = 56.0;   //Ba, barium, 56.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 4.0;    //Be, beryllium, 4.0
    else if(sym[1] == 'h' || sym[1] == 'H') oz = 107.0;  //Bh, bohrium, 107.0
    else if(sym[1] == 'i' || sym[1] == 'I') oz = 83.0;   //Bi, bismuth, 83.0
    else if(sym[1] == 'k' || sym[1] == 'K') oz = 97.0;   //Bk, berkelium, 97.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 35.0;   //Br, bromine, 35.0
    else exit(1);
}else if(sym[0] == 'C' || sym[0] == 'c'){
    if(sym[1] == '\0') oz = 6.0;                         //C, carbon, 6.0
    else if(sym[1] == 'a' || sym[1] == 'A') oz = 20.0;   //Ca, calcium, 20.0
    else if(sym[1] == 'd' || sym[1] == 'D') oz = 48.0;   //Cd, cadmium, 48.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 58.0;   //Ce, cerium, 58.0
    else if(sym[1] == 'f' || sym[1] == 'F') oz = 98.0;   //Cf, califonium, 98.0
    else if(sym[1] == 'l' || sym[1] == 'L') oz = 17.0;   //Cl, chlorine, 17.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 96.0;   //Cm, curium, 96.0
    else if(sym[1] == 'o' || sym[1] == 'O') oz = 27.0;   //Co, cobalt, 27.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 24.0;   //Cr, chromium, 24.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 55.0;   //Cs, cesium, 55.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 29.0;   //Cu, copper, 29.0
    else exit(1);
}else if(sym[0] == 'D' || sym[0] == 'd'){
    if(sym[1] == 'b' || sym[1] == 'B') oz = 105.0;       //Db, dubnium, 105.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 110.0;  //Ds, darmstadtium, 110.0
    else if(sym[1] == 'y' || sym[1] == 'Y') oz = 66.0;   //Dy, dysprosium, 66.0
    else exit(1);
}else if(sym[0] == 'E' || sym[0] == 'e'){
    if(sym[1] == 'r' || sym[1] == 'R') oz = 68.0;        //Er, erbium, 68.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 99.0;   //Es, einsteinium, 99.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 63.0;   //Eu, europium, 63.0
    else exit(1);
}else if(sym[0] == 'F' || sym[0] == 'f'){
    if(sym[1] == '\0') oz = 9.0;                         //F, fluorine, 9.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 26.0;   //Fe, iron, 26.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 100.0;  //Fm, fermium, 100.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 87.0;   //Fr, francium, 87.0
    else exit(1);
}else if(sym[0] == 'G' || sym[0] == 'g'){
    if(sym[1] == 'a' || sym[1] == 'A') oz = 31.0;        //Ga, gallium, 31.0
    else if(sym[1] == 'd' || sym[1] == 'D') oz = 64.0;   //Gd, gadolinium, 64.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 32.0;   //Ge, germanium, 32.0
    else exit(1);
}else if(sym[0] == 'H' || sym[0] == 'h'){
    if(sym[1] == '\0') oz = 1;                           //H, hydrogen, 1.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 2.0;    //He, helium, 2.0
    else if(sym[1] == 'f' || sym[1] == 'F') oz = 72.0;   //Hf, hafnium, 72.0
    else if(sym[1] == 'g' || sym[1] == 'G') oz = 80.0;   //Hg, mercury, 80.0
    else if(sym[1] == 'o' || sym[1] == 'O') oz = 67.0;   //Ho, holmium, 67.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 108.0;  //Hs, hassium, 108.0
    else exit(1);
}else if(sym[0] == 'I' || sym[0] == 'i'){
    if(sym[1] == '\0') oz = 53.0;                        //I, iodine, 53.0
    else if(sym[1] == 'n' || sym[1] == 'N') oz = 49.0;   //In, indium, 49.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 77.0;   //Ir, iridium, 77.0
    else exit(1);
}else if(sym[0] == 'K' || sym[0] == 'k'){
    if(sym[1] == '\0') oz = 19.0;                        //K, potassium, 19.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 36.0;   //Kr, krypton, 36.0
    else exit(1);
}else if(sym[0] == 'L' || sym[0] == 'l'){
    if(sym[1] == 'a' || sym[1] == 'A') oz = 57.0;        //La, lanthanum, 57.0
    else if(sym[1] == 'i' || sym[1] == 'I') oz = 3.0;    //Li, lithium, 3.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 103.0;  //Lr, lawrencium, 103.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 71.0;   //Lu, lutetium, 71.0
    else exit(1);
}else if(sym[0] == 'M' || sym[0] == 'm'){
    if(sym[1] == 'd' || sym[1] == 'D') oz = 101.0;       //Md, mendelejevium, 101.0
    else if(sym[1] == 'g' || sym[1] == 'G') oz = 12.0;   //Mg, magnesium, 12.0
    else if(sym[1] == 'n' || sym[1] == 'N') oz = 25.0;   //Mn, manganese, 25.0
    else if(sym[1] == 'o' || sym[1] == 'O') oz = 42.0;   //Mo, molybdenum, 42.0
    else if(sym[1] == 't' || sym[1] == 'T') oz = 109.0;  //Mt, meitnerium, 109.0
    else exit(1);
}else if(sym[0] == 'N' || sym[0] == 'n'){
    if(sym[1] == '\0') oz = 7.0;                         //N, nitrogen, 7.0
    else if(sym[1] == 'a' || sym[1] == 'A') oz = 11.0;   //Na, sodium, 11.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 41.0;   //Nb, niobium, 41.0
    else if(sym[1] == 'd' || sym[1] == 'D') oz = 60.0;   //Nd, neodymium, 60.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 10.0;   //Ne, neon, 10.0
    else if(sym[1] == 'i' || sym[1] == 'I') oz = 28.0;   //Ni, nickel, 28.0
    else if(sym[1] == 'o' || sym[1] == 'O') oz = 102.0;  //No, nobelium, 102.0
    else if(sym[1] == 'p' || sym[1] == 'P') oz = 93.0;   //Np, neptunium, 93.0
    else exit(1);
}else if(sym[0] == 'O' || sym[0] == 'o'){
    if(sym[1] == '\0') oz = 8.0;                         //O, oxygen, 8.0
    else if(sym[1] == 's' || sym[1] == 'S') oz = 76.0;   //Os, osmium, 76.0
    else exit(1);
}else if(sym[0] == 'P' || sym[0] == 'p'){
    if(sym[1] == '\0') oz = 15.0;                        //P, phosphorus, 15.0
    else if(sym[1] == 'a' || sym[1] == 'A') oz = 91.0;   //Pa, protactinium, 91.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 82.0;   //Pb, lead, 82.0
    else if(sym[1] == 'd' || sym[1] == 'D') oz = 46.0;   //Pd, palladium, 46.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 61.0;   //Pm, prometium, 61.0
    else if(sym[1] == 'o' || sym[1] == 'O') oz = 84.0;   //Po, polonium, 84.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 59.0;   //Pr, praseodymium, 59.0
    else if(sym[1] == 't' || sym[1] == 'T') oz = 78.0;   //Pt, platinum, 78.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 94.0;   //Pu, plutonium, 94.0
    else exit(1);
}else if(sym[0] == 'R' || sym[0] == 'r'){
    if(sym[1] == 'a' || sym[1] == 'A') oz = 88.0;        //Ra, radium, 88.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 37.0;   //Rb, rubidium, 37.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 75.0;   //Re, rhenium, 75.0
    else if(sym[1] == 'f' || sym[1] == 'F') oz = 104.0;  //Rf, rutherfordium, 104.0
    else if(sym[1] == 'g' || sym[1] == 'G') oz = 111.0;  //Rg, roentgenium, 111.0
    else if(sym[1] == 'h' || sym[1] == 'H') oz = 45.0;   //Rh, rhodium, 45.0
    else if(sym[1] == 'n' || sym[1] == 'N') oz = 86.0;   //Rn, radon, 86.0
    else if(sym[1] == 'u' || sym[1] == 'U') oz = 44.0;   //Ru, ruthenium, 44.0
    else exit(1);
}else if(sym[0] == 'S' || sym[0] == 's'){
    if(sym[1] == '\0') oz = 16.0;                        //S, sulphur, 16.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 51.0;   //Sb, antimony, 51.0
    else if(sym[1] == 'c' || sym[1] == 'C') oz = 21.0;   //Sc, scandium, 21.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 34.0;   //Se, selenium, 34.0
    else if(sym[1] == 'g' || sym[1] == 'G') oz = 106.0;  //Sg, seaborgium, 106.0
    else if(sym[1] == 'i' || sym[1] == 'I') oz = 14.0;   //Si, silicon, 14.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 62.0;   //Sm, samarium, 62.0
    else if(sym[1] == 'n' || sym[1] == 'N') oz = 50.0;   //Sn, tin, 50.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 38.0;   //Sr, strontium, 38.0
    else exit(1);
}else if(sym[0] == 'T' || sym[0] == 't'){
    if(sym[1] == 'a' || sym[1] == 'A') oz = 73.0;        //Ta, tantalum, 73.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 65.0;   //Tb, terbium, 65.0
    else if(sym[1] == 'c' || sym[1] == 'C') oz = 43.0;   //Tc, technetium, 43.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 52.0;   //Te, tellurium, 52.0
    else if(sym[1] == 'h' || sym[1] == 'H') oz = 90.0;   //Th, thorium, 90.0
    else if(sym[1] == 'i' || sym[1] == 'I') oz = 22.0;   //Ti, titanium, 22.0
    else if(sym[1] == 'l' || sym[1] == 'L') oz = 81.0;   //Tl, thallium, 81.0
    else if(sym[1] == 'm' || sym[1] == 'M') oz = 69.0;   //Tm, thulium, 69.0
    else exit(1);
}else if(sym[0] == 'U' || sym[0] == 'u'){
    if(sym[1] == '\0') oz = 92.0;                        //U, uranium, 92.0
    else exit(1);
}else if(sym[0] == 'V' || sym[0] == 'v'){
    if(sym[1] == '\0') oz = 23.0;                        //V, vanadium, 23.0
    else exit(1);
}else if(sym[0] == 'W' || sym[0] == 'w'){
    if(sym[1] == '\0') oz = 74.0;                        //W, tungsten, 74.0
    else exit(1);
}else if(sym[0] == 'X' || sym[0] == 'x'){
    if(sym[1] == '\0') oz = 0.0;                         //X, dummy atom, 0.0
    else if(sym[1] == 'e' || sym[1] == 'E') oz = 54.0;   //Xe, xenon, 54.0
    else exit(1);
}else if(sym[0] == 'Y' || sym[0] == 'y'){
    if(sym[1] == '\0') oz = 39.0;                        //Y, yttrium, 39.0
    else if(sym[1] == 'b' || sym[1] == 'B') oz = 70.0;   //Yb, ytterbium, 70.0
    else exit(1);
}else if(sym[0] == 'Z' || sym[1] == 'z'){
    if(sym[1] == 'n' || sym[1] == 'N') oz = 30.0;        //Zn, zinc, 30.0
    else if(sym[1] == 'r' || sym[1] == 'R') oz = 40.0;   //Zr, zirconium, 40.0
    else exit(1);
}else{
exit(1);
} 
return(oz);
}

