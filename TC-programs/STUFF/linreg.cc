#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Usage: ./linreg <nr of points> <input file>\n";
	exit(0);
    }
    int nrop = atoi(argv[1]);             // number of points
    double* xvals = new double[nrop];     // x values
    double* yvals = new double[nrop];     // y values
    double  xysum = 0.;
    double  xsum  = 0.;
    double  ysum  = 0.;
    double  xsqsum = 0.;
    char infile[256];                     // input file
    sprintf(infile, "%s", argv[2]);
    ifstream inf;
    inf.open(infile);
    for(int i = 0; i < nrop; i++){
	inf >> xvals[i] >> yvals[i];
	xsum += xvals[i];
	ysum += yvals[i];
	xysum += (xvals[i]*yvals[i]);
	xsqsum += (xvals[i]*xvals[i]);
    }
    double xsumsq = xsum*xsum;
    double nenner = (nrop * xsqsum) - xsumsq;
    double xav = xsum/nrop;
    double yav = ysum/nrop;
    double slope = ((nrop * xysum) - (xsum * ysum))/nenner;
    double invslope = 1./slope;
    double ordin = ((xsqsum * ysum) - (xsum * xysum))/nenner;
    double xmxav = 0.;
    double ymyav = 0.;
    double xymav = 0.;
    double xmavs = 0.;
    double ymavs = 0.;
    double ydev = 0.;
    for(int i = 0; i < nrop; i++){
	xmxav = xvals[i] - xav;
	ymyav = yvals[i] - yav;
	xymav += (xmxav * ymyav);
	xmavs += (xmxav * xmxav);
	ymavs += (ymyav * ymyav);
	ydev += pow(yvals[i] - (ordin + (slope * xvals[i])),2.);
    }
    double r   = xymav / sqrt(xmavs*ymavs);
    double rsq = r*r;
    double devslope = sqrt((nrop * ydev)/((nrop - 2)*nenner));
    double reldevslope = 100.*devslope/slope;
    double devordin = sqrt((ydev * xsqsum)/((nrop - 2)*nenner));
    double reldevordin = 100.*devordin/ordin;
    double devinvslope = devslope/pow(slope,2.);
    double reldevinvslope = 100.*devinvslope/invslope;
    cout << "Result of linear regression analysis: \n";
    cout << "Slope of the line: " << slope << ", deviation: " << devslope << " (" << reldevslope << "%)\n";
    cout << "Inverted slope: " << invslope << ", deviation: " << devinvslope << " (" << reldevinvslope << "%)\n";
    cout << "Y-intercept of the line: " << ordin << ", deviation: " << devordin << " (" << reldevordin << "%)\n";
    cout << "R^2: " << rsq << "\n";
}

