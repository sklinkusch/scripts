#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void carttosph(double* cart, double* sph);
double signum(double zahl);
void sphtocart(double* sph, double* cart);
void pi_calc(double* dip, double* fld, double sigma);
    
int main(void){
    double* dip_cart = new double[3];
    double* dip_sph  = new double[3];
    double* fld_cart = new double[3];
    double* fld_sph  = new double[3];
    double sigma;
    cout << "Enter x component of dipole moment (a.u.):";
    cin >> dip_cart[0];
    cout << "Enter y component of dipole moment (a.u.):";
    cin >> dip_cart[1];
    cout << "Enter z component of dipole moment (a.u.):";
    cin >> dip_cart[2];
    cout << "Enter FWHM of pulse (a.u.):";
    cin >> sigma;
    carttosph(dip_cart,dip_sph);
    pi_calc(dip_sph,fld_sph,sigma);
    sphtocart(fld_sph,fld_cart);
    cout << "x component of electric field (a.u.): " << fld_cart[0] << "\n";
    cout << "y component of electric field (a.u.): " << fld_cart[1] << "\n";
    cout << "z component of electric field (a.u.): " << fld_cart[2] << "\n";
}

void carttosph(double* cart, double* sph){
    sph[0] = sqrt(pow(cart[0],2) + pow(cart[1],2) + pow(cart[2],2));
    if(sph[0] == 0.){
	sph[1] = 0.;
    }else{
	sph[1] = acos(cart[2]/sph[0]);
    }
    if(cart[0] > 0.){
	sph[2] = atan(cart[1]/cart[0]);
    }else if(cart[0] == 0.){	
	sph[2] = signum(cart[1])*M_PI/2.;
    }else{
	if(cart[1] < 0.){
	    sph[2] = atan(cart[1]/cart[0]) - M_PI;
	}else{
	    sph[2] = atan(cart[1]/cart[0]) + M_PI;
	}
    }
}

double signum(double zahl){
    int x;
    if(zahl > 0){
	x = 1.;
    }else if(zahl < 0){
	x = -1.;
    }else{
	x = 0.;
    }
    return(x);
}

void sphtocart(double* sph, double* cart){
    cart[0] = sph[0] * sin(sph[1]) * cos(sph[2]);
    cart[1] = sph[0] * sin(sph[1]) * sin(sph[2]);
    cart[2] = sph[0] * cos(sph[1]);
}

void pi_calc(double* dip, double* fld, double sigma){
    fld[0] = M_PI/(dip[0]*sigma);
    fld[1] = dip[1];
    fld[2] = dip[2];
}

