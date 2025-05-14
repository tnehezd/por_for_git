#include <iostream>
#include "funcs.h"
#include "define.h"
#include "math.h"

float KepFreq(float r) {
	return sqrt((G*MSTAR)/pow(r,3.0));
}


float SoundOfSpeed(float r) {
	return KepFreq(r) * ASPRATIO;

}

double Viscosity(float r){

/* 
  	double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r;

  	r_dze_i=1.0;
  r_dze_o=10.0;
  Dr_dze_i=0.5;
  Dr_dze_o=0.5;
  a_mod=0.01;
  
*/

	double nu;
	double alpha_r;
	double a_mod = 0.01;
	double r_dze_i = 1.0;
	double r_dze_o = 20.;
	double Dr_dze_i = 2.0 * r_dze_i * ASPRATIO;
	double Dr_dze_o = 2.0 * r_dze_o * ASPRATIO;

	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
  
	nu = ALPHA*ASPRATIO*ASPRATIO*G*sqrt(r)*alpha_r;
  
  	return nu;
  
}

void Kepler() {

	unsigned int i;
	float r;

	for(i = 1; i <= NGRID+1; i++) {
    
		r = RMIN + (i-1) * DD;
//		kep_freq(r);
		std::cout << KepFreq(r)<< " r: " << r << " cs: " << KepFreq(r) << " nu: " << Viscosity(r) << std::endl;
    	}

}


void CalculateSigma(double *r, double *u, double *sigma) {

	unsigned int i;
	for(i = 1; i <= NGRID+1; i++) {
	
		sigma[i] = u[i] / Viscosity(r[i]);
	    
    	}
}
