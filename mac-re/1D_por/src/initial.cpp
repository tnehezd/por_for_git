#include <iostream>
#include "funcs.h"
#include "define.h"
#include <math.h>

double LoadR(int i) {

	return RMIN + (i-1)*DD;

}


long double InitialSigma(int i, double r) {

	return SIGMA_0 * pow(r, SDEXP);
	
}


double Initialize(double *r, double *sigma, double *u) {

	for (int i = 1; i < NGRID+2; i++) {

		r[i] = LoadR(i);
		sigma[i] = InitialSigma(i,r[i]);
		u[i] = sigma[i] * Viscosity(r[i]);

	}	

	u[0] = u[(int)NGRID+1] = 0.0;
	sigma[0] = sigma[(int)NGRID+1] = 0.0;
	r[0] = r[1] - DD;
	r[(int)NGRID+1] = r[(int)NGRID] + DD;
	
}
