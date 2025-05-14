#include "numerics.h"
#include "define.h"
#include "funcs.h"
#include <algorithm>
#include <iostream>

double TimeStep(double *r) {					/*	a diffúziós egyenlet megoldásához szükséges minimális időlépés	*/

	double A_max, stepping;
	int i;

	A_max = -10000.0;
	
	for(i = 0; i < NGRID; i++) {
		if(Coeff_1(r[i]) > A_max) {
			A_max = Coeff_1(r[i]);
		}
	}

	stepping = DD * DD / (2.0 * A_max);

	return stepping;

}


double Coeff_1(double r) {

	return 3.0*Viscosity(r);

}


double Coeff_2(double r) {

	return 9.0*Viscosity(r)/(2*r);

}



void NumDeriv(double *r, double *u, double *w, double dt) {

		u[0] = u[(int)NGRID+1] = 0.0;
		w[0] = w[(int)NGRID+1] = 0.0;


		for (int i = 1; i < NGRID+2; i++) {

			double du = Coeff_1(r[i])*(u[i+1]-2*u[i]+u[i-1])/(DD*DD) + Coeff_2(r[i])*(u[i+1]-u[i])/DD;
			w[i] = u[i] + dt * du;
		}

//		std::swap(w,u);	

		for (int i = 0; i < NGRID+2; i++) {

			u[i] = w[i];
		}


}
