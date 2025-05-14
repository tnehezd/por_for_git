#ifndef NUMERICS_H
#define NUMERICS_H

double TimeStep(double *r);
double Coeff_1(double r);
double Coeff_2(double r);
void NumDeriv(double *r, double *u, double *w, double dt);

#endif
