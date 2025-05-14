#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	min	-1.0
#define	max	1.0

FILE *distrfil;

/*	random general szamokat min és max kozott	*/
double random_gen() {

	return ((fabs(min) + fabs(max)) * rand()/((double)RAND_MAX + fabs(max)) - fabs(min));

}

/*	a random szamokat betoltjuk a normal-eloszlasba	*/
double distr(double sigma, double m, double x) {
	
	double eloszl;
	eloszl = 1.0 / (sigma * sqrt(2.0 * M_PI)) * exp(-((x-m) * (x-m)) / (2.0 * sigma * sigma));
	return eloszl;

}

int main() {
	srand(time(NULL));
	double distrib[10], sigma, m, x[10];
	int i;
	
	sigma = 1.0;
	m = 0;


	distrfil = fopen("distr.dat","w");	
	
	for (i = 0; i < 10; i++) {
		x[i] = random_gen();			/*	random szamok min es max kozott		*/
		printf("%lg\n",x[i]);
		distrib[i] = distr(sigma,m,x[i]);	/* 	a normal-eloszlas letrehozasa		*/
		fprintf(distrfil,"%lf\t%lf\n",x[i],distrib[i]);
	}

	fclose(distrfil);

	return 0;
}

