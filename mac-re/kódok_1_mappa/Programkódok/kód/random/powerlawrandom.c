#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	min	0.0
#define	max	1.0

FILE *distrfil;

/*	random general szamokat min és max kozott	*/
double random_gen() {

	return (rand()/((double)RAND_MAX));

}

/*	a random szamokat betoltjuk a normal-eloszlasba	*/
double distr(double sigma, double mu, double prob) {
	
	double x;
	x = pow(prob,(1.0/(-7.0)));
	return x;

}

int main() {
	srand(time(NULL));
	double distrib[100], sigma, m, prob[100],x[100], dgrid;
	int i;
	
	sigma = 1.0;
	m = 0;

	dgrid = 2.0 / 100.0;

	distrfil = fopen("distr.dat","w");	
	
	for (i = 0; i < 100; i++) {
		prob[i] = random_gen();			/*	random szamok min es max kozott		*/
		distrib[i] = distr(sigma,m,prob[i]);	/* 	a normal-eloszlas letrehozasa		*/
		x[i] = -1.0 + (i-1) * dgrid;
		printf("%lg  %lg\n",prob[i],distrib[i]);
		fprintf(distrfil,"%lf\t%lf\t%lg\n",prob[i],distrib[i], x[i]);
	}

	fclose(distrfil);

	return 0;
}

