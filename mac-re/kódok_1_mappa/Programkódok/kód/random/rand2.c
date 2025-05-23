#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

FILE *distrfil;

double RAND_NORMAL(double mu, double sigma) {


	return (mu + (rand()%2 ? -1.0 : 1.0)*sigma*pow(-log(0.99999*((double) rand()/RAND_MAX)), 0.5));

}


int main () {

	srand(time(NULL));
	double r[10], sigma, mu;
	int i;

	mu = 0.0;
	sigma = 1.0;
	
	distrfil = fopen("distr2.dat","w");	
	
	for (i = 0; i < 10; i++) {
		r[i] = RAND_NORMAL(mu,sigma);
		fprintf(distrfil,"%lf\n",r[i]);
		printf("%lg\n",r[i]);		
	}

	fclose(distrfil);

	return 0;

}

