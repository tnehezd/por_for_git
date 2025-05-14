#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	min	100.0
#define	max	500.0
#define pwind	-7			/*	distribution power			*/
#define N	200			/*	number of random generated numbers	*/
#define chunk 	10.0
#define dr	((max-min)/chunk)

FILE *distrfil, *histfil;

/*	random general szamokat locmin és locmax kozott --> ez egyenletes eloszlast ad	*/
double random_gen() {

	double locmin, locmax;
	locmin = 0.0;
	locmax = 1.0;
	return ((fabs(locmin) + fabs(locmax)) * rand()/((double)RAND_MAX + fabs(locmax)) - fabs(locmin));

}


/*		egyenletes eloszlasbol (pl. 0 es 1 kozott) hatvanyfuggveny eloszlas:
					    x = [(x1^(n+1) - x0^(n+1))*y + x0^(n+1)]^(1/(n+1))   
				http://mathworld.wolfram.com/RandomNumber.html		*/
double uniform_to_power_law() {

	double y, temp1, temp2, temp3, temp4, temp5;	/*	az egyenletes eloszlasvol hatvanyfuggveny eloszlashoz szukseges lokalis valtozok: temp1,temp2,temp3,temp4,temp5		*/
 
	y = random_gen();				/*	random szamok generator fuggveny meghivasa (es az y valtozo feltoltese a generalt szammal)	*/
	temp1 = pwind+1.0;			
	temp2 = pow(max,temp1);
	temp3 = pow(min,temp1);
 	temp4 = 1.0/temp1;
	temp5 = (temp2 - temp3)*y + temp3;
	
	return pow(temp5,temp4);
  
}

/*	sorting the elements of the random generated pow-law distributed array		*/
void sort(double *x) {
  
	int i, j;
	double temp;
    
	for(i = 1; i < N; i++) {
        	for(j=0;j < N-i;j++) {
	    		if(x[j] > x[j+1]) { 
	        		temp = x[j];
				x[j] = x[j+1];
				x[j+1] = temp;
	   		}
		}
	} 
  
}


int main() {

	double x[N],xsort[N], count,size,temp,counttemp,countloc,xrange[(int)dr];
	int i,h;
	srand(time(NULL));				/*	a randomszam generalast cpu idohoz koti, igy valoban random szamokat kapunk		*/
	
	size = 0.;
	count = 0.;
	temp = 0.;
	count = 0.;
	
	distrfil = fopen("distr.dat","w");
	histfil = fopen("hist.dat","w");

	for(i = 0; i < N; i++) {	
		x[i] = uniform_to_power_law();
	}

	sort(x);

	for(i = 0; i < N; i++) {	
		
		fprintf(distrfil,"%lg\n",x[i]);
	}

	fclose(distrfil);	

	countloc = min + chunk;
	counttemp = 0.0;
	

	for(h = 1; h < dr; h++) {
	  
		counttemp = 0.;
		for(i = 0; i < N; i++) {
			if ((countloc - chunk) < x[i] && x[i] < countloc) {
				counttemp++;

			}

		}
		
		xrange[h] = counttemp;
		fprintf(histfil,"%lg\n",xrange[h]);
		countloc = countloc + chunk;
	}
	fclose(histfil);

 	printf("%lg  %lg  %lg\n",(double)sizeof(x)/sizeof(double), size, chunk);
	
	return 0;

}
