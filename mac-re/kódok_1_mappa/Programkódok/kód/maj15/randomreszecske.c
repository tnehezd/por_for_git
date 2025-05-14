#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	min	10.0
#define	max	500.0
#define pwind	-7			/*	distribution power			*/
#define N	2000			/*	number of random generated numbers	*/
#define chunk 	10.0
#define dr	((max-min)/chunk)

FILE *distrfil;

/*	random general szamokat locmin és locmax kozott --> ez egyenletes eloszlast ad	*/
double random_gen(double locmin,double locmax) {

	return ((fabs(locmin) + fabs(locmax)) * rand()/((double)RAND_MAX + fabs(locmax)) - fabs(locmin));

}


/*		egyenletes eloszlasbol (pl. 0 es 1 kozott) hatvanyfuggveny eloszlas:
					    x = [(x1^(n+1) - x0^(n+1))*y + x0^(n+1)]^(1/(n+1))   
				http://mathworld.wolfram.com/RandomNumber.html		*/
double uniform_to_power_law() {

	double y, temp1, temp2, temp3, temp4, temp5;	/*	az egyenletes eloszlasvol hatvanyfuggveny eloszlashoz szukseges lokalis valtozok: temp1,temp2,temp3,temp4,temp5		*/
 	double locmin, locmax;
	locmin = 0.0;
	locmax = 1.0;

	y = random_gen(locmin,locmax);				/*	random szamok generator fuggveny meghivasa (es az y valtozo feltoltese a generalt szammal)	*/
	temp1 = pwind+1.0;			
	temp2 = pow(max,temp1);
	temp3 = pow(min,temp1);
 	temp4 = 1.0/temp1;
	temp5 = (temp2 - temp3)*y + temp3;
	
	return pow(temp5,temp4);
  
}

/*	sorting the elements of the random generated pow-law distributed array		*/
void sort(double *mass) {
  
	int i, j;
	double temp;
    
	for(i = 1; i < N; i++) {
        	for(j=0;j < N-i;j++) {
	    		if(mass[j] > mass[j+1]) { 
	        		temp = mass[j];
				mass[j] = mass[j+1];
				mass[j+1] = temp;
	   		}
		}
	} 
  
}




int main() {

	double radius[N],xsort[N], count,size,temp,counttemp,countloc,xrange[(int)dr],x[N],y[N],z[N],vx[N],vy[N],vz[N],mass[N];
	int i,h;
	srand(time(NULL));				/*	a randomszam generalast cpu idohoz koti, igy valoban random szamokat kapunk		*/
	
	size = 0.;
	count = 0.;
	temp = 0.;
	count = 0.;


 	double locmin1, locmax1, locmin2, locmax2, locmin3, locmax3;	/*	min max ertekek a random hely- es sebessegvektor generalashoz	*/
	locmin1 = -10.0;		/*	csillagaszati egysegben, x es y koordinatakhoz	*/
	locmax1 = 10.0;
	locmin2 = -2.0;
	locmax2 = 2.0;		/*	csillagaszati egysegben a z koordinatakhoz	*/
	locmin3 = -1.0e-2;
	locmax3 = 1.0e-2;

	

	
	distrfil = fopen("distr.dat","w");


	for(i = 0; i < N; i++) {	
		radius[i] = uniform_to_power_law();
		x[i] = random_gen(locmin1,locmax1);
		y[i] = random_gen(locmin1,locmax1);
		z[i] = random_gen(locmin2,locmax2);
		vx[i] = random_gen(locmin3,locmax3);
		vy[i] = random_gen(locmin3,locmax3);
		vz[i] = random_gen(locmin3,locmax3);
	}

	sort(mass);

	for(i = 0; i < N; i++) {	
		
		fprintf(distrfil,"%e  %e  %e  %e  %e  %e  %e\n",radius[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]);
	}

	fclose(distrfil);	


	return 0;

}
