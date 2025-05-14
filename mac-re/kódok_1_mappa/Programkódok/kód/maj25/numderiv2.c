#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define ngrid 		1000				/*	a gridpontok szama			*/
#define rmin 		0.1				/*	0.1 AU-tol szamol			*/
#define rmax 		100.0				/*	15 AU-ig szamol				*/
#define dgrid 		(rmax-rmin)/(ngrid-1)		/*	dr 					*/

#define sigma0	 	0.0001				/*	init. sigma ertek			*/	
#define sdexp 		-0.5 				/*	surface density profile exponent 	*/
#define alpha_visc 	0.01				/*	alpha parameter				*/
#define asp_ratio 	0.05				/*	aspect ratio				*/

#define M_star 		1.0				/*	classical T Tauri mass			*/

//#define G		1.0 				//ekkor r=1-nél a periódus 2*pi

#define G 		0.01720209895			/*	Gauss-grav const			*/
#define AU2CM		1.496e13			/*	1 AU --> 1 cm				*/
#define GCM2MSAU	1.125211e-7			/*	g/cm2 --> M_Sun/AU2			*/
#define MS2AUDAY	5.7754e-7			/*	m/s --> AU/day				*/
#define dt		0.1				/*	timestep				*/
#define timemax		1000.0				/*	meddig fusson a program			*/
	
FILE *fmo;


/*	angular keplerian frequency in AU/?		*/
double kep_freq(double r){					
  
  	double omega;
  	omega = sqrt((G * G * M_star) / (r * r * r));
  	return omega;
  
}


/*	viscosity in the disk	*/
double visc(double r){
 
  	double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r, kepler;

	kepler = kep_freq(r);

  	r_dze_i = 2.0;					/*	a belso deadzone belso hatara CSE-ben			*/
 	r_dze_o = 10.0;					/*	a kulso deadzone belso hatara CSE-ben			*/
  	Dr_dze_i = 0.1;					/*	a belso deadzonehatar szelessege CSE-ben		*/
  	Dr_dze_o = 0.5;					/*	a kulso deadzonehatar szelessege CSE-ben		*/
  	a_mod = 0.01;

  	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));	/* 	10 AU-nal a_mod-dal megvaltoztatja a viszkozitas merteket		*/
  	nu = alpha_visc * alpha_r * asp_ratio * asp_ratio * r * r * kepler;
  
  	return nu;
  
}

/*	local speed of sound in AU/?	*/
double speed_of_sound(double r) {

	double kepler;
	kepler = kep_freq(r);
	return asp_ratio * r * kepler;

}


/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
double Coeff_1(double r){					
  
  	double A;
  	A = 3.0 * visc(r);
  	return A;
  
}


/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/		
double Coeff_2(double r){					
  
  	double B;
  	B = 9.0 * visc(r) / (2.0 * r);
  	return B;
    
}


/*	boundary condition for sigma	*/
void Perem(double *sigmavec) {					
  
  	sigmavec[0] = (sigmavec[1] / visc(rmin)) * visc(rmin - dgrid);
  	sigmavec[ngrid+1]=(sigmavec[ngrid] / visc(rmin + (ngrid) * dgrid)) * visc(rmin + (ngrid+1) * dgrid);
  
}


/*	initial surface density profile in M_Sol / AU / AU	*/ 
void Initial_Profile(double *sigmavec, double *r){		

  	int i;
  
  	for(i = 1; i <= ngrid; i++) {
    		sigmavec[i] = sigma0 * pow(r[i],sdexp);		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}

  	Perem(sigmavec);
  
}


/*	density in the mid-plane in M_sol / AU / AU / AU	*/
double rho_midpl(double r, double sigma) {

	return sigma/(sqrt(2.0 * M_PI) * asp_ratio * r);

}


/*	pressure if we assume isothermal disk --> M_SOL / AU / ? / ?	*/
double pressure(double r, double sigma) {

	double cs, rho;
	cs = speed_of_sound(r);
	rho = rho_midpl(r,sigma);
	return cs * cs * rho;

}


void pressure_der(double *p, double *dp) {
/* i 1 - n-1 és p(0) = p(1) és p(n+1) = p(n) */

	int i;
	
	for(i = 1; i <= ngrid; i++) {	
		dp[i] = (p[i+1] - p[i-1]) / (2.0 * dgrid);
//		printf("%lg \n", dp[i]);
	}

	dp[0] = dp[1];
	dp[ngrid+1] = dp[ngrid];

}


int main(){

  	int i;
  	double time, sigmavec[ngrid+2], rvec[ngrid+2], r, pressvec[ngrid+2], dpressvec[ngrid+2];
	char dens_name[64];

  
  	time = 0.0;

	for(i = 0; i <= ngrid+1; i++) {							/*	load an array of radii	*/
 		rvec[i] = rmin + (i-1) * dgrid;
	}

  	Initial_Profile(sigmavec,rvec);							/*	load the inital surface density profile	*/

	for(i = 0; i <= ngrid+1; i++) {							/*	load an array of pressure	*/
 		pressvec[i] = pressure(rvec[i],sigmavec[i]);
	}

	pressure_der(pressvec,dpressvec);

  	do {



		if(fmod(time, (int)(timemax/10.0)) < dt || time == 0) {			/*	a futasi ido alatt 10-szer irja ki, hogy hol tart	*/
			printf("Az eltelt ido: %lg nap \n",time);	
			snprintf(dens_name,64,"dens_%i.dat", (int)(time));		/*	puts the timestep in the name of the file that contains selected information	*/
  			
			fmo = fopen(dens_name,"w");
 			for(i = 1; i <= ngrid; i++) {
	   			fprintf(fmo, "%lg  %lg  %lg  %lg\n", rvec[i], sigmavec[i],pressvec[i],dpressvec[i]); 	/*	print sigma in a file in each selected timestep	*/
   
 			}
			
			fclose(fmo);  

		}

  
  		for(i = 1; i <= ngrid; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	*/
    			sigmavec[i] = (sigmavec[i] * visc(rvec[i]) + dt * (Coeff_1(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - 2.0 * sigmavec[i] * visc(rvec[i]) + sigmavec[i-1] * visc(rvec[i-1])) / (dgrid * dgrid) + Coeff_2(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - sigmavec[i-1] * visc(rvec[i-1])) / (2.0 * dgrid))) / visc(rvec[i]);
	 		pressvec[i] = pressure(rvec[i],sigmavec[i]);
    
  		}
	
		pressvec[0] = pressvec[1];
		pressvec[ngrid+1] = pressvec[ngrid];
  		Perem(sigmavec);						/*	loading boundary condition in each timestep	*/
		pressure_der(pressvec,dpressvec);
  		time = time + dt;						/*	timestepping	*/



  	} while(time < timemax);

	return 0;

}

