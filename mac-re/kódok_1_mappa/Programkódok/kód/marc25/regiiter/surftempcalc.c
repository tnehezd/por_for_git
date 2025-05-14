#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define ngrid 		1000				/*	a gridpontok szama				*/
#define rmin 		0.1				/*	0.1 AU-tol szamol				*/
#define rmax 		15.0				/*	15 AU-ig szamol					*/
#define dgrid 		(rmax-rmin)/(ngrid-1)		/*	dr 						*/

#define stconst 	5.67e-5				/* 	erg/cm2/s/K4, Stefan-B const			*/
#define sigma0	 	0.0001				/*	init. sigma ertek				*/	
#define sdexp 		-0.5 				/*	surface density profile exponent 		*/
#define alpha_visc 	0.01				/*	alpha parameter					*/
#define asp_ratio 	0.05				/*	aspect ratio					*/

#define M_star 		1				/*	classical T Tauri mass				*/
#define R_star		0.0116275			/*	2.5 AU (T Tauri) radius				*/
#define T_star		4000				/*	K (classical T Tauri temp)			*/

#define C 		6.938685e18			/*	homerseklet kiszamolasahoz a ket oldal dimenziojanak osszeegyeztetese (dimenziok: nap, AU, M_Nap, Kelvin) * 9/4	*/

#define G 		0.01720209895			/*	Gauss-grav const				*/
#define k_boltzman	1.38e-23
#define adiabat		1.4				/*	diatomic ideal gass				*/
#define mH		1.67262e-27 * 2.0	 	/*	mass of hydrogen molecule in kg			*/
#define AU2CM		1.496e13			/*	1 AU --> 1 cm					*/
#define GCM2MSAU	1.125211e-7			/*	g/cm2 --> M_Sun/AU2				*/
#define MS2AUDAY	5.7754e-7			/*	m/s --> AU/day					*/
#define dt		0.5				/*	timestep					*/
#define timemax		100.0			/*	meddig fusson a program				*/
#define WO		10.				/*	hány lépésenként írja ki az adatokat file-ba	*/

FILE *surfinitdat, *surftimedat, *tempinitdat, *temptimedat, *f_dens, *filgnu, *filplot, *f_temp, *filgnu2, *filplot2;
 


double sigmainit(double r) {

	return sigma0*pow(r,sdexp);					/*	initial surface density profile		*/

}


double speed_of_sound(double temp){					/*	local speed of sound			*/
	
	double c;
	c = sqrt(temp * k_boltzman * adiabat / mH);
  	c = c * 5.77547828e-7; 						/*	m/s to AU/day				*/
	return c;

}


double kep_freq(double r){						/*	angular keplerian frequency		*/
  	
  	double omega;
  	omega = sqrt((G * M_star) / (r * r * r));
  	return omega;
  
}


double scale_height(double temp, double r){				/*	local scale height			*/

	return speed_of_sound(temp) / (kep_freq(r) * sqrt(adiabat));

}


double visc(double temp, double r){
 
  	double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r, c_s;

  	r_dze_i = 2.0;					/*	a belso deadzone belso hatara CSE-ben		*/
 	r_dze_o = 10.0;					/*	a kulso deadzone belso hatara CSE-ben		*/
  	Dr_dze_i = 0.1;					/*	a belso deadzonehatar szelessege CSE-ben		*/
  	Dr_dze_o = 0.5;					/*	a kulso deadzonehatar szelessege CSE-ben		*/
  	a_mod = 0.01;

  	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));	/* 	10 AU-nal a_mod-dal megvaltoztatja a viszkozitas merteket		*/
//	alpha_r = 1.0;
	c_s = speed_of_sound(temp);
	nu = alpha_visc * alpha_r * c_s * c_s / kep_freq(r);


//  	nu = alpha_visc * alpha_r * asp_ratio * asp_ratio * G * sqrt(r);						/*	AU2/D		*/
  
  	return nu;
  
}


double rho0(double temp, double r, double localsigma) {		/* midplane density	*/

	double rho_mp;
	rho_mp = localsigma / GCM2MSAU / (sqrt(2.0 * M_PI) * scale_height(temp,r) * AU2CM);
	return rho_mp;

}


void Perem(double *sigmavec, double *tempvec) {					/*	boundary condition for sigma		*/
  
  	sigmavec[0] = (sigmavec[1] / visc(tempvec[1],rmin)) * visc(tempvec[0],rmin - dgrid);
  	sigmavec[ngrid+1]=(sigmavec[ngrid] / visc(tempvec[ngrid],rmin + (ngrid) * dgrid)) * visc(tempvec[ngrid+1],rmin + (ngrid+1) * dgrid);
  
}


void Initial_Profile(double *sigmavec, double *r, double *tempvec){		/*	initial profile of sigma		*/

  	int i;
  
	surfinitdat = fopen("surfinit.dat","w");

  	for(i = 1; i <= ngrid; i++) {
    		sigmavec[i] = sigma0 * pow(r[i],sdexp);		/*	sigma0*r^x (x could be eg. -1/2)	*/
		fprintf(surfinitdat,"%lg   %lg  \n", r[i], sigmavec[i]);
  	}

  	Perem(sigmavec,tempvec);
	fclose(surfinitdat);  

}


double Coeff_1(double temp, double r){					/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
  
  	double A;
  	A = 3.0 * visc(temp,r);
  	return A;
  
}


double Coeff_2(double temp, double r){					/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/		
  
  	double B;
  	B = 9.0 * visc(temp,r) / (2.0 * r);
  	return B;
    
}


double Coeff_m(){						/*	for solving 2 * (stefan-boltzman_const) * T^4 = (3 * tau / 8 + sqrt(3) / 4 + 1 / (4 * tau) + 2 * (stefan-boltzman_const) * T_background^4 --> Coeff_m = 2 * (stefan-boltzman_const))	*/
 
  	double m_sigma;
  	m_sigma = 2.0 * stconst;
  	return m_sigma;
  
}


double kappa(double temp, double r, double localsigma) {

	double T1,T2,T3,T4,T5,T6;
  	double kk,a,b,kap,l,h,k;
    
	l=4.6e3 * pow(rho0(temp,r,localsigma), 1.0 / 15.0);
	h=1.1e4 * pow(rho0(temp,r,localsigma), 1.0 / 21.0);
	k=3.0e4 * pow(rho0(temp,r,localsigma), 4.0 / 75.0);
     
  	kk = 0.0;
  	T1 = 170.0;
  	T2 = 210.0;
  	T3 = l;		/*	1347.1765	*/
  	T4 = 3000.0;
  	T5 = h;		/*	4575.518379	*/
  	T6 = k;		/*	11231.93517	*/
	      
	if(temp <= T1){ 
		kap = 2.0e-4;
		a = 0.0;
		b = 2.0;	
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}
			      
	else if((temp > T1) && (temp < T2)){
		kap = 2.0e16;
		a = 0.0; 
		b = -7.0;	 
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}

      
	else if((temp > T2) && (temp < T3)){
		kap = 5.0e-3;
		a = 0.0; 
		b = 1.0;
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}
		      
	else if((temp > T3) && (temp < T4)){
		kap = 2.0e34;
		a = 2.0/3.0; 
		b = -9.0;	
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}

			      
	else if((temp > T4) && (temp < T5)){
		kap = 2.0e-8;
		a = 2.0/3.0; 
		b = 3.0;	
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}
			      
	else if((temp > T5) && (temp < T6)){
		kap = 1.0e-36;
		a = 1.0/3.0; 
		b = 10.0;	
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}
			      
	else if(temp > T6){
		kap = 1.5e20;
		a = 1.0; 
		b = -5.0/2.0;	
		kk = kap * pow(temp,b) * pow((rho0(temp,r,localsigma)),a);
	}
	
	return kk;
	
}


double opt_dept(double temp, double r, double localsigma) {

	double tau;
	tau = kappa(temp,r,localsigma) * localsigma / GCM2MSAU / 2.0;
	return tau;

}


double t_eff(double temp, double r, double localsigma) {

	double taueff;
	taueff = (3.0 * opt_dept(temp,r,localsigma) / 8.0) + (sqrt(3.0) / 4.0) + (1.0 / (4.0 * opt_dept(temp,r,localsigma)));
	return taueff;	

}


double star_heat(double r){
	
	return ((1.0 / 10.0) * ((R_star * R_star * R_star) / (r * r * r)) * T_star * T_star * T_star * T_star);

}


double eq_rhs(double temp, double r, double localsigma) {

	double right;
	right = (C * localsigma * kep_freq(r) * kep_freq(r) * visc(temp,r) * t_eff(temp,r,localsigma) + star_heat(r) * 2.0 * stconst) / (2.0 * stconst);
	return pow(right,0.25);

}


double calc_temp(double r, double localsigma) {

	double T,eps,a,b,x,x1;
	int it = 1,n;
	eps = 1.0e-6;					
	
	T = 10.0;
	n = 100;

	a = r + dgrid;
	b = r - dgrid;
	x = (a+b) / 2.0;
	it++;

	while(it<n) {
		if(eq_rhs(T,a,localsigma) * eq_rhs(T,x,localsigma) < 0) //checking for the signs
			b = x;   //new interval (a,x)//
		else
			a = x;   //new interval (x,b)//
			x = (a+b) / 2.0;
		if(fabs(x1-x) < eps) {
		break; 
		}
		x1 = x;
		it++;
	}
	
	if(it >= n) {
		printf("\nsolution does not exist as iterations are not sufficient");
	}
	
	T = eq_rhs(T,x1,localsigma);
	return T;	
}


void Peremtemp(double *tempvec){				/*	boundary condition for temperature	*/
  
  	tempvec[0] = tempvec[1];
  	tempvec[ngrid+1] = tempvec[ngrid];

}


void Initial_Temp(double *tempvec, double *rvec){			/*	initial temperature profile	*/

	int i;

	tempinitdat = fopen("tempinit.dat","w");

	for(i = 1; i <= ngrid; i++) {
		tempvec[i] = calc_temp(rvec[i],sigmainit(rvec[i]));
//		tempvec[i] = 100.0 / rvec[i] ;
		fprintf(tempinitdat,"%lg   %lg  \n", rvec[i], tempvec[i]);

	}

	Peremtemp(tempvec);
	printf("temp: %lg %lg \n",tempvec[0],tempvec[ngrid+1]);

	fclose(tempinitdat);

}


void Gnuplot(const char * gnucommand){
  
	char syscommand[1024];

	printf("... A gnuplot megnyitasa es az animaciok elkeszitese ...\n\n\n");
 
	sprintf(syscommand, "echo \"%s\" | gnuplot", gnucommand);
	system(syscommand);

}



int main(){

	clock_t start = clock();

  	int i;
  	double time, sigmavec[ngrid+2], rvec[ngrid+2], tempvec[ngrid+2];
	char dens_name[64], animplot[50], gifnev[50], timechar[50], filekinev2[100], anim2plot[50], gif2nev[50];

  	time = 0.0;

	sprintf(timechar,"t=%%i");
  	surftimedat = fopen("surface.dat","w");
	sprintf(animplot,"anim.gnuplot");
	sprintf(filekinev2, "dens_%%i.dat");
	sprintf(anim2plot,"anim2.gnuplot");
	sprintf(gifnev,"surf.gif");
	sprintf(gif2nev,"temp.gif");


	for(i = 0; i <= ngrid+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = rmin + (i-1) * dgrid;
	}

	Initial_Temp(tempvec,rvec);
  	Initial_Profile(sigmavec,rvec,tempvec);


  	do {
  
  		for(i = 1; i <= ngrid; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	*/
			tempvec[i] = calc_temp(rvec[i],sigmavec[i]);
    			sigmavec[i] = (sigmavec[i] * visc(tempvec[i],rvec[i]) + dt * (Coeff_1(tempvec[i],rvec[i]) * (sigmavec[i+1] * visc(tempvec[i+1],rvec[i+1]) - 2.0 * sigmavec[i] * visc(tempvec[i],rvec[i]) + sigmavec[i-1] * visc(tempvec[i-1],rvec[i-1])) / (dgrid * dgrid) + Coeff_2(tempvec[i],rvec[i]) * (sigmavec[i+1] * visc(tempvec[i+1],rvec[i+1]) - sigmavec[i-1] * visc(tempvec[i-1],rvec[i-1])) / (2.0 * dgrid))) / visc(tempvec[i],rvec[i]);
  		}

  		Perem(sigmavec,tempvec);						/*	loading boundary condition in each timestep	*/
  		time = time + dt;

		if(fmod(time, (int)(timemax/WO)) < dt || time == 0) {		/*	a futasi ido alatt 10-szer irja ki, hogy hol tart	*/
			snprintf(dens_name,64,"dens_%i.dat", (int)(time));

			printf("%lg \n",time);

			filgnu=fopen(animplot,"w");
			fprintf(filgnu,"filename(n) = sprintf('%s', n) \n set xlabel 'Distance [AU]' \n set ylabel 'Surface density [M_Sun/AU/AU]' \n set title 'Surface density profile' \n plot filename(i) using 1:2 title sprintf('%s nap',i) with line  \n i=i+%i \n if (i < n) reread \n ",filekinev2,timechar, (int)(timemax/WO)); 

			fclose(filgnu);

			filgnu2=fopen(anim2plot,"w");
			fprintf(filgnu2,"filename(n) = sprintf('%s', n) \n set xlabel 'Distance [AU]' \n set ylabel 'Temperature [K]' \n set title 'Temperature profile' \n plot filename(i) using 1:3 title sprintf('%s nap',i) with line  \n i=i+%i \n if (i < n) reread \n ",filekinev2,timechar, (int)(timemax/WO)); 

			fclose(filgnu2);

			filplot=fopen("plot.surf","w");

/*	Animacio letrehozasa az elkeszult file-okbok. A szamozott file-okat egy for ciklussal nyitja meg a gnuplot.	*/		
			fprintf(filplot,"set term gif giant animate size 450, 450 \n set output '%s' \n n=%i \n i=%i \n load '%s' \n", gifnev, (int)timemax, (int)(timemax/WO), animplot); 

			fclose(filplot);

			filplot2=fopen("plot.temp","w");

/*	Animacio letrehozasa az elkeszult file-okbok. A szamozott file-okat egy for ciklussal nyitja meg a gnuplot.	*/		
			fprintf(filplot2,"set term gif giant animate size 450, 450 \n set output '%s' \n n=%i \n i=%i \n load '%s' \n", gif2nev, (int)timemax, (int)(timemax/WO), anim2plot); 

			fclose(filplot2);



			f_dens=fopen(dens_name,"w");

  			for(i = 1; i <= ngrid; i++) {							
				fprintf(f_dens,"%lg %lg %lg\n",rvec[i], sigmavec[i],tempvec[i]);
			}

			fclose(f_dens);
		}


  	} while(time < timemax);

 	for(i = 1; i <= ngrid; i++) {
  
   		fprintf(surftimedat, "%lg   %lg\n", rvec[i], sigmavec[i]); 		/*	print sigma in a file in tmax	*/
   
 	}

	surftimedat = fopen("surftime.dat","w");
	temptimedat = fopen("temptime.dat","w");
  
 	for(i = 1; i <= ngrid; i++) {

  		fprintf(temptimedat, "%lg   %lg\n", rvec[i], tempvec[i]); 		
  		fprintf(surftimedat, "%lg   %lg\n", rvec[i], sigmavec[i]);   

 	}

	fclose(temptimedat);
	fclose(surftimedat);

	Gnuplot("load 'plot.surf'");	
	Gnuplot("load 'plot.temp'");	

	clock_t stop = clock();
   	printf("Finished in about %Lg seconds. \n", (long double)(stop - start) / CLOCKS_PER_SEC);

	return 0;

}

