#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define	G		1.0
//#define G 		0.01720209895			/*	Gaussian-grav const in astron. units		*/
#define G2		G*G
#define AU2CM		1.496e13			/*	1 AU --> 1 cm					*/
#define asp_ratio 	0.05				/*	aspect ratio					*/

#define MMOL		3.9e-24				/*	mean molecular weight in g-s? (mu*m_hydrogen)	*/
#define SIGMOL		2.0e-15				/*	collision cross-section in cm*cm		*/
#define PDENSITY    	3.0
#define CBODY		1.0				/*	central body mass in solar units		*/
#define GCM2MSAU	1.125211e-7			/*	g/cm/cm --> M_Sun/AU/AU				*/
#define M_SOL		1.9891e33			/*	mass of Sun in gramms				*/

#define WO		100.0

#define filenev		"dens_0.dat"
#define filenev2	"distr.dat"

FILE *filebe, *fildrag, *out, *histfil;



int Sorokszama(int numout){

	char c;
	filebe = fopen(filenev,"r+");		
	numout = 0;

/*	A feluletisuruseg profilt tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal	*/

		while((c = fgetc(filebe)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 					*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert								*/
	fclose(filebe);	
	return numout;

}


int Reszecskekszama(int numout){

	char c;
	filebe = fopen(filenev2,"r+");		
	numout = 0;

/*	A random generalt reszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal	*/

		while((c = fgetc(filebe)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 						*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert									*/
	fclose(filebe);				
	return numout;

}


/*	A feluletisuruseg profilt tartalmazo file tartalmanak beolvasasa: r a sugar, sigma a feluletisuruseg, p a nyomas, dp a nyomas derivaltja	*/
void load_r_sigma_p_dp(int num, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *DD, double *RMIN, double *RMAX) {
		
	int i;

	filebe = fopen(filenev,"r+");

	for (i = 1; i <= num; i++){ 
/*	file-bol beolvasas: r [AU], sigma [M_SOL / AU / AU], homerseklet [K], p [M_SUN / (4*PI*PI) / AU], dp [M_SUN / (4*PI*PI) / AU / AU]	*/
		fscanf(filebe, "%lf %lf %lf %lf",&rvec[i],&sigmavec[i],&pressvec[i],&dpressvec[i]); 
	}
	fclose(filebe);

	rvec[0] = rvec[1];
	rvec[num+1] = rvec[num];
	sigmavec[0] = sigmavec[1];
	sigmavec[num+1] = sigmavec[num];
	pressvec[0] = pressvec[1];
	pressvec[num+1] = pressvec[num];
	dpressvec[0] = dpressvec[1];
	dpressvec[num+1] = dpressvec[num];

	*DD = rvec[2] - rvec[1];		/*	lepeskoz	*/
	*RMIN = rvec[0];			/*	"rmin"		*/
	*RMAX = rvec[num];			/*	"rmax"		*/
}


/*	a reszecske adatainak betoltese: rad a sugara (cm-ben!!!), x a helyvektor koordinatai (AU)	*/
void load_rad_xyz(int partnum, double *body_rad, double *x) {
		
	int i;

	filebe = fopen(filenev2,"r+");

	for (i = 0; i < partnum; i++){ 
/*	file-bol beolvasas: body_rad [cm], x, y [AU]	*/
		fscanf(filebe, "%lf %lf",&body_rad[i],&x[i]); 
	}
	fclose(filebe);
}


void interpol(double RMIN, double DD, double r, double *rvec,  double *sigmavec, double *pressvec, double *dpressvec, double *sigma_interp, double *press_interp, double *dpress_interp) {

	double coef1, coef2, coef3, rmid, rindex;
	int index;

    	rmid = (r - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
   	index = (int) floor(rmid + 0.5);					/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
    	rindex = rvec[index];       						/* 	The corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

   	coef1 = (sigmavec[index + 1] - sigmavec[index]) / DD;      		/*	Ez az alabbi ket sor a linearis interpolacio 			*/
     	*sigma_interp = sigmavec[index] + coef1 * (r - rindex);          	/*	M_Sun / AU / AU (cgs-ben kell majs hasznalni)			*/

     	coef2 = (dpressvec[index + 1] - dpressvec[index]) / DD; 		/*	Ugyanaz az interpolacio a nyomas derivaltjara			*/
     	*dpress_interp  = dpressvec[index] + coef2 * (r - rindex);		/*	[M_SUN / (yr * yr / 4*PI*PI) / AU / AU]					*/
	
	coef3 = (pressvec[index + 1] - pressvec[index]) / DD;			/*	Ugyanaz az interpolacio a nyomasra				*/
	*press_interp = pressvec[index] + coef3 * (r - rindex);			/*	[M_SUN / (yr * yr / 4*PI*PI) / AU]					*/
	

}


/*	Keplerian velocity of the particles	*/
double kepler_vel(double r) {

	return sqrt(G2 * CBODY / r);		/*	v_kepler in AU / yr / 2*PI	*/

}


/*	density in the mid-plane in M_sol / AU / AU / AU	*/
double rho_midpl(double r, double sigma) {

	return sigma/(sqrt(2.0 * M_PI) * asp_ratio * r);

}


/*	mean-path of the molecules (in cm)	*/
double mean_free_path(double r, double sigma) {

	double rho;
	rho = rho_midpl(r,sigma);				
	return MMOL / SIGMOL / (rho * M_SOL / AU2CM / AU2CM / AU2CM);		/*	lambda in cm-s	*/

}


/*	Local speed of sound in the dimension of v_kepler	*/
double speed_of_sound(double r) {

	double v_kep;
	v_kep = kepler_vel(r);				/*	AU / day / 2*PI		*/
	return v_kep * asp_ratio;			/*	c_s = H * kep_freq = h * r * kep_freq = h * r * v_kep / r = h * v_kep	*/

}


/*	Mach-number: dimensionless if delv and csound are in the same dimension		*/
double Mach(double r, double delv) {

	double cs;
	cs = speed_of_sound(r);
	return delv/cs;

}


/*	Reynolds number: dimensionless if body_radius and lambda in cm	*/
double Reynolds(double r , double sigma, double delv, double body_radius) {

	double cs,lambda,reyn;
	cs = speed_of_sound(r);							/*	in AU / day		*/
	lambda = mean_free_path(r,sigma);					/*	mean free-path in cm	*/
	reyn = 6.0 * body_radius * delv / (sqrt(8.0 / M_PI) * cs * lambda);	/*	delv in AU/day		*/
	return reyn;	

}


/*	Knudsen number: dimensionless if lambda and body_radius in cm	*/
double Knudsen(double r, double sigma, double body_radius) {

	double lam;
	lam = mean_free_path(r,sigma);
	return 0.5 * lam / body_radius;

}


void Eps_Stokes_coeff(double mach, double ren, double *ceps, double *cstokes) {

	double cstk;

	*ceps = 2.0 * sqrt(1.0 + 128.0 / (9.0 * M_PI * mach * mach));		/*	dimensionless	*/

	if (ren <= 500.0) {
		cstk = 24.0 / ren + 3.6 * pow(ren,-0.313);		
	} else if (ren <= 1500.0) {
		cstk = 9.5e-5 * pow(ren,1.397);
	} else {
		cstk = 2.61;
	}

	*cstokes = cstk;							/*	dimensionless	*/
	
}


/*	Asimuthal velocity of the gas. 
   If dp/dr < 0 --> v(phi),gas < v_kepler    ---> "drag force" is negative 	(gas slows down the particles --> inward movement)
   If dp/dr = 0  --> v(phi),gas = v_kepler   ---> "drag force" is zero 		(gas orbits with keplerian velocity)
   If dp/dr > 0  --> v(phi),gas > v_kepler   ---> "drag force" is positive 	(gas accelerates the particles --> outward movement) 	*/
void v_phi_gas(double r, double sigma, double x, double dp, double *dvx) {

	double rho,v_kep,v_g2,v;
	
	rho = rho_midpl(r,sigma);			/*	rho in M_SUN / AU / AU / AU					*/
	v_kep = kepler_vel(r); 				/*	v_kepler in AU / day						*/
	v_g2 = v_kep * v_kep + r * dp / rho;		/*	v(phi),gas^2 / r = G * M_star / r^2 + 1/rho*dp/dr (in cm/s) 	*/

	v = sqrt(v_g2); 				/* 	Azimuthal velocity of the gas in AU / day			*/
	
	*dvx = v_kep - v;					/*	relative velocity of gas and the particle in AU / day		*/
	
}


/*	get the "acceleration" from gas: r in AU, sigma in M_SUN / AU / AU, dvx in AU / day, body_radius in cm: gives adrag in	*/
void get_drag_acc(double r, double sigma, double dvx, double body_radius, double *adrag) {

	double machn,ren,cdrag,knud,delv,ceps,cstokes,rhogas;

	delv = dvx;								/*	relative velocity of gas and the particle in AU / day		*/

	machn = Mach(r,delv);							/*	dimensionless if delv in AU / day				*/
	ren = Reynolds(r,sigma,delv,body_radius);				/*	dimensionless if body_rad in cm, r in AU, delv in AU/day	*/
	knud = Knudsen(r,sigma,body_radius);					/*	dimensionless if body_radius in cm				*/	
	rhogas = rho_midpl(r,sigma) * M_SOL / AU2CM / AU2CM / AU2CM;		/*	density in g/cm/cm/cm						*/

	Eps_Stokes_coeff(machn,ren,&ceps,&cstokes);				/*	dimensionless Epstein & Stokes coeff				*/

	cdrag = (9.0 * knud * knud * ceps + cstokes) / ((3.0 * knud + 1.0) * (3.0 * knud + 1.0));	/*	dimensionless "average" coefficient	*/

	*adrag = -3.0 * rhogas * cdrag * delv / (8.0 * PDENSITY * body_radius / AU2CM);			/*	1/day					*/
	
}


/*	friction time of particles in days if adrag is in 1/day	*/
double fric_time(double r, double x, double delv, double adrag, double mass) {
	
	double drag_acc;
	drag_acc = - G2 * (CBODY + mass) * x / r / r / r + adrag * delv;	/*	acceleration from dragforce of gas in AU / day / day		*/

//	printf("tau_fric: %lg sec\n",fabs(delv)/fabs(drag_acc) * 24.0 * 3600.0);
	
	return  fabs(delv)/fabs(drag_acc);								/*	in days								*/

}


/*	dimensionless if tau_fric is in days and omega_kep is in 1/day	*/
double Stokes_param(double r, double x, double delv, double adrag, double mass){

	double tau_fric,kep_freq;
	tau_fric = fric_time(r,x,delv,adrag,mass);			/*	in days								*/
	kep_freq = kepler_vel(r) / r;					/*	Kepler angular frequency in 1 / day				*/
	return tau_fric * kep_freq;					/*	dimensionless if tau_fric is in days and omega_kep is in 1/day	*/

}

/*	derivates for 4th order Runge-Kutta: x' = st/(1+st*st) * H/r * r/p * dp/dr * c_sound	*/
void eqrhs(double r, double sigma, double x, double press, double dpress, double delv, double body_radius, double mass, double *der_x) {

	double adrag_i,temp_vx,cs,stokes;

	get_drag_acc(r,sigma,delv,body_radius,&adrag_i);		/*	adrag_i in 1/day					*/
	stokes = Stokes_param(r,x,delv,adrag_i,mass);			/*	dimensionless if tau_fric in days and adrag in 1/day	*/
	cs = speed_of_sound(r);						/*	in AU / day						*/

	temp_vx = stokes / (1.0 + stokes * stokes) * asp_ratio * r / press * dpress * cs; 	/*	in AU / day			*/

//	printf("temp_vx: %lg\n",temp_vx);

	*der_x = temp_vx;

}

/*
void Kunge_Kutta(double r, double sigma, double x, double press, double dpress,  double delv, double body_radius, double mass, double step, double *x_new) {

	double xtemp, dy1, dy2, dy3, dy4, xnew;

	eqrhs(r,sigma,x,press,dpress,delv,body_radius,mass,&dy1);
	xtemp = x + 0.5 * step * dy1;

	eqrhs(r,sigma,xtemp,press,dpress,delv,body_radius,mass,&dy2);
	xtemp = x + 0.5 * step * dy2;

	eqrhs(r,sigma,xtemp,press,dpress,delv,body_radius,mass,&dy3);
	xtemp = x + step * dy3;

	eqrhs(r,sigma,xtemp,press,dpress,delv,body_radius,mass,&dy4);
	xnew = x + step * (dy1 + 2.0 * dy2 * 2.0 * dy3 * dy4) / 6.0;

	*x_new = xnew;

}
*/


void Euler(double r, double sigma, double x, double press, double dpress,  double delv, double body_radius, double mass, double step, double *x_new) {

	double xtemp, dy1;

	eqrhs(r,sigma,x,press,dpress,delv,body_radius,mass,&dy1);
	xtemp = x + step * dy1;

	*x_new = xtemp;

}



/*
void histogram(int num, double r, int *hist, double RMAX) {

	int index;
	double rmid, hist_i;
    	rmid = r / RMAX * num;     						/* 	The integer part of this gives at which index is the body	*/
//   	index = (int) floor(rmid + 0.5);					/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
//	hist_i = hist[index];
//    	hist[index] = hist_i + 1;       					/* 	The corresponding r, e.g rd[ind] < r < rd[ind+1]		*/
/*
}
*/


int main() {

	int i,num,linesout,partnum;
	linesout = 0;
	num = 0;

	double t, deltat,t_integration;

	t = 0.0;
	deltat = 0.1;
	t_integration = 1000000.0 / 2.0 * M_PI;

	num = Sorokszama(linesout);							/*	A feluletisuruseg profilt tartalmazo file sorainak szama elmentve egy integerbe	*/	
	partnum = Reszecskekszama(linesout); 						/*	A reszecskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/

	double body_radius[partnum],x[partnum],vx[partnum];
	double rvec[num],sigmavec[num],pressvec[num],dpressvec[num];
	double DD,RMIN,RMAX;
	double rvec_part[partnum];
	double sigma_interp, press_interp, dpress_interp;
	double sigma_interp_vec[partnum], press_interp_vec[partnum], dpress_interp_vec[partnum];
	double dvx, dvx_vec[partnum];
	double x_new[partnum];
	double mass[partnum];								/*	mass of the particles	*/
	
	load_r_sigma_p_dp(num,rvec,sigmavec,pressvec,dpressvec,&DD,&RMIN,&RMAX);	/*	A feluletisuruseg es a nyomasderivalt profil betoltese es eltarolasa	*/
	load_rad_xyz(partnum,body_radius,x);						/*	A porreszecskek adatainak betoltese es eltarolasa			*/

	for(i = 0; i < partnum; i++) {

		rvec_part[i] = x[i];
		vx[i] = kepler_vel(rvec_part[i]);
		interpol(RMIN,DD,rvec_part[i],rvec,sigmavec,pressvec,dpressvec,&sigma_interp,&press_interp,&dpress_interp);
		sigma_interp_vec[i] = sigma_interp;									/*	M_sun / AU / AU 								*/
		press_interp_vec[i] = press_interp;					/*	M_SUN / (yr * yr / 4*PI*PI) / AU   --> M_SUN / AU /  day / day			*/
		dpress_interp_vec[i] = dpress_interp;					/*	M_Sun / AU / AU / AU / (yr * yr / 4*PI*PI) -->  M_SUN / AU / AU /  day / day	*/
		v_phi_gas(rvec_part[i],sigma_interp_vec[i],x[i],dpress_interp_vec[i],&dvx);
		dvx_vec[i] = dvx;
		mass[i] = PDENSITY * 4.0 * M_PI / 3.0 * body_radius[i] * body_radius[i] * body_radius[i] / M_SOL;	/*	mass of the particle in solar mass	*/
		printf("x: %lg sigma: %lg  p: %lg  dp: %lg delv: %lg  mass: %lg \n",x[i],  sigma_interp_vec[i], press_interp_vec[i], dpress_interp_vec[i], dvx_vec[i], mass[i]);

	}



	do {

//		for (i = 0; i < num; i++) hist[i] = 0;

		for (i = 0; i < partnum; i++) {
			rvec_part[i] = x[i];
						
			if(rvec_part[i] >= RMIN && rvec_part[i] <= RMAX)	{
				
				Euler(rvec_part[i],sigma_interp_vec[i],x[i],press_interp_vec[i],dpress_interp_vec[i],dvx_vec[i],body_radius[i],mass[i],deltat,&x_new[i]);

				x[i] = x_new[i];

				interpol(RMIN,DD,rvec_part[i],rvec,sigmavec,pressvec,dpressvec,&sigma_interp,&press_interp,&dpress_interp);
				sigma_interp_vec[i] = sigma_interp;								
				press_interp_vec[i] = press_interp;
				dpress_interp_vec[i] = dpress_interp;			
				v_phi_gas(rvec_part[i],sigma_interp_vec[i],x[i],dpress_interp_vec[i],&dvx);
				dvx_vec[i] = dvx;

				vx[i] = kepler_vel(rvec_part[i] = x[i]);
				printf("%i time: %lg r: %lg  dvx: %e\n",i,t/(2.0 * M_PI),rvec_part[i],dvx_vec[i]);
//				histogram(num,rvec_part[i],hist,RMAX);	

			} else {

				x[i] = 0.0;
				vx[i] = 0.0;

			}		

		}
		
	

/*		if(fmod(t, (int)(t_integration/WO)) < deltat || t == 0) {		/*	a futasi ido alatt WO-szor irja ki, hogy hol tart	*/
/*			snprintf(hist_name,64,"hist_%i.dat", (int)(t));
			snprintf(filname,64,"out_%i.dat",(int)(t));
			histfil = fopen(hist_name,"w");
			out = fopen(filname,"w"); 
	
			printf("%lg\n",t);

  			for(i = 0; i < num; i++) {							
				fprintf(histfil,"%lg %i\n",rvec[i], hist[i]);		/*	file for histogram to find out the flux of particles	*/
/*			}

			fclose(histfil);

			for(i = 0; i < partnum; i++) {
				fprintf(out,"%lg %lg %lg %lg %lg\n",rvec_part[i],v_kep_vec[i],sq_v_phi,sigma_interp,dpress_interp);
			}
			fclose(out);
		}
*/
		t = t + deltat;

	} while (t < t_integration);	



	
		
	return 0;

}
