#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//#define	G		1.0
#define G 		0.01720209895			/*	Gaussian-grav const in astron. units		*/
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

#define filenev		"dens_47476.dat"
#define filenev2	"distr.dat"

#define PARTGASRATIO	0.01

FILE *filebe, *fildrag, *out, *histfil, *mass_growth;



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

	for (i = 0; i < num; i++){ 
/*	file-bol beolvasas: r [AU], sigma [M_SOL / AU / AU], p [M_SUN / day / day / AU], dp [M_SUN / day / day / AU / AU]	*/
		fscanf(filebe, "%lf %lf %lf %lf",&rvec[i],&sigmavec[i],&pressvec[i],&dpressvec[i]); 
	}
	fclose(filebe);
/*
	rvec[0] = rvec[1];
	rvec[num+1] = rvec[num];
	sigmavec[0] = sigmavec[1];
	sigmavec[num+1] = sigmavec[num];
	pressvec[0] = pressvec[1];
	pressvec[num+1] = pressvec[num];
	dpressvec[0] = dpressvec[1];
	dpressvec[num+1] = dpressvec[num];
*/
	*DD = rvec[1] - rvec[0];		/*	lepeskoz	*/
	*RMIN = rvec[0];			/*	"rmin"		*/
	*RMAX = rvec[num-1];			/*	"rmax"		*/

}


/*	a reszecske adatainak betoltese: rad a sugara (cm-ben!!!), x a helyvektor koordinatai (AU)	*/
void load_rad_xyz(int partnum, double *body_rad, double *x, double *y) {
		
	int i;
	double temp;

	filebe = fopen(filenev2,"r+");

	for (i = 0; i < partnum; i++){ 
/*	file-bol beolvasas: body_rad [cm], x, y [AU]	*/
		fscanf(filebe, "%lf %lf %lf",&body_rad[i],&x[i], &y[i]); 
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
     	*dpress_interp  = dpressvec[index] + coef2 * (r - rindex);		/*	[M_SUN / day / day / AU / AU]					*/
	
	coef3 = (pressvec[index + 1] - pressvec[index]) / DD;			/*	Ugyanaz az interpolacio a nyomasra				*/
	*press_interp = pressvec[index] + coef3 * (r - rindex);			/*	[M_SUN / day / day / AU]					*/
	

}


/*	Keplerian velocity of the particles	*/
double kepler_vel(double r) {

	return sqrt(G2 * CBODY / r);		/*	v_kepler in AU / day	*/

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
	v_kep = kepler_vel(r);				/*	AU / day 		*/
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
void v_phi_gas(double r, double sigma, double x, double y, double vx, double vy, double dp, double *dvx, double *dvy) {

	double rho,v_kep,v_g2,v,ux, uy,phi,truan;

	rho = rho_midpl(r,sigma);
	v_kep = kepler_vel(r); 				/*	v_kepler in AU / day						*/
	v_g2 = v_kep * v_kep + r * dp / rho;		/*	v(phi),gas^2 / r = G * M_star / r^2 + 1/rho*dp/dr (in cm/s) 	*/

//	v = 0.996*v_kep; 				/* 	Azimuthal velocity of the gas in AU / day			*/

	v = sqrt(v_g2);

//	printf("vkep: %lg vx: %lg vy: %lg x: %lg  y: %lg",v_kep,vx,vy,x,y);


	truan = atan2(y,x);
	
	
	ux =-v*y/r; // a gaz sebessege  
	uy = v*x/r;
	

//	phi = acos(x / r);
//	ux = -sin(phi) * v; // a gaz sebessege  
//	uy = cos(phi) * v;

//	printf("",);

	*dvx = vx - ux;					/*	relative velocity of gas and the particle in AU / day		*/
	*dvy = vy - uy;
	
}


/*	get the "acceleration" from gas: r in AU, sigma in M_SUN / AU / AU, dvx in AU / day, body_radius in cm: gives adrag in	*/
void get_drag_acc(double r, double sigma, double delv, double body_radius, double *adrag) {

	double machn,ren,cdrag,knud,ceps,cstokes,rhogas;

	machn = Mach(r,delv);							/*	dimensionless if delv in AU / day				*/
	ren = Reynolds(r,sigma,delv,body_radius);				/*	dimensionless if body_rad in cm, r in AU, delv in AU/day	*/
	knud = Knudsen(r,sigma,body_radius);					/*	dimensionless if body_radius in cm				*/	
	rhogas = rho_midpl(r,sigma) * M_SOL / AU2CM / AU2CM / AU2CM;		/*	density in g/cm/cm/cm						*/

	Eps_Stokes_coeff(machn,ren,&ceps,&cstokes);				/*	dimensionless Epstein & Stokes coeff				*/

	cdrag = (9.0 * knud * knud * ceps + cstokes) / ((3.0 * knud + 1.0) * (3.0 * knud + 1.0));	/*	dimensionless "average" coefficient	*/

	*adrag = -3.0 * rhogas * cdrag * delv / (8.0 * PDENSITY * body_radius / AU2CM);	/*	1/day					*/
	
}


/*	derivates for 4th order Runge-Kutta: x' = vx, y' = vy, z' = vz, vx' = ax, vy' = ay, vz' = az	*/
void eqrhs(double r, double sigma, double *y, double press, double dpress, double dvx, double dvy, double body_radius, double mass, double *dy) {

	double adrag_i,r3,delv;

	delv = sqrt(dvx * dvx + dvy * dvy);

	get_drag_acc(r,sigma,delv,body_radius,&adrag_i);		/*	adrag_i in 1/day					*/

	r3 = r * r * r;

	dy[0] = y[2];
	dy[1] = y[3];
	dy[2] = -G2 * (CBODY + mass) * y[0] / r3 + adrag_i * dvx;
	dy[3] = -G2 * (CBODY + mass) * y[1] / r3 + adrag_i * dvy;

}


void Runge_Kutta(double r, double sigma, double x, double y, double vx, double vy, double press, double dpress,  double dvx, double dvy, double body_radius, double mass, double step, double *x_new, double *y_new, double *vx_new, double *vy_new) {

	double temp[4], dy1[4], posvel[4], dy2[4], dy3[4], dy4[4], newpos[4];
	int i;

	posvel[0] = x;
	posvel[1] = y;
	posvel[2] = vx;
	posvel[3] = vy;

	eqrhs(r,sigma,posvel,press,dpress,dvx,dvy,body_radius,mass,dy1);
	for(i = 0; i < 4; i++) temp[i] = posvel[i] + 0.5 * step * dy1[i];

	eqrhs(r,sigma,temp,press,dpress,dvx,dvy,body_radius,mass,dy2);
	for(i = 0; i < 4; i++) temp[i] = posvel[i] + 0.5 * step * dy2[i];

	eqrhs(r,sigma,temp,press,dpress,dvx,dvy,body_radius,mass,dy3);
	for(i = 0; i < 4; i++) temp[i] = posvel[i] + step * dy3[i];

	eqrhs(r,sigma,temp,press,dpress,dvx,dvy,body_radius,mass,dy4);
	for(i = 0; i < 4; i++) newpos[i] = posvel[i] + step * (dy1[i] + 2.0 * dy2[i] + 2.0 * dy3[i] + dy4[i]) / 6.0;

	*x_new = temp[0];
	*y_new = temp[1];
	*vx_new = temp[2];
	*vy_new = temp[3];

}


/*
void Euler(double r, double sigma, double x, double y, double vx, double vy, double press, double dpress,  double dvx, double dvy, double body_radius, double mass, double step, double *x_new, double *y_new, double *vx_new, double *vy_new) {

	double temp[4], dy1[4], posvel[4];
	int i;

	posvel[0] = x;
	posvel[1] = y;
	posvel[2] = vx;
	posvel[3] = vy;

	eqrhs(r,sigma,posvel,press,dpress,dvx,dvy,body_radius,mass,dy1);
	for(i = 0; i < 4; i++) temp[i] = posvel[i] + step * dy1[i];

	*x_new = temp[0];
	*y_new = temp[1];
	*vx_new = temp[2];
	*vy_new = temp[3];

}
*/


/*	mass of gas at a given distance	*/
double gas_mass(double r, double sigma) {

	return M_PI * r * sigma;

}


/*	mass of the particles at a given distance	*/
double part_mass(double r, double sigma) {

	return gas_mass(r,sigma) * PARTGASRATIO;

}

void histogram(int num, double r, int *hist, double RMAX) {

	int index;
	double rmid, hist_i;
    	rmid = r / RMAX * num;     						/* 	The integer part of this gives at which index is the body	*/
   	index = (int) floor(rmid + 0.5);					/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
	hist_i = hist[index];
    	hist[index] = hist_i + 1;       					/* 	The corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

	

}


int main() {

	int i,num,linesout,partnum;
	linesout = 0;
	num = 0;

	double t, deltat,t_integration;

	t = 0.0;
	deltat = 0.1;
	t_integration = 100.0 * 365.2422;

	num = Sorokszama(linesout);							/*	A feluletisuruseg profilt tartalmazo file sorainak szama elmentve egy integerbe	*/	
	partnum = Reszecskekszama(linesout); 						/*	A reszecskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/

	double body_radius[partnum],x[partnum],vx[partnum],y[partnum],vy[partnum];
	double rvec[num],sigmavec[num],pressvec[num],dpressvec[num];
	double DD,RMIN,RMAX;
	double rvec_part[partnum];
	double sigma_interp, press_interp, dpress_interp;
	double sigma_interp_vec[partnum], press_interp_vec[partnum], dpress_interp_vec[partnum];
	double dvx, dvx_vec[partnum], dvy, dvy_vec[partnum];
	double x_new[partnum], vx_new[partnum],y_new[partnum], vy_new[partnum];
	double mass[partnum];								/*	mass of the particles	*/
	double v_kep,phi;
	double mass_particle_in, mass_particle_out;
	int num_particle_in, num_particle_out;
	double ind_ii, ind_io, ind_oi, ind_oo;
	double rmid;
	int ind, ind_min, ind_max;


	int hist[num];

	char filname[64], hist_name[64], mass_name[64];
	
	load_r_sigma_p_dp(num,rvec,sigmavec,pressvec,dpressvec,&DD,&RMIN,&RMAX);	/*	A feluletisuruseg es a nyomasderivalt profil betoltese es eltarolasa	*/
	load_rad_xyz(partnum,body_radius,x,y);						/*	A porreszecskek adatainak betoltese es eltarolasa			*/

	for(i = 0; i < partnum; i++) {

		mass[i] = PDENSITY * 4.0 * M_PI / 3.0 * body_radius[i] * body_radius[i] * body_radius[i] / M_SOL;	/*	mass of the particle in solar mass	*/
		rvec_part[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
		v_kep = kepler_vel(rvec_part[i]);

//		printf("%lg %lg %lg \n",rvec_part[i],x[i],y[i]);
		
		phi = acos(x[i] / rvec_part[i]);
		vx[i] = sin(phi) * v_kep;
		vy[i] = - cos(phi) * v_kep;
		
		interpol(RMIN,DD,rvec_part[i],rvec,sigmavec,pressvec,dpressvec,&sigma_interp,&press_interp,&dpress_interp);
		sigma_interp_vec[i] = sigma_interp;									/*	M_sun / AU / AU 								*/
		press_interp_vec[i] = press_interp;					/*	M_SUN / (yr * yr / 4*PI*PI) / AU   --> M_SUN / AU /  day / day			*/
		dpress_interp_vec[i] = dpress_interp;					/*	M_Sun / AU / AU / AU / (yr * yr / 4*PI*PI) -->  M_SUN / AU / AU /  day / day	*/
		v_phi_gas(rvec_part[i],sigma_interp_vec[i],x[i],y[i],vx[i],vy[i],dpress_interp_vec[i],&dvx,&dvy);
		dvx_vec[i] = dvx;
		dvy_vec[i] = dvy;

	}

	do {

		for (i = 0; i < num; i++) hist[i] = 0;
		mass_particle_in = 0.0;
		mass_particle_out = 0.0;
		num_particle_in = 0.0;
		num_particle_out = 0.0;
		ind_ii = 0;
		ind_io = 0;
		ind_io = 0;
		ind_oo = 0;
		ind_min = 0;
		ind_max = 0;

		if(fmod(t, (int)(t_integration/WO)) < deltat || t == 0) {		/*	a futasi ido alatt WO-szor irja ki, hogy hol tart	*/
			snprintf(hist_name,64,"hist_%i.dat", (int)(t/365.2422));
			snprintf(filname,64,"out_%i.dat",(int)(t/365.2422));
			snprintf(mass_name,64,"mass_%i.dat",(int)(t/365.2422));

			histfil = fopen(hist_name,"w");
			out = fopen(filname,"w"); 
			mass_growth = fopen(mass_name,"w"); 
	
			printf("%lg  \n",t/365.2422);

			for(i = 0; i < partnum; i++) {
	
				histogram(num,rvec_part[i],hist,RMAX);			/*	counts the number of particles at a given grid	*/
		
				fprintf(out,"%lg %lg %lg %lg %lg %lg %lg\n",rvec_part[i],x[i],y[i],kepler_vel(rvec_part[i]),sqrt(kepler_vel(rvec_part[i]) * kepler_vel(rvec_part[i]) + rvec_part[i] * dpress_interp_vec[i] / rho_midpl(rvec_part[i],sigma_interp_vec[i])),sigma_interp_vec[i],dpress_interp_vec[i]);

/*	If the distance of a particle is around the pressure maximum with 0.2 AU, it counts the "mass" of the particles: counts in the 0.4 AU annulus the mass of the gas, then by 0.01 particle to gas mass ratio, it counts the mass of the particles. Then later the program counts the number of the particles in the annulus, and dividing the mass of the particles with their number, we can get a fictive mass of each particle (how mach mass they "represent")	*/
				if((rvec_part[i] >= 1.86 - 0.2) && (rvec_part[i] <= 1.86 + 0.2)) {
					mass_particle_in += part_mass(rvec_part[i],sigma_interp_vec[i]);
				}

/*	The same like earlier	*/
				if((rvec_part[i] >= 9.4923 - 0.5) && (rvec_part[i] <= 9.4923 + 0.5)) {
					mass_particle_out += part_mass(rvec_part[i],sigma_interp_vec[i]);
				}

			}

			fclose(out);


/*	Counting the number of the particles in a given annulus with the pressure maxima	*/
			for(i = 0; i < num; i++) {
				if(rvec[i] > 1.48 - DD / 2.0 && rvec[i] < 1.48 + DD / 2.0) {
				    	rmid = (rvec[i] - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
					ind_ii = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/ 
				}

				if(rvec[i] > 1.88 - DD / 2.0 && rvec[i] < 1.88 + DD / 2.0) {
				    	rmid = (rvec[i] - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
					ind_io = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
				}

				if(rvec[i] > 8.9923 - DD / 2.0 && rvec[i] < 8.9923 + DD / 2.0) {
				    	rmid = (rvec[i] - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
					ind_oi = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
				}	

				if(rvec[i] > 9.9923 - DD / 2.0 && rvec[i] < 9.9923 + DD / 2.0) {
				    	rmid = (rvec[i] - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
					ind_oo = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
				}

			}


			ind_min = (int)ind_ii;
			ind_max = (int)ind_io;

/*	Between the upper and lower limit (index) it counts the number of the particles		*/
			for (i = ind_min; i <= ind_max; i++) {
				num_particle_in = num_particle_in + hist[i];
//				printf("%i %lg\n",hist[i],rvec[i]);
//				getchar();
			}

			ind_min = (int)ind_oi;
			ind_max = (int)ind_oo;

/*	Between the upper and lower limit (index) it counts the number of the particles		*/
			for (i = ind_min; i <= ind_max; i++) {
				num_particle_out = num_particle_out + hist[i];
			}
	
			fprintf(mass_growth,"%lg %lg\n",mass_particle_in / (double)num_particle_in, mass_particle_out / (double)num_particle_out);

			fclose(mass_growth);

 			for(i = 0; i < num; i++) {							
				fprintf(histfil,"%lg %i\n",rvec[i], hist[i]);		/*	file for histogram to find out the flux of particles	*/
			}

			fclose(histfil);

		}

		for (i = 0; i < partnum; i++) {

			rvec_part[i] = sqrt(x[i] * x[i] + y[i] * y[i]);		/* 	Counts r = sqrt(x*x + y*y) for each particle	*/

/*	If RMIN < r < RMAX, it does the following:	*/
			if((rvec_part[i] > RMIN && rvec_part[i] < RMAX)) {

/*	If ((RMIN < ABS(x) < RMAX) AND 0 < ABS(y) < RMAX) OR (RMIN < ABS(y) < RMAX) AND 0 < ABS(x) < RMAX) it does the following:	*/
				if((fabs(x[i]) > RMIN && fabs(x[i]) <= RMAX && fabs(y[i]) >= 0 && fabs(y[i]) <= RMAX) || (fabs(y[i]) >= RMIN && fabs(y[i]) <= RMAX && fabs(x[i]) >= 0 && fabs(x[i]) <= RMAX)) {

					interpol(RMIN,DD,rvec_part[i],rvec,sigmavec,pressvec,dpressvec,&sigma_interp,&press_interp,&dpress_interp);	/*	interpolation of the surface density, pressure and derivate of the pressure for each particle at a given distance	*/
					sigma_interp_vec[i] = sigma_interp;										/*	Loads a vector with the interpolated surface density ...								*/
					press_interp_vec[i] = press_interp;										/*	... pressure ...													*/					
					dpress_interp_vec[i] = dpress_interp;										/*	... and the derivate of the pressure											*/
					v_phi_gas(rvec_part[i],sigma_interp_vec[i],x[i],y[i],vx[i],vy[i],dpress_interp_vec[i],&dvx,&dvy);		/*	Counts the velocity of the gas												*/
					dvx_vec[i] = dvx;												/*	Loads a vector with the x component of the deltav (vkepler - vgas) 							*/
					dvy_vec[i] = dvy;												/*	Loads a vector with the y component of the deltav (vkepler - vgas) 							*/

					Runge_Kutta(rvec_part[i],sigma_interp_vec[i],x[i],y[i],vx[i],vy[i],press_interp_vec[i],dpress_interp_vec[i],dvx_vec[i],dvy_vec[i],body_radius[i],mass[i],deltat,&x_new[i],&y_new[i],&vx_new[i],&vy_new[i]);	/*	timestepping by RK4	*/
					x[i] = x_new[i];		/*	With the RK4 method we get the new position x ...						*/
					y[i] = y_new[i];		/*	... y, and new velocity ...									*/
					vx[i] = vx_new[i];		/*	... vx and ...											*/
					vy[i] = vy_new[i];		/*	... vy components, and it loads the new components in the "old" position and velocity vectors	*/

				}

/*	If RMIN > r OR r > RMAX, it gives all the position and vector components zero, so it steps out of the "system"	*/
			} else {

				x[i] = 0.0;
				y[i] = 0.0;
				vx[i] = 0.0;
				vy[i] = 0.0;
			}	

		}
		
		t = t + deltat;		/*	timestepping	*/

	} while (t < t_integration);	


	return 0;

}
