#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define G 		0.01720209895			/*	Gaussian-grav const in astron. units		*/
#define G2		G*G
#define AU2CM		1.496e13			/*	1 AU --> 1 cm					*/
#define asp_ratio 	0.05				/*	aspect ratio					*/

#define MMOL		3.9e-24				/*	mean molecular weight ??? in cm-s?		*/
#define SIGMOL		2.0e-15				/*	collision cross-section in cm*cm		*/
#define PDENSITY    	3.0
#define CBODY		1.0				/*	central body mass in solar units		*/
#define GCM2MSAU	1.125211e-7			/*	g/cm/cm --> M_Sun/AU/AU				*/
#define	DDENSASTR2CGS	1.1906e-3			/*	deriv p: M_Sun/day/day/AU/AU --> g/s/s/m/m	*/
#define M_SOL		1.9891e33			/*	mass of Sun in gramms				*/

#define WO		100.0

#define filenev		"dens_559000.dat"
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


/*	A feluletisuruseg profilt tartalmazo file tartalmanak beolvasasa: r a sugar, sigma a feluletisuruseg, dp a nyomas derivaltja	*/
void load_r_sigma_dp(int num, double *rvec, double *sigmavec, double *dpressvec, double *DD, double *RMIN, double *RMAX) {
		
	int i;
	double temp[num];					/*	homerseklet tarolasa, egyelore nem hasznaljuk	*/

	filebe = fopen(filenev,"r+");

	for (i = 1; i < num; i++){ 
/*	file-bol beolvasas: r,sigma,homerseklet,dp	*/
		fscanf(filebe, "%lf %lf %lf %lf",&rvec[i],&sigmavec[i],&temp[i],&dpressvec[i]); 
	}
	fclose(filebe);

	rvec[0] = rvec[1];
	rvec[num] = rvec[num-1];
	sigmavec[0] = sigmavec[1];
	sigmavec[num] = sigmavec[num-1];
	dpressvec[0] = dpressvec[1];
	dpressvec[num] = dpressvec[num-1];

	*DD = rvec[2] - rvec[1];		/*	lepeskoz	*/
	*RMIN = rvec[1];			/*	"rmin"		*/
	*RMAX = rvec[num-1];			/*	"rmax"		*/
}


/*	a reszecske adatainak betoltese: rad a sugara (cm-ben!!!), x,y,z a hely-, vx,vy,vz a sebessegvektor koordinatai (AU es AU/day)	*/
void load_rad_xyz_vxvyvz(int partnum, double *rad, double *x, double *y, double *z, double *vx, double *vy, double *vz) {
		
	int i;

	filebe = fopen(filenev2,"r+");

	for (i = 0; i < partnum; i++){ 
/*	file-bol beolvasas: r,sigma,homerseklet,dp	*/
		fscanf(filebe, "%lf %lf %lf %lf %lf %lf %lf",&rad[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]); 
	}
	fclose(filebe);
}


/*	density in the mid-plane in M_sol / AU / AU / AU	*/
double rho_midpl(double r, double sigma) {

	return sigma/(sqrt(2.0 * M_PI) * asp_ratio * r);

}


/*	density "out of" the mid-plane	*/
/*double rho_gas(double r, double sigma, double z) {

	double rho_0, rho;
	rho_0 = rho_midpl(r,sigma);
	rho = rho_0 * exp(- z * z / (2.0 * asp_ratio * asp_ratio * r * r));
	return rho;

}
*/

/*	mean-path of the molecules	*/
double mean_free_path(double r, double sigma, double z) {

	double rho;
//	rho = rho_gas(r,sigma,z);
	rho = rho_midpl(r,sigma);				
	return MMOL / SIGMOL / (rho * M_SOL / AU2CM / AU2CM / AU2CM);		/*	lambda in cm-s	*/

}


/*	Keplerian angular velocity of the particles	*/
double v_kepler(double r) {

	double v_kep;
	v_kep = sqrt(G2 * CBODY / r);			/*	v_kepler in AU / day	*/
	return v_kep;		

}


/*	Local speed of sound	*/
double speed_of_sound(double r) {

	double v_kep;
	v_kep = v_kepler(r);				/*	AU / day		*/
	return v_kep * asp_ratio;

}


void interpol(double RMIN, double DD, double r, double *sigmavec, double *dpressvec, double *rvec, double *sigma_interp, double *dpress_interp) {

	double coef1, coef2, rmid, rindex;
	int index;

    	rmid = (r - RMIN)/DD;     						/* 	The integer part of this gives at which index is the body	*/
   	index = (int) floor(rmid + 0.5);					/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
    	rindex = rvec[index];       						/* 	The corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

   	coef1 = (sigmavec[index + 1] - sigmavec[index]) / DD;      		/*	Ez az alabbi ket sor a linearis interpolacio 			*/
     	*sigma_interp = sigmavec[index] + coef1 * (r - rindex);          	/*	M_Sun / AU / AU (cgs-ben kell majs hasznalni)	*/

     	coef2 = (dpressvec[index + 1] - dpressvec[index]) / DD; 		/*	Ugyanaz az interpolacio a nyomas derivaltjara			*/
//	printf("%lg %lg %lg\n",r,dpressvec[index],dpressvec[index] + coef2 * (r - rindex));
     	*dpress_interp  = dpressvec[index] + coef2 * (r - rindex);

}


/*	Mach-number: dimensionless if delv and csound in AU/day	*/
double Mach(double r, double delv) {

	double cs;
	cs = speed_of_sound(r);
	return delv/cs;

}


/*	Reynolds number: dimensionless if body_radius and lambda in cm	*/
double Reynolds(double body_radius, double r, double delv, double z, double sigma) {

	double cs,lambda,reyn;
	cs = speed_of_sound(r);							/*	in AU / day		*/
	lambda = mean_free_path(r,sigma,z);					/*	mean free-path in cm	*/
	reyn = 6.0 * body_radius * delv / (sqrt(8.0 / M_PI) * cs * lambda);	/*	delv in AU/day		*/
	return reyn;	

}


/*	Knudsen number	*/
double Knudsen(double r, double sigma, double z, double body_radius) {

	double lam;
	lam = mean_free_path(r,sigma,z);
	return 0.5 * lam / body_radius;

}


void Eps_Stokes_coeff(double mach, double ren, double *ceps, double *cstokes) {

	double cstk;

	*ceps = 2.0 * sqrt(1.0 + 128.0 / (9.0 * M_PI * mach * mach));

	if (ren <= 500.0) {
		cstk = 24.0 / ren + 3.6 * pow(ren,-0.313);
	} else if (ren <= 1500.0) {
		cstk = 9.5e-5 * pow(ren,1.397);
	} else {
		cstk = 2.61;
	}

	*cstokes = cstk;
	
}


/*	Angular velocity of the gas. 
   If dp/dr < 0 --> v(phi),gas < v_kepler    ---> "drag force" is negative 	(gas slows down the particles --> inward movement)
   If dp/dr = 0  --> v(phi),gas = v_kepler   ---> "drag force" is zero 		(gas orbits with keplerian velocity)
   If dp/dr > 0  --> v(phi),gas > v_kepler   ---> "drag force" is positive 	(gas accelerates the particles --> outward movement) 	*/
void v_gas(double r, double sigma, double x, double y, double z, double vx, double vy, double vz, double dp, double *ux, double *uy, double *uz, double *dvx, double *dvy, double *dvz) {

	double rho,v_g2,v_kep,v,sigma_cgs,truan,vgx,vgy,vgz;
	
//	sigma_cgs = sigma / GCM2MSAU;
	rho = rho_midpl(r,sigma) * M_SOL / AU2CM / AU2CM / AU2CM;	/*	rho in g/cm/cm/cm		*/
	v_kep = v_kepler(r) * AU2CM / 86400.0;				/*	v_kepler in cm/s		*/
	v_g2 = v_kep * v_kep + r * AU2CM * dp  / rho;			/*	v(phi),gas^2 / r = G * M_star / r^2 + 1/rho*dp/dr	*/

	truan = atan2(y,x);
	
	v = sqrt(v_g2) * (86400.0 / AU2CM);				/* 	Azimuthal velocity of the gas --> conversion back to AU / DAY		*/
//	printf("%lg %lg %lg \n",r,v_kep , r * dp / rho);
	vgx = -v * sin(truan);
	vgy = v * cos(truan);
	vgz = 0.0;
	
	*dvx = vx - vgx;
	*dvy = vy - vgy;
	*dvz = vz - vgz;

	*ux = vgx;
	*uy = vgy;
	*uz = vgz;
	
}


void get_drag_acc(double r, double sigma, double z, double body_radius, double dvx, double dvy, double dvz, double *adrag) {

	double machn,ren,cdrag,knud,delv,ceps,cstokes,rhogas;

	delv = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

	machn = Mach(r,delv);
	ren = Reynolds(body_radius,r,delv,z,sigma);
	knud = Knudsen(r,sigma,z,body_radius);
	rhogas = rho_midpl(r,sigma) * M_SOL / AU2CM / AU2CM / AU2CM;		/*	density in g/cm/cm/cm	*/

	Eps_Stokes_coeff(machn,ren,&ceps,&cstokes);

	cdrag = (9.0 * knud * knud * ceps + cstokes) / ((3.0 * knud + 1.0) * (3.0 * knud + 1.0));

	*adrag = -3.0 * rhogas * cdrag * delv / (8.0 * PDENSITY * body_radius / AU2CM);		/*	1/day	*/
	
}


/*	derivates for 4th order Runge-Kutta: x' = vx, y' = vy, z' = vz; vx' = ax, vy' = ay, vz' = az	*/
void eqrhs(double r, double sigma, double *pos_vel, double body_radius, double mass, double *del_vel, double *der) {

	double adrag_i,r3;
	get_drag_acc(r,sigma,pos_vel[2],body_radius,del_vel[0],del_vel[1],del_vel[2],&adrag_i);
	r3 = r * r * r;

	der[0] = pos_vel[3];
	der[1] = pos_vel[4];
	der[2] = pos_vel[5];
	der[3] = -G2 * (CBODY + mass) * pos_vel[0] / r3 + adrag_i * del_vel[0];			/*	+ adrag_i * dvx comes from "braking force" of gas	*/
	der[4] = -G2 * (CBODY + mass) * pos_vel[1] / r3 + adrag_i * del_vel[1];
	der[5] = -G2 * (CBODY + mass) * pos_vel[2] / r3 + adrag_i * del_vel[2];	

}


void Kunge_Kutta(double r, double sigma, double step, double x, double y, double z, double body_radius, double mass, double dvx, double dvy, double dvz, double vx, double vy, double vz, double *x_new, double *y_new, double *z_new, double *vx_new, double *vy_new, double *vz_new) {

	int n;
	double ytemp[6], dy1[6], dy2[6], dy3[6], dy4[6], pos_vel[6], del_vel[3], ynew[6];

	pos_vel[0] = x;
	pos_vel[1] = y;
	pos_vel[2] = z;
	pos_vel[3] = vx;
	pos_vel[4] = vy;
	pos_vel[5] = vz;
	del_vel[0] = dvx;
	del_vel[1] = dvy;
	del_vel[2] = dvz;

	eqrhs(r,sigma,pos_vel,body_radius,mass,del_vel,dy1);
	for(n = 0; n < 6; n++) ytemp[n] = pos_vel[n] + 0.5 * step * dy1[n];

	eqrhs(r,sigma,ytemp,body_radius,mass,del_vel,dy2);
	for(n = 0; n < 6; n++) ytemp[n] = pos_vel[n] + 0.5 * step * dy2[n];

	eqrhs(r,sigma,ytemp,body_radius,mass,del_vel,dy3);
	for(n = 0; n < 6; n++) ytemp[n] = pos_vel[n] + step * dy3[n];

	eqrhs(r,sigma,ytemp,body_radius,mass,del_vel,dy4);
	for(n = 0; n < 6; n++) ynew[n] = pos_vel[n] + step * (dy1[n] + 2.0 * dy2[n] * 2.0 * dy3[n] * dy4[n]) / 6.0;

	*x_new = ynew[0];
	*y_new = ynew[1];
	*z_new = ynew[2];
	*vx_new = ynew[3];
	*vy_new = ynew[4];
	*vz_new = ynew[5];

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

	double t, deltat, t_integration, DD, RMIN, RMAX, DD_cgs, RMIN_cgs, RMAX_cgs,rho,adrag_i,dvx_i,dvy_i,dvz_i;
	int i,num,linesout,count,partnum;
	linesout = 0;
	num = 0;
	count = 0;
	partnum = 0;

	t = 0.0;
	deltat = 0.1;
	t_integration = 1000.0*365.25;

	num = Sorokszama(linesout);						/*	A feluletisuruseg profilt tartalmazo file sorainak szama elmentve egy integerbe	*/
	linesout = 0;
	partnum = Reszecskekszama(linesout); 					/*	A reszecskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	
	double body_radius[partnum],x[partnum],y[partnum],z[partnum],vx[partnum],vy[partnum],vz[partnum],mass[partnum];	/*	A reszecskek adatait tarolo vektorok: body_radius a reszecske sugara, x,y,z a reszecsek hely-, vx,vy,vz a reszecskek sebessegvektorainak koordinatai	*/
	double sigmavec[num+2],rvec[num+2],dpressvec[num+2];			/*	Sigmavec: feluletisuruseget, rvec: az ehhez tartalmazo x pontokat, dpressvec: a nyomas derivaltat tartalmazo vektorok	*/
	double adrag, dvx, dvy, dvz;						/*	A reszecskekre hato surlodasbol szarmazo gyorsulas, es a relativ sebessegkomponensek	*/	
	double sigma_interp,dpress_interp,sigma_interp_vec[partnum+2],dpress_interp_vec[partnum+2];
	double v_kep_vec[partnum],rvec_part[partnum];
	double ux, uy, uz, sq_v_phi;
	double x_new[partnum],y_new[partnum],z_new[partnum],vx_new[partnum],vy_new[partnum],vz_new[partnum];
	double rmid, hist_i;
	int index;
	int hist[num];

	char hist_name[64], filname[64];

	load_r_sigma_dp(num,rvec,sigmavec,dpressvec,&DD,&RMIN,&RMAX);		/*	A feluletisuruseg es a nyomasderivalt profil betoltese es eltarolasa	*/
	load_rad_xyz_vxvyvz(partnum,body_radius,x,y,z,vx,vy,vz);		/*	A porreszecskek adatainak betoltese es eltarolasa	*/

	out = fopen("out.dat","w");

	for(i = 1; i < partnum; i++) {
		rvec_part[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
		interpol(RMIN,DD,rvec_part[i],sigmavec,dpressvec,rvec,&sigma_interp,&dpress_interp);
		sigma_interp_vec[i] = sigma_interp;					/*	M_sun / AU / AU --> g / cm / cm		*/
		dpress_interp_vec[i] = dpress_interp * DDENSASTR2CGS;			/*	M_Sun / AU / AU / AU / s / s --> g / cm / cm / cm / s / s	*/
		mass[i] = PDENSITY * 4.0 * M_PI / 3.0 * body_radius[i] * body_radius[i] * body_radius[i] / M_SOL;	/*	mass od the particle in solar mass	*/
		v_kep_vec[i] = v_kepler(rvec_part[i]);
		v_gas(rvec_part[i],sigma_interp_vec[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],dpress_interp_vec[i],&ux,&uy,&uz,&dvx,&dvy,&dvz);
		sq_v_phi = sqrt(ux * ux + uy * uy);
		get_drag_acc(rvec_part[i],sigma_interp_vec[i],z[i],body_radius[i],dvx,dvy,dvz,&adrag);
//		printf("%lg\n",adrag);
		fprintf(out,"%lg %lg %lg %lg %lg\n",rvec_part[i],v_kep_vec[i],sq_v_phi,sigma_interp,dpress_interp);
	}

	fclose(out);



	do {

		for (i = 0; i < num; i++) hist[i] = 0;

		for (i = 0; i < partnum; i++) {
			rvec_part[i] = sqrt(x[i] * x[i] + y[i] * y[i]);
						
			if(rvec_part[i] > RMIN)	{
				
				Kunge_Kutta(rvec_part[i],sigma_interp_vec[i],deltat,x[i],y[i],z[i],body_radius[i],mass[i],dvx,dvy,dvz,vx[i],vy[i],vz[i],&x_new[i],&y_new[i],&z_new[i],&vx_new[i],&vy_new[i],&vz_new[i]);
				x[i] = x_new[i];
				y[i] = y_new[i];
				z[i] = z_new[i];
				vx[i] = vx_new[i];
				vy[i] = vy_new[i];
				vz[i] = vz_new[i];
//				printf("%i time: %lg %lg %lg %lg\n",i,t,x[i],y[i],rvec_part[i]);
				histogram(num,rvec_part[i],hist,RMAX);	

			} else {

				x[i] = 0.0;
				y[i] = 0.0;
				z[i] = 0.0;
				vx[i] = 0.0;
				vy[i] = 0.0;
				vz[i] = 0.0;

			}		

		}
		
	

		if(fmod(t, (int)(t_integration/WO)) < deltat || t == 0) {		/*	a futasi ido alatt WO-szor irja ki, hogy hol tart	*/
			snprintf(hist_name,64,"hist_%i.dat", (int)(t));
			snprintf(filname,64,"out_%i.dat",(int)(t));
			histfil = fopen(hist_name,"w");
			out = fopen(filname,"w"); 
	
			printf("%lg\n",t);

  			for(i = 0; i < num; i++) {							
				fprintf(histfil,"%lg %i\n",rvec[i], hist[i]);		/*	file for histogram to find out the flux of particles	*/
			}

			fclose(histfil);

			for(i = 0; i < partnum; i++) {
				fprintf(out,"%lg %lg %lg %lg %lg\n",rvec_part[i],v_kep_vec[i],sq_v_phi,sigma_interp,dpress_interp);
			}
			fclose(out);
		}

		t = t + deltat;

	} while (t < t_integration);	

	return 0;

}
