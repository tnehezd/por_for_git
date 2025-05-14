#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NGRID 		1500
#define RMIN		1.
#define	RMAX		50.
#define	DD		(RMAX-RMIN)/((double)NGRID - 1)
#define SIGMA0		1.25e-5
#define	SIGMAP_EXP	-1.5
#define TWOPI        	2.0 * M_PI
#define G		1.0	
#define G2		G * G

/*	makro a korong geometriajara: aspect ratio	*/
#define HASP		0.05

/*	atvaltasi ertekek				*/
#define SUN2GR	    	1.989e33
#define SDCONV      	1.12521e-7
#define AU2CM		1.496e13
#define CMPSECTOAUPYRP2PI	3.35725e-07	//conversion between the velocities

#define STAR		1.0
#define FLIND		0.0
#define HASP		0.05
#define	PDENSITYDIMLESS	1.6/ SUN2GR * AU2CM * AU2CM * AU2CM

#define EPS		0.01			// por/gaz arany

double	pradius	= 20./ AU2CM;
double	fm = 0.97;

/*	r vektor (gridcellák) inicializálása	*/
void load_R(double *rvec) {
	
	int i;
 	for(i = 0; i <= NGRID+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = RMIN + (i-1) * DD;
	}

}


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma) {		/*	in the Epstein drag regime	*/

	return PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);

} 


/*	lokális kepleri sebesség	*/
double v_kep(double r) {

	return sqrt(G2 * STAR / r);
	
}


/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r) {

	return sqrt(G2 * STAR / r / r / r);			/*	v_kepler in AU / (yr/2pi)	*/

}


/*	local scale height	*/
double scale_height(double r) {

	return pow(r,1.+FLIND) * HASP;

}


/*	local sound speed		*/
double c_sound(double r) {

	return kep_freq(r) * scale_height(r);

}


/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r) {

	double alpha_r;
	double a_mod = 0.01;
	double r_dze_i = 2.7;
	double r_dze_o = 35.;
	double Dr_dze_i = 1.5;
	double Dr_dze_o = 1.5;
	double alpha_visc = 0.01;
	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
//	alpha_r = 1.;
   	return alpha_r * alpha_visc;

}




/*	Lokalis viszkozitas erteke	*/
double visc(double r){
 
  	double nu, alpha_r;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);

/*	alpha_r: a redukcio merteke	*/
 	alpha_r = alpha_turb(r);

    	nu = alpha_r * cs * H;
  
  	return nu;
  
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


/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {

	return 1. / (2.0 * M_PI) * sigma / scale_height(r);

}



/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd) {

	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = RMIN + (i1-1) * dd;
	x2 = RMIN + (i2-1) * dd;
	x3 = RMIN + (i3-1) * dd;
 
	y1 = vec[i1];
	y2 = vec[i2];
	y3 = vec[i3];

	av = ((y1 - y3) / (x1 - x3) - (y1 - y2) / (x1 - x2)) / (x3 - x2);
	bv = (y1 - y2) / (x1 - x2) - av * (x1 + x2);
	cv = y1 - av * x1 * x1 - bv * x1;

	*a = av;
	*b = bv;
	*c = cv;

}


/*	A peremen parabolat illeszt	*/
void Perem(double *vec) {					/*	boundary condition for sigma, p, dp...	*/

	double a, b, c; 

	Parabola(vec, 1, 2, 3, &a, &b, &c, DD);
	vec[0] =  a * (RMIN - DD) * (RMIN - DD) + b * (RMIN - DD) + c;

	Parabola(vec, NGRID - 2, NGRID - 1, NGRID, &a, &b, &c, DD);
	vec[NGRID+1] = a * (RMAX + DD) * (RMAX + DD) + b * (RMAX + DD) + c;

}



/* 	local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r) {

	return rho_mp(sigma,r) * c_sound(r) * c_sound(r);

}


void Initial_Press(double *pressvec, double *sigmavec, double *rvec){		/*	initial profile of pressure		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		pressvec[i] = press(sigmavec[i],rvec[i]);
  	}
  
  	Perem(pressvec);

}



/*	a nyomas derivaltja	*/
void dpress(double *dp, double *p) {

	int i;
	double ptemp, pvec[NGRID+1];

	for(i = 1; i <= NGRID; i++) {
		ptemp = (p[i+1] - p[i-1]) / (2.0 * DD);
		pvec[i] = ptemp;
	}
	
	for(i = 1; i <= NGRID; i++) {
		dp[i] = pvec[i];

	}

}


void Initial_dPress(double *dpressvec, double *pressvec){		/*	initial profile of pressure		*/

	dpress(dpressvec,pressvec);
   	Perem(dpressvec);

}


/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
double Coeff_3(double sigma, double r){

	return -1.0 * (3.0 / (sigma * sqrt(r)));

}



/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(double *sigmavec, double *rvec, double *ug) {

	double tempug, ugvec[NGRID+2],ugvectemp[NGRID+1];
	int i;

	for(i = 0; i <= NGRID+1; i++) {
		ugvec[i] = sigmavec[i] * visc(rvec[i]) * sqrt(rvec[i]);
	}

	for(i = 1; i <= NGRID; i++) {
		tempug = (ugvec[i+1] - ugvec[i-1]) / (2.0 * DD);
		ugvectemp[i] = Coeff_3(sigmavec[i],rvec[i]) * tempug;
	}
	
	for(i = 1; i <= NGRID; i++) {
		ug[i] = ugvectemp[i];
	}

}


/*	ug vektor feltoltese az u_gas ertekevel	*/
void Initial_Ugas(double *sigmavec, double *rvec, double *ug){		/*	initial profile of pressure		*/

	u_gas(sigmavec,rvec,ug);
  	Perem(ug);
}



/*	a sigmara kezdeti profil betoltese	*/
void Initial_Profile(double *sigmavec, double *r){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		sigmavec[i] = SIGMA0 * pow(r[i],SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}
  
  	Perem(sigmavec);

}


/*	a sigmara kezdeti profil betoltese	*/
void Initial_DustProfile(double *sigmavec, double *sigmadvec){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 0; i <= NGRID+1; i++) {
    		sigmadvec[i] = sigmavec[i] * EPS;		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}
  

}


double time_step(double *rvec) {					/*	a diffúziós egyenlet megoldásához szükséges minimális időlépés	*/

	double A_max, stepping;
	int i;

	A_max = -10000.0;
	
	for(i = 0; i < NGRID; i++) {
		if(Coeff_1(rvec[i]) > A_max) {
			A_max = Coeff_1(rvec[i]);
		}
	}

	stepping = DD * DD / (2.0 * A_max);

	return stepping;

}


/*	Fuggveny a sigma, p, dp kiszamolasara	*/
void Get_Sigma_P_dP(double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double deltat) {

	double u, u_bi, u_fi, sigma_temp[NGRID+2], uvec[NGRID+2];
	int i;

	sigma_temp[0] = sigmavec[0];
	sigma_temp[NGRID+1] = sigmavec[NGRID+1];
	uvec[0] = sigmavec[0] * visc(rvec[0]);
	uvec[NGRID+1] = sigmavec[NGRID+1] * visc(rvec[NGRID+1]);

	for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
		u = sigmavec[i] * visc(rvec[i]);			/*	az elozo lepesbol szarmazo sigma*nu elmentese	*/
		u_bi = sigmavec[i-1] * visc(rvec[i-1]);
		u_fi = sigmavec[i+1] * visc(rvec[i+1]);
		uvec[i] = u;
/*	az adott idolepesben az uj sigma*nu kiszamolasa a diffuzios egyenletbol	*/
		double temp = Coeff_1(rvec[i]) * (u_fi - 2.0 * u + u_bi) / (DD * DD) + Coeff_2(rvec[i]) * (u_fi - u_bi) / (2.0 * DD);
		sigma_temp[i] = uvec[i] + deltat * temp;		/*	regi + dt*uj	*/
	}

	for(i = 1; i <= NGRID; i++) {
		sigmavec[i] = sigma_temp[i]/visc(rvec[i]);		/*	nekunk sigma kell, nem sigma*nu		*/
		pressvec[i] = press(sigmavec[i],rvec[i]);		/*	nyomas kiszamolasa az uj sigmaval	*/
	}

	Perem(sigmavec);						/*	loading boundary condition in each timestep		*/
	dpress(dpressvec,pressvec);					/*	loading boundary condition in each timestep		*/
	Perem(pressvec);						/*	loading boundary condition in each timestep		*/
	Perem(dpressvec);						/*	loading boundary condition in each timestep		*/
	
}



void Get_udust(double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugas, double *udust) {

	int i;
	double dlnpdlnr;
	double St, cs;

	for(i = 1; i<NGRID+1; i++) {
		St = Stokes_Number(pradius,sigmavec[i]);
		cs = c_sound(rvec[i]);
		dlnpdlnr = rvec[i] / pressvec[i] * dpressvec[i];
		udust[i] = ugas[i] / (1. + St * St) + 2.0 / (St + (1. / St)) * cs * cs / (2.0 * v_kep(rvec[i])) * dlnpdlnr;
	}	

	Perem(udust);

}


void Get_Sigma_dust(double *rvec, double *sigmadvec, double *udust, double deltat) {

	int i;
	double u, sigmad_temp[NGRID+2], u_bi, u_fi, coeff, tempvec[NGRID+2];


	sigmad_temp[0] = sigmadvec[0];
	sigmad_temp[NGRID+1] = sigmadvec[NGRID+1];
	tempvec[0] = sigmadvec[0] * udust[0];
	tempvec[NGRID+1] = sigmadvec[NGRID+1] * udust[NGRID+1];	

	printf("r0: %lg sd0: %lg   sdveg:%lg\n",rvec[0],sigmadvec[0],sigmadvec[1]);
	getchar();

	for(i = 1; i < NGRID+1; i++) {

		u = sigmadvec[i] * udust[i];
		u_bi = sigmadvec[i-1] * udust[i-1];		
		u_fi = sigmadvec[i+1] * udust[i+1];
//		coeff = u / rvec[i];
		tempvec[i] = u;

		double temp = - 1.0 * (u + rvec[i] * (u_fi - u_bi) / (2.0 * DD));
		sigmad_temp[i] = tempvec[i] + deltat * temp;		/*	regi + dt*uj	*/

	}

	for(i = 1; i < NGRID+1; i++) {
		sigmadvec[i] = sigmad_temp[i] / udust[i];
		printf("sdtemp: %lg  s: %lg   u: %lg\n",sigmad_temp[i],sigmadvec[i], udust[i]);					
	}

//	Perem(sigmadvec);

	sigmadvec[0] = sigmadvec[1] + (sigmadvec[2] - sigmadvec[1])/ (rvec[2] - rvec[1]);  

/*	double a, b, y;
	a = (sigmadvec[2] - sigmadvec[1])/ (rvec[2] - rvec[1]);
	b = sigmadvec[1] - a * rvec[1];
	
	sigmadvec[0] = a * rvec[0] + b;
	sigmadvec[NGRID+1] = sigmadvec[NGRID]/(rvec[NGRID+1]);
	
	printf("r0: %lg sd0: %lg   sdveg:%lg\n",rvec[0],sigmadvec[0],sigmadvec[NGRID+1]);

	getchar();
		
*/

}


/*	Itt vegzi el az integralast, ha szukseg van ra	*/
//void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec) {
void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec, double *sigmadvec) {
	double TMAX = 150000.;
	double L = 0;
	double deltat = time_step(rvec);
	double WO = 150.;
//	double deltat = 0.5;
	double t_integration = TMAX * 2.0 * M_PI;
	WO = TMAX / WO;
	double t = 0;

	double udustvec_mic[NGRID+2], udustvec_large[NGRID+2], udustvec[NGRID+2];

	int i;

	char dens_name[1024];

   	do {

		double time = t / 2.0 / M_PI;


		if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){

			snprintf(dens_name,1024,"surface_%i.dat",(int)time);
			FILE *mikell;
			mikell=fopen(dens_name,"w");
			for(i = 1; i < NGRID+1; i++) fprintf(mikell,"%lg   %lg  %lg  %lg  %lg  %lg  %lg\n",rvec[i],sigmavec[i],sigmadvec[i],pressvec[i],dpressvec[i],ugvec[i], udustvec[i]);
			fclose(mikell);
			L = L+(double)(TMAX/WO);
		}
		
		


		Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
		u_gas(sigmavec,rvec,ugvec);
//		Get_udust(rvec,sigmavec,pressvec,dpressvec,ugvec,udustvec_mic);
		Get_udust(rvec,sigmavec,pressvec,dpressvec,ugvec,udustvec_large);

//		for(i = 1; i < NGRID+2; i++) udustvec[i] = (1. - fm) * udusvec_mic[i] + fm * udustvec_large[i];
		for(i = 1; i < NGRID+2; i++) {
			udustvec[i] = fm * udustvec_large[i];
//			printf("r: %lg  u: %lg   ul: %lg\n",rvec[i], udustvec_large[i], udustvec[i]);
		}
		Perem(udustvec);

		Get_Sigma_dust(rvec,sigmadvec,udustvec,deltat);


		t = t + deltat;	

   	} while (t <= t_integration);

}


int main(int argc, const char **argv) {


   	double sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], ugvec[NGRID+2], sigmadvec[NGRID+2];
	load_R(rvec);
	Initial_Profile(sigmavec,rvec);
	Initial_DustProfile(sigmavec,sigmadvec);
	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	Initial_Ugas(sigmavec,rvec,ugvec);

	char nev[1024];

	sprintf(nev,"izemiedummy.dat");

	FILE *izeke;

/*	izeke=fopen(nev,"w");

	int i;
	for(i=0; i<NGRID+2; i++) fprintf(izeke,"%lg 	%lg  \n",rvec[i],sigmavec[i]);

	fclose(izeke);
*/

//	tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec);
	tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec,sigmadvec);


	return 0;

}
