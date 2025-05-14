#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TWOPI        	2.0 * M_PI
//#define G            	0.01720209895  //Gaussian grav. constant
#define G		1.0	
#define G2		G * G
#define STAR		1.0
#define HASP		0.05
//#define SIGMAP_EXP  	-0.5
//#define SIGMA0      	2.37838e-5 			/*	M_SUN / AU / AU				*/
//#define alpha_visc 	0.01				/*	alpha parameter				*/

#define SUN2GR	    	1.989e33
#define SDCONV      	1.12521e-7
#define AU2CM		1.496e13
//#define PDENSITY    	4.0				/*	a reszecskek atlag surusege		*/
#define CMPSECTOAUPYRP2PI	3.35725e-07	//conversion between the velocities

//#define PDENSITYDIMLESS	PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM 


#define NGRID       	2000
/*#define RMIN		1.0
#define RMAX		100.0
#define DGRID 		(RMAX - RMIN) / (NGRID - 1)
*/
#define STEPWRITE	10000				/*	hanyadik lepesenkent irassa ki az adatokat	*/

#define filenev		"init_data.dat"

double RMIN, RMAX;
int PARTICLE_NUMBER;
//double DD = DGRID;	/*	valamiert ha az interpolalasnal csak siman DGRID-et hasznalok, akkor nem jol szamol!	*/
double SIGMAP_EXP, SIGMA0;

double r_dze_i; // = 5.0;					
double r_dze_o; // = 20.0;						/*	az atmenet 2*H	*/
double Dr_dze_i; // = 2. * scale_height(r_dze_i);
double Dr_dze_o; // = 2. * scale_height(r_dze_o);
double a_mod; // = 0.01;
double PDENSITY;
double PDENSITYDIMLESS;// = PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;
double alpha_visc;

FILE *fin1, *fin2, *fout, *fmo, *fil, *fout2, *histfil, *massfil, *fout3; 

/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int reszecskek_szama(int numout){

	char c;
	fin1 = fopen(filenev,"r+");		
	numout = 0;

/*	A porreszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal				*/

		while((c = fgetc(fin1)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 							*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert										*/
	fclose(fin1);	

	return numout;

}


/*	A korong parametereinek beolvasasa	*/
void disk_param_be(double *sigma0, double *sdexp, double *Rmin, double *Rmax, double *r_dzei, double *r_dzeo, double *dr_dzei, double *dr_dzeo, double *alph_mod, double *rho_p, double *rho_p_dimless, double *alphav) {

	double dummy, rmin, rmax, drdzei, drdzeo, rdzei, rdzeo, amod, rhop, alpha;
	int dummy2;

	double sig0, exp;

	fin2 = fopen("disk_param.dat","r");

           	if(fscanf(fin2,"%lg  %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&rmin,&rmax,&dummy2,&exp,&sig0,&dummy,&rdzei,&rdzeo,&drdzei,&drdzeo,&amod,&rhop,&alpha) == 13) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 5, mert 5 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  sigma0: %lg, sdexp: %lg\n\n\n", sig0,exp);  
			*sigma0 = sig0;
			*sdexp = exp; 
			*Rmin = rmin;
			*Rmax = rmax;       
			*r_dzei = rdzei;
			*r_dzeo = rdzeo;
			*dr_dzei = drdzei;
			*dr_dzeo = drdzeo;
			*alph_mod = amod;
			*rho_p = rhop;
			*rho_p_dimless = rhop / SUN2GR * AU2CM * AU2CM * AU2CM;
			*alphav = alpha;
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
//			getchar();
			exit(EXIT_FAILURE);
   	        }
	fclose(fin2);
	
}


/*	A porreszecskek adatainak beolvasasa	*/
void por_be(double radius[][2], double *mass) {

	int i, dummy;
	double dummy2;
	double distance, particle_radius, reprmass;
	
   	fin1 = fopen(filenev,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%d %lg %lg %lg ",&dummy,&distance,&reprmass,&particle_radius) == 4) {//,&dummy2,&dummy2,&dummy2) == 7) {	
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 4, mert 4 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/

            		radius[i][0] = distance;
	   		radius[i][1] = particle_radius / AU2CM;		/* a részecske mérete AU-ban	*/
//			radius[i][1] = 1./AU2CM;
			mass[i] = reprmass;

	    	} else {

/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/					
			printf("\n\n*******************     ERROR!     *********************\n\n  Nem sikerult a %i-ik sort beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n",dummy);
//			getchar();
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(fin1);
	
	printf("\n\n *******   A file beolvasasa sikerult!   ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");	
//	getchar();

}


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma, double PDENSITYDIMLESS) {		/*	in the Epstein drag regime	*/

	return PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);

} 


/*	local scale height	*/
double scale_height(double r) {

	return r * HASP;

}


/*	lokális kepleri sebesség	*/
double v_kep(double r) {

	return sqrt(G2 * STAR / r);
	
}


/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r) {

	return sqrt(G2 * STAR / r / r / r);			/*	v_kepler in AU / (yr/2pi)	*/

}


/*	local sound speed		*/
double c_sound(double r) {

	return kep_freq(r) * scale_height(r);

}


/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r, double alpha_visc) {

	double alpha_r;
	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
   	return alpha_r ;

}


/*	Lokalis viszkozitas erteke	*/
double visc(double r, double alphvisc){
 
  	double nu, alpha_r;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);

/*	alpha_r: a redukcio merteke	*/
 	alpha_r = alpha_turb(r,alphvisc);

    	nu = alpha_r * cs * H;
  
  	return nu;
  
}


/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
double Coeff_1(double r, double alphvisc){					
  
  	double A;
  	A = 3.0 * visc(r,alphvisc);
  	return A;
  
}


/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/
double Coeff_2(double r, double alphvisc){							
  
  	double B;
  	B = 9.0 * visc(r,alphvisc) / (2.0 * r);
  	return B;
    
}


/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double DD, double RMIN) {

	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = RMIN + (i1-1) * DD;
	x2 = RMIN + (i2-1) * DD;
	x3 = RMIN + (i3-1) * DD;
 
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


void Perem(double *vec, double DD, double RMIN, double RMAX) {					/*	boundary condition for sigma		*/

	double a, b, c; 

	Parabola(vec, 1, 2, 3, &a, &b, &c, DD, RMIN);
	vec[0] =  a * (RMIN - DD) * (RMIN - DD) + b * (RMIN - DD) + c;

	Parabola(vec, NGRID - 2, NGRID - 1, NGRID, &a, &b, &c, DD, RMIN);
	vec[NGRID+1] = a * (RMAX + DD) * (RMAX + DD) + b * (RMAX + DD) + c;

	
/*  	sigmavec[0] = (sigmavec[1] / visc(RMIN,alphvisc)) * visc(RMIN - DD,alphvisc);
  	sigmavec[NGRID+1] = (sigmavec[NGRID] / visc(RMIN + (NGRID) * DD,alphvisc)) * visc(RMIN + (NGRID+1) * DD,alphvisc);
*/  
}


void Initial_Profile(double *sigmavec, double *r, double DD, double RMIN, double RMAX){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		sigmavec[i] = SIGMA0 * pow(r[i],SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}
  
  	Perem(sigmavec,DD,RMIN,RMAX);
}


double time_step(double *rvec, double DD, double alphvisc) {					/*	a diffúziós egyenlet megoldásához szükséges minimális időlépés	*/

	double A_max, stepping;
	int i;

	A_max = -10000.0;
	
	for(i = 0; i < NGRID; i++) {
		if(Coeff_1(rvec[i],alphvisc) > A_max) {
			A_max = Coeff_1(rvec[i],alphvisc);
		}
	}

	stepping = DD * DD / (2.0 * A_max);

	return stepping;

}


/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {

	return 1. / (2.0 * M_PI) * sigma / (HASP * r);

}


/* 	local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r) {

	return rho_mp(sigma,r) * c_sound(r) * c_sound(r);

}


void Initial_Press(double *pressvec, double *sigmavec, double *rvec, double DD, double RMIN, double RMAX){		/*	initial profile of pressure		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		pressvec[i] = press(sigmavec[i],rvec[i]);		
  	}
  
  	Perem(pressvec,DD,RMIN,RMAX);
}


/*	a nyomas derivaltja	*/
void dpress(double *dp, double *p, double DGRID) {

	int i;

	for(i = 1; i <= NGRID; i++) {	
		dp[i] = (p[i+1] - p[i-1]) / (2.0 * DGRID);
	}


}	


void Initial_dPress(double *dpressvec, double *pressvec, double DD, double RMIN, double RMAX){		/*	initial profile of pressure		*/

	dpress(dpressvec,pressvec,DD);
   	Perem(dpressvec,DD,RMIN,RMAX);

}


/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double r, double pdensdim, double *drdt) {
     
     	double P, H, dPdr, St, csound;
    
	St = Stokes_Number(pradius,sigma,pdensdim);
	H = scale_height(r);  
	P = press(sigma,r);
	dPdr = dp;
	csound = c_sound(r); 

     	*drdt = St / (1 + St * St) * H / P * dPdr * csound;
     
}


/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double DD, double *rvec, double pos, double *out) {

	double rmid, rindex, coef1, temp;
	int index; 

     	rmid = pos - RMIN;
	rmid = rmid / DD;     						/* 	the integer part of this gives at which index is the body	*/
	index = (int) floor(rmid+0.5);					/* 	ez az rmid egesz resze	(kerekites 0.5-tol)			*/
	rindex = rvec[index];       					/* 	the corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

 	coef1 = (invec[index + 1] - invec[index]) / DD; 		/*	ez az alabbi ket sor a linearis interpolacio - remelem, jo!!!	*/
	temp = invec[index] + coef1 * (pos - rindex);          		/*	a beerkezo dimenzionak megfelelo mertekegysegben		*/

	*out = temp;

}


//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret
/*double a_drift() {

	s_drift = f_drift * 2.0 / M_PI * Sigmad_cgs/rho_p * v_kep2 / c_s2 * fabs(1.0 / dlnPdlnr);

}
*/


// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret
double a_turb(double sigma, double r, double rho_p, double alphvisc) {

	double s_frag, f_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

  	f_frag  = 0.37;
  	u_frag = 1000.0; 	// cm/s 
	u_frag = u_frag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
  	u_frag2 = u_frag * u_frag;
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r); // / CMPSECTOAUPYRP2PI;
	c_s2 = c_s * c_s;

	s_frag = f_frag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(r,alphvisc)) * u_frag2 / c_s2;
	
	return s_frag;

}


// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret
double a_df(double sigma, double r, double p, double dp, double rho_p) {

	double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

  	u_frag = 1000.0; 	// cm/s 
	u_frag = u_frag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r);
	c_s2 = c_s * c_s;
	dlnPdlnr = r / p * dp;
	vkep = v_kep(r);

	s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);
	
	return s_df;

}


/*	minimum megkeresese harom elem kozul	*/
double find_min(double s1, double s2, double s3) {

	double min;
	
	if (s1 < s2) {
		if (s1 < s3) {
			min = s1;
		} else {
			min = s3;
		}
	} else {
		if (s2 < s3) {
			min = s2;	
		} else {
			min = s3;
		}
		
	}  

	return min;
}


//megiscsak vektorokkal kell majd csinalni? -- mit? :D -- ja! :D
/*	Runge-Kutta4 integrator	*/
void int_step(double time, double prad, double *pressvec, double *dpressvec, double *sigmavec, double *rvec, double DD, double step, double y, double pdensdim, double alphvisc, double *ynew, double *pradnew) {
	double dy1,dy2,dy3,dy4;
	double ytemp;
	double sigma, dpress; 
	double sturb, sdf, smin, pdens, p;

	interpol(sigmavec,DD,rvec,y,&sigma);
	interpol(dpressvec,DD,rvec,y,&dpress);

	if(time != 0.) {

		interpol(pressvec,DD,rvec,y,&p);
		pdens = pdensdim * SUN2GR / AU2CM / AU2CM / AU2CM;
		sturb = a_turb(sigma,y,pdens,alphvisc);
		sdf = a_df(sigma,y,p,dpress,pdens);
		smin = find_min(sturb,sdf,HUGE_VAL);
		prad = smin / AU2CM;
	}

	*pradnew = prad;

	eqrhs(prad, dpress, sigma, y, pdensdim, &dy1);
	
	ytemp = y + 0.5 * step * dy1;
	eqrhs(prad, dpress, sigma, ytemp, pdensdim, &dy2);
		
	ytemp = y + 0.5 * step * dy2;
	eqrhs(prad, dpress, sigma, ytemp, pdensdim, &dy3);
	
	ytemp = y + step * dy3;
	eqrhs(prad, dpress, sigma, ytemp, pdensdim, &dy4);

	*ynew = y + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;

}


/*	megkeresi egy tomb maximumat	*/
double find_max(double r[][2]) {

	int i;
	double maxim = -1.;

	for(i = 0; i < PARTICLE_NUMBER; i++) {

		if (r[i][0] > maxim) {
			maxim = r[i][0];
		}
	
	}


	return maxim;
}


/*	counting the number of zero points of the pressure gradient function (and also pressure maxima)	*/
int find_num_zero(double *rvec, double *dp) {

	int i,count;
	count = 0;

	for(i = 0; i < NGRID-1; i++) {
		if(((dp[i] * dp[i+1]) <= 0.)  && (dp[i] > dp[i+1])) {	/*	Osszeszorozza a ket ponton a nyomas derivaltjanak erteket, ahol a szorzat negativ, ott elojelvaltas tortenik --> negativbol pozitivba, vagy pozitivbol negativba valt --> nyomasi maximum. Maximum pedig ott talalhato, ahol a fuggveny pozitivbol negativba valt at (ezt keresi a masodik feltetel).	*/
			count++;
		} 

	}

	return count;
}


/*	solving a*x + b = y (here a = r, y = dp)	*/
double find_r_zero(double r1, double r2, double dp1, double dp2) {

	double a, b, r_zero;
	a = (dp2 - dp1) / (r2 - r1);
	b = dp1 - a * r1;
	r_zero = - b / a;

	return r_zero;

}


/*	this function counts where (which r) the pressure maximum is	*/
double find_zero(int i, double *rvec, double *dp) {

	double r;
	
	if(((dp[i] * dp[i+1]) <= 0.) && (dp[i] > dp[i+1])) {
		r = find_r_zero(rvec[i],rvec[i+1],dp[i],dp[i+1]);
	} else {
		r = 0.0;
	}

	return r;

}


/*	A nyomasi maximum korul 2H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(double DD, double *rvec, double rout, double *ind_oi, double *ind_oo) {

	int i;
	double rmid, rtemp;
	double H = HASP * rout;

	for(i = 1; i <= NGRID; i++) {
	  
		if(rvec[i] > (rout - H) - DD / 2.0 && rvec[i] < (rout - H) + DD / 2.0) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oi = rtemp;
		}
				
		if(rvec[i] > (rout + H) - DD / 2.0 && rvec[i] < (rout + H) + DD / 2.0) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oo = rtemp;
		}

	}

}


int main() {
   
   	int i, L, linesout, dim, j;
   	double y, y_out, t, t_integration, deltat, particle_radius, sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2];
	double max;
	char dens_name[64],size_name[64];
   
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	PARTICLE_NUMBER = reszecskek_szama(linesout);  		

/*	A korong parametereit beolvassa: a sigma profil kitevojet es a sigma0-t					*/
	disk_param_be(&SIGMA0, &SIGMAP_EXP, &RMIN, &RMAX, &r_dze_i, &r_dze_o, &Dr_dze_i, &Dr_dze_o, &a_mod, &PDENSITY, &PDENSITYDIMLESS, &alpha_visc);

	printf("sig0: %lg, sdex: %lg, rmin: %lg, rmax: %lg, r_dze_i: %lg, r_dze_o: %lg, drdzei: %lg, drdzeo: %lg, amod: %lg\n", SIGMA0, SIGMAP_EXP, RMIN, RMAX, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod);	
//	getchar();

	double DD = (RMAX - RMIN) / (NGRID - 1);

   	double radius[PARTICLE_NUMBER][2];
	double drdt[PARTICLE_NUMBER], massvec[PARTICLE_NUMBER];


/*	A reszecskek adatainak beolvasasa file-bol								*/
	por_be(radius,massvec);				


 	for(i = 0; i <= NGRID+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = RMIN + (i-1) * DD;
	}

/*	Kezdeti profilok betoltese a megadott vektorokba	*/
	Initial_Profile(sigmavec,rvec,DD,RMIN,RMAX);
	Initial_Press(pressvec,sigmavec,rvec,DD,RMIN,RMAX);
	Initial_dPress(dpressvec,pressvec,DD,RMIN,RMAX);
	
  	t = 0.0;
   	t_integration = 500000.0 * 2.0 * M_PI; 					/*	numerikus integralas idotartama		*/
	deltat = time_step(rvec,DD,alpha_visc)/5.;						/*	idolepes	*/
	printf("dt: %lg\n",deltat);	
	

   	L = 0;

   	fout = fopen("pormozgas.dat","w");
   	massfil = fopen("mass.dat","w");

	double partmassind[PARTICLE_NUMBER][3];

   	do {

		if(t == 0.0) {							/*	A t=0 idopillanatban elkesziti az alabbi file-okat:	*/

			double rmid;
			int rindex;
		
  			fmo = fopen("surface.dat","w");				/*	surface.dat, ez tartalmazza a korong kezdeti idopillanati adatait	*/

 				for(i = 1; i <= NGRID; i++) {
   					fprintf(fmo, "%lg   %lg   %lg   %lg  %lg\n", rvec[i], sigmavec[i],pressvec[i],dpressvec[i],c_sound(rvec[i])); 				}

			fclose(fmo);

			fout2 = fopen("timescale.dat","w");			
/*	timescale.dat, ez tartalmazza a reszecskek bearamlasanak idoskalajat (evben)				*/
/*	A timescale.dat-hoz a t_drift = r / v(r,peak) kepletet a kovetkezo modon oldja meg a program:		*/
/*	Kiszamolja, hogy a kovetkezo idopillanatban hol lenne minden reszecske, majd az uj koordinatabol (x(1))	*/
/*	a regit (x(0)) levonva, majd az egeszet elosztva dt-vel megkapjuk a bearamlas sebesseget (drdt)		*/
/*	ez utan kiszamolhato az eredeti t_drift keplet, ezt iratom ki a file-ba evben (ezert osztok le 2PI-vel	*/

			fout3 = fopen("size.0.dat","w");

			for (i = 0; i < PARTICLE_NUMBER; i++) {	

				fprintf(fout3,"%lg  %lg  %lg\n",t / (2.0 * M_PI), radius[i][0], radius[i][1]*AU2CM);
				fprintf(fout,"%lg %d %lg\n",t / (2.0 * M_PI),i,radius[i][0]);
				fflush(fout);

  	  			rmid = (radius[i][0] - RMIN) / DD;     						/* 	The integer part of this gives at which index is the body	*/
				rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
	
				partmassind[i][0] = (double)i;							/*	index of the particles						*/
				partmassind[i][1] = massvec[i];		/*	mass of the particles						*/
				partmassind[i][2] = rindex;							/*	initial distance of the particles				*/

				if (radius[i][0] > RMIN && radius[i][0] < RMAX) { 

					double temp;

					y = radius[i][0];
		     			particle_radius = radius[i][1]; 

		     			int_step(t,particle_radius,pressvec,dpressvec,sigmavec,rvec,DD,deltat,y,PDENSITYDIMLESS,alpha_visc,&y_out,&temp);

					drdt[i] = (fabs(y_out - y)/(deltat));
	
				} else {

					y_out = 0.0;
				}

				fprintf(fout2,"%lg %lg\n",radius[i][0], (radius[i][0] / drdt[i])/2.0/M_PI);

			}
	
			fclose(fout2);
			fclose(fout3);

			double ind_oi, ind_oo;			

			find_r_annulus(DD,rvec,r_dze_o,&ind_oi,&ind_oo);		/*	A nyomasi maximum korul 2H tavolsagban keres korgyurut, a fuggveny visszaadja a cellak indexet	*/


			int index_i;
			double mass = 0.;

			for (i = 0; i < NGRID; i++) {

				index_i = partmassind[i][2];

				double masstemp = 0.;
				if ((index_i >= (int)ind_oi) && (index_i <= (int)ind_oo)) {
					
					masstemp = partmassind[i][1];
					mass = mass + masstemp;

				}
			}


			fprintf(massfil,"%lg %lg\n",t/2./M_PI,mass);
			fflush(massfil);


		}

		max = find_max(radius);					/*	Megkeresi, hogy melyik a legtavolabbi reszecske a kozponti csillagtol	*/

		if(max >= RMIN) {							/*	Ha a legtavolabbi porreszke merete nagyobb RMIN-nel, akkor folytatodik az integralas, kulonben kilep	*/

  			for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
  				sigmavec[i] = (sigmavec[i] * visc(rvec[i],alpha_visc) + deltat * (Coeff_1(rvec[i],alpha_visc) * (sigmavec[i+1] * visc(rvec[i+1],alpha_visc) - 2.0 * sigmavec[i] * visc(rvec[i],alpha_visc) + sigmavec[i-1] * visc(rvec[i-1],alpha_visc)) / (DD * DD) + Coeff_2(rvec[i],alpha_visc) * (sigmavec[i+1] * visc(rvec[i+1],alpha_visc) - sigmavec[i-1] * visc(rvec[i-1],alpha_visc)) / (2.0 * DD))) / visc(rvec[i],alpha_visc);
				pressvec[i] = press(sigmavec[i],rvec[i]);
    
  			}

  			Perem(sigmavec,DD,RMIN,RMAX);				/*	loading boundary condition in each timestep		*/
			dpress(dpressvec,pressvec,DD);					/*	loading boundary condition in each timestep		*/
			Perem(pressvec,DD,RMIN,RMAX);						/*	loading boundary condition in each timestep		*/
			Perem(dpressvec,DD,RMIN,RMAX);					/*	loading boundary condition in each timestep		*/
	
			t = t + deltat;							/*	Idoleptetes	*/
			L++;								/*	Idolepes szamlalo	*/

			double prad_new;
			for (i = 0; i < PARTICLE_NUMBER; i++) {	
		
				if (radius[i][0] >= RMIN+DD && radius[i][0] <= RMAX-DD) { 	/*	Ha a reszecske benne RMIN es RMAX kozott van, akkor megcsinalja az integralast	*/

					y = radius[i][0];
			     		particle_radius = radius[i][1]; 
//					particle_radius = 5.0 / AU2CM;

			     		int_step(t,particle_radius,pressvec,dpressvec,sigmavec,rvec,DD,deltat,y,PDENSITYDIMLESS,alpha_visc,&y_out,&prad_new);/*	Integralas	*/
	
					radius[i][1] = prad_new;
			     		radius[i][0] = y_out;
	
				} else {						/*	Ha a reszecske nincs RMIN es RMAX kozott, akkor a program kinullazza a reszecske tavolsagat	*/

					radius[i][0] = 0.0;

				}
		     
			     	if (L%STEPWRITE==0) {						/*	Megadott idolepesenkent kiirja a reszecske tavolsagat es az idot a pormozgas.dat file-ba	*/					
			         	fprintf(fout,"%lg %d %lg\n",t / (2.0 * M_PI),i,radius[i][0]);
					fflush(fout);

			     	}
			

			}
		
			if (L%STEPWRITE==0) {							/*	Megadott idolepesenkent kiirja a korong parametereit egy file-ba. A file neve tartalmazza az idot evben	*/
				snprintf(dens_name,64,"dens_%lg.dat", t / (2.0 * M_PI));
				snprintf(size_name,64,"size.%lg.dat",t / (2.0 * M_PI));
				fout3 = fopen(size_name,"w"); 
				fil = fopen(dens_name,"w");

				dim = find_num_zero(rvec,dpressvec);

				double r_count[dim];
				double temp_new = 0.;
				double temp = 0.;
				double rout = r_dze_o;
				double ind_oi, ind_oo;

				j = 0;
	
				if(dim != 0) {
					for(i = 0; i < NGRID; i++) {

						temp_new = find_zero(i,rvec,dpressvec);
						if(temp != temp_new && i > 3 && temp_new != 0.0) {
							r_count[j] = temp_new;
							j++;
						}
						
						if(temp_new > 0.) {
							temp = temp_new;
							rout = temp;
						} 
					}

				}

				find_r_annulus(DD,rvec,rout,&ind_oi,&ind_oo);		/*	A nyomasi maximum korul 2H tavolsagban keres korgyurut, a fuggveny visszaadja a cellak indexet	*/

				int index_i;
				double mass = 0.;
				
				for(i = 0; i < PARTICLE_NUMBER; i++) {
	
					if (radius[i][0] > RMIN) fprintf(fout3,"%lg %lg %lg\n",t / (2.0 * M_PI), radius[i][0],radius[i][1] * AU2CM);

					index_i = partmassind[i][2];

					double masstemp = 0.;
					if ((index_i >= (int)ind_oi) && (index_i <= (int)ind_oo)) {

						masstemp = partmassind[(int)index_i][1];
//						mass_growthvec_r2[count] += mass;
						mass = mass + masstemp;

					}
					

				}

				fclose(fout3);
			
				double rmid;
				int rindex;

				for(i = 0; i < PARTICLE_NUMBER; i++) {
					rmid = (radius[i][0] - RMIN) / DD;     					/* 	The integer part of this gives at which index is the body	*/
					rindex = (int) floor(rmid);						/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
					partmassind[i][2] = rindex;

				}

				if (mass != 0.) {

					fprintf(massfil,"%lg %lg\n",t/2./M_PI,mass);
					fflush(massfil);
				}




		
		  		for(i = 1; i <= NGRID; i++) {							
					fprintf(fil,"%lg %lg %lg %lg %lg\n",rvec[i], sigmavec[i], pressvec[i], dpressvec[i],c_sound(rvec[i]));
				}
	
				fclose(fil);

				printf("t: %lg\n",t / (2.0 * M_PI));
		
			}


		} else {				/*	Ha a legmesszebbi reszecske tavolsaga mar nem nagyobb, vagy egyenlo, mint RMIN, akkor a program "figyelmezteto szoveg" mellett sikeresen kilep, nem fut "feleslegesen" tovabb.	*/
			printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",t/2./M_PI);
//			getchar();
			exit(EXIT_SUCCESS);
		}

   	} while (t < t_integration);
  

/*	Az idoleptetes leteltevel a program sikeresen kilep	*/
	printf("\n\nA program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n",max); 
//	getchar();


	return 0;

}
