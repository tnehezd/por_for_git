#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TWOPI        	2.0 * M_PI
#define G		1.0	
#define G2		G * G
#define HASP		0.05

#define SUN2GR	    	1.989e33
#define SDCONV      	1.12521e-7
#define AU2CM		1.496e13
#define CMPSECTOAUPYRP2PI	3.35725e-07	//conversion between the velocities

#define NGRID       	2000

#define filenev1	"init_data.dat"
#define filenev2	"disk_param.dat"
#define filenev3	"time.dat"

double RMIN, RMAX;
int PARTICLE_NUMBER;
double SIGMAP_EXP, SIGMA0;

double r_dze_i; // = 5.0;					
double r_dze_o; // = 20.0;						/*	az atmenet 2*H	*/
double Dr_dze_i; // = 2. * scale_height(r_dze_i);
double Dr_dze_o; // = 2. * scale_height(r_dze_o);
double a_mod; // = 0.01;
double PDENSITY;
double PDENSITYDIMLESS;// = PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;
double alpha_visc;
double STAR;
double TMAX;
double WO;
double TCURR;
double FLIND;
double DD;

FILE *fin1, *fin2, *fout, *fmo, *fil, *fout2, *histfil, *massfil, *fout3, *jelfut; 

/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int reszecskek_szama(int numout){

	char c;
	fin1 = fopen(filenev1,"r+");		
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
void timePar(double *tMax, double *step, double *current) {

	double tmax,stepping,curr;

	fin2 = fopen(filenev3,"r");

           	if(fscanf(fin2,"%lg  %lg %lg",&tmax,&stepping,&curr) == 3) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 14, mert 14 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  tmax: %lg, a program %lg evenkent irja ki a file-okat\n\n\n",tmax,stepping);  
			*tMax = tmax;
			stepping = tmax/stepping;
			*step = stepping; 
			*current = curr;
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			exit(EXIT_FAILURE);
   	        }
	fclose(fin2);
	
}


/*	A korong parametereinek beolvasasa	*/
void disk_param_be(double *sigma0, double *sdexp, double *Rmin, double *Rmax, double *r_dzei, double *r_dzeo, double *dr_dzei, double *dr_dzeo, double *alph_mod, double *rho_p, double *rho_p_dimless, double *alphav, double *mStar, double *gamma) {

	double dummy, rmin, rmax, drdzei, drdzeo, rdzei, rdzeo, amod, rhop, alpha, mstar, flind;
	int dummy2;

	double sig0, exp;

	fin2 = fopen(filenev2,"r");

           	if(fscanf(fin2,"%lg  %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&rmin,&rmax,&dummy2,&exp,&sig0,&dummy,&rdzei,&rdzeo,&drdzei,&drdzeo,&amod,&rhop,&alpha,&mstar,&flind) == 15) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 14, mert 14 oszlopot szeretnénk beolvasni.	*/
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
			*mStar = mstar;
			*gamma = flind;			

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
	double distance, particle_radius, reprmass;
	
   	fin1 = fopen(filenev1,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%d %lg %lg %lg ",&dummy,&distance,&reprmass,&particle_radius) == 4) {	
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 4, mert 4 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/

            		radius[i][0] = distance;
	   		radius[i][1] = particle_radius / AU2CM;		/* a részecske mérete AU-ban	*/
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


/*	local scale height	*/
double scale_height(double r) {

	return pow(r,1.+FLIND) * HASP;

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
double alpha_turb(double r) {

	double alpha_r;
	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
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


/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c) {

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


void Perem(double *vec) {					/*	boundary condition for sigma		*/

	double a, b, c; 

	Parabola(vec, 1, 2, 3, &a, &b, &c);
	vec[0] =  a * (RMIN - DD) * (RMIN - DD) + b * (RMIN - DD) + c;

	Parabola(vec, NGRID - 2, NGRID - 1, NGRID, &a, &b, &c);
	vec[NGRID+1] = a * (RMAX + DD) * (RMAX + DD) + b * (RMAX + DD) + c;

}


void Initial_Profile(double *sigmavec, double *r){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		sigmavec[i] = SIGMA0 * pow(r[i],SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}
  
  	Perem(sigmavec);
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


/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {

	return 1. / (2.0 * M_PI) * sigma / scale_height(r);

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


double Coeff_3(double sigma, double r){

	return -1.0 * (3.0 / (sigma * sqrt(r)));

}

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

void Initial_Ugas(double *sigmavec, double *rvec, double *ug){		/*	initial profile of pressure		*/

	u_gas(sigmavec,rvec,ug);
  	Perem(ug);
}



/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt) {
     
     	double P, H, dPdr, St, csound;
    
	St = Stokes_Number(pradius,sigma);
	H = scale_height(r);  
	P = press(sigma,r);
	dPdr = dp;
	csound = c_sound(r); 

     	*drdt = ug / (1. + St * St) + St / (1. + St * St) * H / P * dPdr * csound;
     
}


/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out) {

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
double a_turb(double sigma, double r, double rho_p) {

	double s_frag, f_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

  	f_frag  = 0.37;
  	u_frag = 1000.0; 	// cm/s 
	u_frag = u_frag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
  	u_frag2 = u_frag * u_frag;
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r); // / CMPSECTOAUPYRP2PI;
	c_s2 = c_s * c_s;

	s_frag = f_frag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(r)) * u_frag2 / c_s2;
	
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
void int_step(double time, double prad, double *pressvec, double *dpressvec, double *sigmavec, double *rvec, double *ugvec, double step, double y, double *ynew, double *pradnew) {
	double dy1,dy2,dy3,dy4;
	double ytemp;
	double sigma, dpress, ugas; 
	double sturb, sdf, smin, pdens, p;

	interpol(sigmavec,rvec,y,&sigma);
	interpol(dpressvec,rvec,y,&dpress);
	interpol(ugvec,rvec,y,&ugas);

	if(time != 0.) {

		interpol(pressvec,rvec,y,&p);
		pdens = PDENSITY; //DIMLESS * SUN2GR / AU2CM / AU2CM / AU2CM;
		sturb = a_turb(sigma,y,pdens);
		sdf = a_df(sigma,y,p,dpress,pdens);
		smin = find_min(sturb,sdf,HUGE_VAL);
		prad = smin / AU2CM;
	}

	*pradnew = prad;

	eqrhs(prad, dpress, sigma, ugas, y, &dy1);
	
	ytemp = y + 0.5 * step * dy1;
	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy2);
		
	ytemp = y + 0.5 * step * dy2;
	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy3);
	
	ytemp = y + step * dy3;
	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy4);

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
void find_r_annulus(double *rvec, double rout, double *ind_oi, double *ind_oo) {

	int i;
	double rmid, rtemp;
	double roimH = (rout - scale_height(rout)) - DD / 2.0;
	double roipH = (rout - scale_height(rout)) + DD / 2.0;
	double roomH = (rout + scale_height(rout)) - DD / 2.0;
	double roopH = (rout + scale_height(rout)) + DD / 2.0;

	for(i = 1; i <= NGRID; i++) {

		double H = scale_height(rvec[i]);

		if(rvec[i] > roimH && rvec[i] < roipH) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oi = rtemp;
		}
				
		if(rvec[i] > roomH && rvec[i] < roopH) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oo = rtemp;
		}

	}

}


void Print_Mass(double step, double *rvec, double partmassind[][3], double t, double *dpressvec) {

	double ind_oi, ind_oo, tav;	
	
	if(t==0) {
		tav = r_dze_o;		
	} else {

		int dim = find_num_zero(rvec,dpressvec);		// megnezi, hogy hany 0 pont van a derivaltban
		double r_count[dim];
		double temp_new = 0.;
		double temp = 0.;
		double rout = r_dze_o;

		int j, i;
		j=0;
	
		if(dim != 0) {						// ha van nullpont, akkor megkeresi, hogy hol
			for(i = 0; i < NGRID; i++) {
				temp_new = find_zero(i,rvec,dpressvec);	// ha van nullpont, akkor a temp_new valtozoba tarolja el --> ehhez vegig megy az egesz r-en vegig megy
// ez akkor kell, ha kulso es belso bump-ot is vizsgalunk, ekkor elmenti mindket nyomasi maximum helyet
/*				if(temp != temp_new && i > 3 && temp_new != 0.0) {
					printf("nullpont: %lg, j: %i\n",temp_new,j);
					r_count[j] = temp_new;
					j++;

				}
*/				
				if(temp_new > 0.) {			//	most csak a kulso bump-ot vizsgaljuk, elmentjuk a bump helyet
					temp = temp_new;
					rout = temp;
				} 
			}
		}

		tav = rout;
	}

	find_r_annulus(rvec,tav,&ind_oi,&ind_oo);		/*	A nyomasi maximum korul 2H tavolsagban keres korgyurut, a fuggveny visszaadja a cellak indexet	*/

	int index_i,i;
	double mass = 0.;

	for (i = 0; i < PARTICLE_NUMBER; i++) {
		index_i = partmassind[i][2];
		double masstemp = 0.;
		if ((index_i >= (int)ind_oi) && (index_i <= (int)ind_oo)) {
	
			masstemp = partmassind[i][1];
			mass = mass + masstemp;

		}
	}

	fprintf(massfil,"%lg %lg %lg\n",step,tav,mass);
	fflush(massfil);

}



void Print_Sigma(char *dens_name, double *rvec, double *sigmavec, double *pressvec, double *dpressvec) {

	int i;
	fmo = fopen(dens_name,"w");				/*	surface.dat, ez tartalmazza a korong kezdeti idopillanati adatait	*/

 	for(i = 1; i <= NGRID; i++) {
   		fprintf(fmo, "%lg   %lg   %lg   %lg\n", rvec[i], sigmavec[i],pressvec[i],dpressvec[i]);
	}

	fclose(fmo);

}


void Print_Pormozg_Size(char *size_name, int step, double rad[][2]){

	int i;
	fout3 = fopen(size_name,"w");

	for(i=0;i<PARTICLE_NUMBER;i++){
		fprintf(fout,"%lg %d %lg\n",(double)step,i,rad[i][0]);
		fprintf(fout3,"%lg  %lg  %lg\n",(double)step, rad[i][0], rad[i][1]*AU2CM);
	}

	fclose(fout3);

}


void Get_Sigma_P_dP(double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double deltat) {

	double u, u_bi, u_fi, sigma_temp[NGRID+2], uvec[NGRID+2];
	int i;

	sigma_temp[0] = sigmavec[0];
	sigma_temp[NGRID+1] = sigmavec[NGRID+1];
	uvec[0] = sigmavec[0] * visc(rvec[0]);
	uvec[NGRID+1] = sigmavec[NGRID+1] * visc(rvec[NGRID+1]);

	for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
		u = sigmavec[i] * visc(rvec[i]);
		u_bi = sigmavec[i-1] * visc(rvec[i-1]);
		u_fi = sigmavec[i+1] * visc(rvec[i+1]);
		uvec[i] = u;
		double temp = Coeff_1(rvec[i]) * (u_fi - 2.0 * u + u_bi) / (DD * DD) + Coeff_2(rvec[i]) * (u_fi - u_bi) / (2.0 * DD);
		sigma_temp[i] = uvec[i] + deltat * temp;

    
	}

	for(i = 1; i <= NGRID; i++) {
		sigmavec[i] = sigma_temp[i]/visc(rvec[i]);
		pressvec[i] = press(sigmavec[i],rvec[i]);
	}

	Perem(sigmavec);						/*	loading boundary condition in each timestep		*/
	dpress(dpressvec,pressvec);					/*	loading boundary condition in each timestep		*/
	Perem(pressvec);						/*	loading boundary condition in each timestep		*/
	Perem(dpressvec);						/*	loading boundary condition in each timestep		*/
	
}

void Get_Radius(char *nev, double radius[][2], double *pressvec, double *dpressvec, double *sigmavec, double *rvec, double *ugvec, double deltat, double t) {

	int i;
	double y, y_out, prad_new, particle_radius;	
	char scout[1024];

	for (i = 0; i < PARTICLE_NUMBER; i++) {	

		double drdt[PARTICLE_NUMBER];

		if (radius[i][0] > RMIN && radius[i][0] < RMAX) { 

			y = radius[i][0];
     			particle_radius = radius[i][1]; 

     			int_step(t,particle_radius,pressvec,dpressvec,sigmavec,rvec,ugvec,deltat,y,&y_out,&prad_new);

			if(t == 0) {
				drdt[i] = (fabs(y_out - y)/(deltat));
				radius[i][1] = radius[i][1];
				radius[i][0] = radius[i][0];
			} else {
				radius[i][1] = prad_new;
		     		radius[i][0] = y_out;
			}	
	
		} else {
			y_out = 0.0;
			radius[i][0] = 0.0;
		}

		
		if(t==0) {
	
			sprintf(scout,"%s/timescale.dat",nev);
			fout2 = fopen(scout,"w");			
/*	timescale.dat, ez tartalmazza a reszecskek bearamlasanak idoskalajat (evben)				*/
/*	A timescale.dat-hoz a t_drift = r / v(r,peak) kepletet a kovetkezo modon oldja meg a program:		*/
/*	Kiszamolja, hogy a kovetkezo idopillanatban hol lenne minden reszecske, majd az uj koordinatabol (x(1))	*/
/*	a regit (x(0)) levonva, majd az egeszet elosztva dt-vel megkapjuk a bearamlas sebesseget (drdt)		*/
/*	ez utan kiszamolhato az eredeti t_drift keplet, ezt iratom ki a file-ba evben (ezert osztok le 2PI-vel	*/

			fprintf(fout2,"%lg %lg\n",radius[i][0], (radius[i][0] / drdt[i])/2.0/M_PI);
			fclose(fout2);
		}
	}

}


void Count_Mass(double step, double radius[][2], double partmassind[][3], double *massvec, double *dpressvec, double *rvec, double t) {

	int i, rindex;
	double rmid;	

	for (i = 0; i < PARTICLE_NUMBER; i++) {	

  		rmid = (radius[i][0] - RMIN) / DD;     						/* 	The integer part of this gives at which index is the body			*/
		rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/

		partmassind[i][0] = (double)i;							/*	index of the particles				*/
		partmassind[i][1] = massvec[i];							/*	mass of the particles				*/
		partmassind[i][2] = rindex;							/*	initial distance of the particles				*/
	
	}

	Print_Mass(step,rvec,partmassind,t,dpressvec);

}


void Mk_Dir(char *nev) {

	int i, step, mappa;
	char comm[1024];

/*	A kimeneti file-ok szamara egy tarolo mappa letrehozasa (system(mkdir ...)) paranccsal	*/
	mappa=system("mkdir output");

/*	A ciklus ellenorzi, hogy letezik-e mar ez adott mappa, ha nem, letrehozza output neven...	*/
		if (mappa == 0) {
			printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa elkeszult ...\n\n\n");
					sprintf(nev,"output");				/*	A mappa nevenek eltarolasa azert, hogy a file-okat az eppen aktualis mappaba tudja kesobb kiirni a program	    */

/*	...ha letezik az output mappa, akkor egy do-while ciklussal szamozott mappat hoz letre (pl. output.0), es addig fut a ciklus, mig aztan nem talalja mar a soron kovetkezo szamozott mappat. Ekkor letre hozza azt	*/
		} else {
			printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa letrehozasa meghiusult ...\n");
			printf("... A mappa valoszinuleg mar letezik, uj mappa letrehozasa ...\n\n\n");

			do{
				for(i=0;i<=step;i++){
					sprintf(nev,"output.%i",i);			/*	A mappa nevenek eltarolasa azert, hogy a file-okat az eppen aktualis mappaba tudja kesobb kiirni a program	    */
					sprintf(comm,"mkdir output.%i",i);	
				}

				mappa=system(comm);
				step++;					

			} while (mappa!=0);

			printf("... A kimeneti file-ok tarolasahoz a(z) %s mappa elkeszult ...\n\n\n",nev);
		}

}






/*	Itt vegzi el az integralast, ha szukseg van ra	*/
void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec) {

   	int linesout;
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	PARTICLE_NUMBER = reszecskek_szama(linesout);  		

	char dens_name[1024],size_name[1024], porout[1024], massout[1024], mv[1024];
   	double radius[PARTICLE_NUMBER][2],radius_rec[PARTICLE_NUMBER][2];
	double massvec[PARTICLE_NUMBER];
	double max, min;

/*	A reszecskek adatainak beolvasasa file-bol		*/
	por_be(radius,massvec);				

	int dummy,i;
	sprintf(mv,"mv %s %s/",filenev1,nev);
	dummy = system(mv);	


	snprintf(porout,1024,"%s/pormozgas.dat",nev);
	snprintf(massout,1024,"%s/mass.dat",nev);
   	fout = fopen(porout,"w");
   	massfil = fopen(massout,"w");

  	double t = 0.0;
   	double t_integration = TMAX * 2.0 * M_PI; 						/*	numerikus integralas idotartama		*/
	double deltat = time_step(rvec)/5.;							/*	idolepes	*/
	double partmassind[PARTICLE_NUMBER][3];
	double L = 0.;
	double lep = 0.;


   	do {



		for (i=0; i < PARTICLE_NUMBER; i++) {
			if (radius[i][0] > 0. && radius[i][0] > RMIN) {
				radius_rec[i][0] = 1. / radius[i][0];
			} else {
				radius_rec[i][0] = 0.;
			}
		}

		max = find_max(radius);					/*	Megkeresi, hogy melyik a legtavolabbi reszecske a kozponti csillagtol	*/
		min = find_max(radius_rec);
		min = 1. / min;

		if(max >= RMIN && min != max) {						/*	Ha a legtavolabbi porreszke merete nagyobb RMIN-nel, akkor folytatodik az integralas, kulonben kilep	*/

			double time = t / 2.0 / M_PI;
		
			if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){	

				printf("max: %lg, min: %lg\n", max, min);						

				if (t==0) {
					snprintf(dens_name,1024,"%s/surface.dat",nev);
				} else {
					snprintf(dens_name,1024,"%s/dens.%d.dat",nev,(int)L);

				}

				snprintf(size_name,1024,"%s/size.%d.dat",nev,(int)L);
	
				Print_Sigma(dens_name, rvec, sigmavec, pressvec, dpressvec);
				Print_Pormozg_Size(size_name,(int)L,radius);
				Count_Mass(L,radius,partmassind,massvec,dpressvec,rvec,t);

				fflush(fout);
				L = L+(double)(TMAX/WO);
		
			}

			Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
			Get_Radius(nev,radius,pressvec,dpressvec,sigmavec,rvec,ugvec,deltat,t); 
			
			t = t + deltat;						/*	Idoleptetes		*/

		} else {				/*	Ha a legmesszebbi reszecske tavolsaga mar nem nagyobb, vagy egyenlo, mint RMIN, akkor a program "figyelmezteto szoveg" mellett sikeresen kilep, nem fut "feleslegesen" tovabb.	*/
			printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",L);
			exit(EXIT_SUCCESS);
		}


   	} while (t <= t_integration);
  

/*	Az idoleptetes leteltevel a program sikeresen kilep	*/
	printf("\n\nA program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n",max); 


}


/*	Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok	*/
void infoCurrent(char *nev) {


	char out[1024];

	sprintf(out,"run_%i.dat",(int)TCURR);
	jelfut = fopen(out,"w");
	
	fprintf(jelfut,"A jelenlegi futás a %s mappaban taláható!\n",nev);
	fprintf(jelfut,"\n\nA korong paraméterei:\nRMIN: %lg, RMAX: %lg\nSIGMA0: %lg, SIGMA_EXP: %lg, flaring index: %lg\nALPHA_VISC: %lg, ALPHA_MOD: %lg\nR_DZE_I: %lg, R_DZE_O: %lg, DR_DZEI: %lg, DR_DZE_O: %lg  (*** R_DZE_I/O = 0, akkor azt a DZE-t nem szimulálja a futás! ***)\n\n\n",RMIN,RMAX,SIGMA0,SIGMAP_EXP,FLIND,alpha_visc,a_mod,r_dze_i,r_dze_o,Dr_dze_i,Dr_dze_o);
	fprintf(jelfut,"A központi csillag tömege: %lg M_Sun\n",STAR);
	fclose(jelfut);

}


int main() {
   
   	double sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], ugvec[NGRID+2];
	char nev[1024], mv[1024];
	char dens_name[1024];
   

/*	A korong parametereit beolvassa: a sigma profil kitevojet es a sigma0-t					*/
	disk_param_be(&SIGMA0, &SIGMAP_EXP, &RMIN, &RMAX, &r_dze_i, &r_dze_o, &Dr_dze_i, &Dr_dze_o, &a_mod, &PDENSITY, &PDENSITYDIMLESS, &alpha_visc,&STAR,&FLIND);

	DD = (RMAX - RMIN) / (NGRID - 1);						/*	rácsfelbontás	*/

/*	Kezdeti profilok betoltese a megadott vektorokba	*/
	load_R(rvec);
	Initial_Profile(sigmavec,rvec);
	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	Initial_Ugas(sigmavec,rvec,ugvec);
	timePar(&TMAX,&WO,&TCURR);

	Mk_Dir(nev);								/*	mappa letrehozasa az adatok eltarolasahoz	*/

	int dummy;
	sprintf(mv,"mv %s %s/",filenev2,nev);
	dummy = system(mv);	
	sprintf(mv,"mv %s %s/",filenev3,nev);
	dummy = system(mv);

	infoCurrent(nev);

	double opt = 1.;



	if(opt == 0.) {

/*	t=0-ban kiirja a sigma-t, a nyomast es a nyomasderivaltjat	*/
		snprintf(dens_name,1024,"%s/surface.dat",nev);
		Print_Sigma(dens_name,rvec,sigmavec,pressvec,dpressvec);
	} 
	
	if(opt == 1.) {

		tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec);
	}

	return 0;

}
