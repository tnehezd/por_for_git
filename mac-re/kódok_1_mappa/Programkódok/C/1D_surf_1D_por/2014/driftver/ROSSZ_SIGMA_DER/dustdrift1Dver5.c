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
#define alpha_visc 	0.001				/*	alpha parameter				*/

#define SUN2GR	    	1.989e33
#define SDCONV      	1.12521e-7
#define AU2CM		1.496e13
#define PDENSITY    	4.0				/*	a reszecskek atlag surusege		*/
#define PDENSITYDIMLESS	PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM 


#define NGRID       	2000
#define RMIN		1.0
#define RMAX		100.0
#define DGRID 		(RMAX - RMIN) / (NGRID - 1)

#define STEPWRITE	5000				/*	hanyadik lepesenkent irassa ki az adatokat	*/

#define filenev		"init_data.dat"

int PARTICLE_NUMBER;
double DD = DGRID;	/*	valamiert ha az interpolalasnal csak siman DGRID-et hasznalok, akkor nem jol szamol!	*/
double SIGMAP_EXP, SIGMA0;


FILE *fin1, *fin2, *fout, *fmo, *fil, *fout2; 


/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int reszecskek_szama(int numout){

	char c;
	fin1 = fopen(filenev,"r+");		
	numout = 0;

/*	A p	orreszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal				*/

		while((c = fgetc(fin1)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 							*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert										*/
	fclose(fin1);	

	return numout;

}


/*	A korong parametereinek beolvasasa	*/
void disk_param_be(double *sigma0, double *sdexp) {

	double dummy;
	int dummy2;

	double sig0, exp;

	fin2 = fopen("disk_param.dat","r");

           	if(fscanf(fin2,"%lg  %lg %d %lg %lg",&dummy,&dummy,&dummy2,&exp,&sig0) == 5) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 5, mert 5 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  sigma0: %lg, sdexp: %lg\n\n\n", sig0,exp);  
			*sigma0 = sig0;
			*sdexp = exp;        
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			getchar();
			exit(EXIT_FAILURE);
   	        }
	fclose(fin2);
	
}


/*	A porreszecskek adatainak beolvasasa	*/
void por_be(double radius[][2]) {

	int i, dummy;
	double distance, particle_radius, reprmass;
	
   	fin1 = fopen(filenev,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%d %lg %lg %lg",&dummy,&distance,&reprmass,&particle_radius) == 4) {	
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 4, mert 4 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/

            		radius[i][0] = distance;
	   		radius[i][1] = particle_radius / AU2CM;		/* a részecske mérete AU-ban	*/
	    	} else {

/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/					
			printf("\n\n*******************     ERROR!     *********************\n\n  Nem sikerult a %i-ik sort beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n",dummy);
			getchar();
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(fin1);
	
	printf("\n\n *******   A file beolvasasa sikerult!   ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");	
	getchar();

}


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma) {		/*	in the Epstein drag regime	*/

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


/*	Lokalis viszkozitas erteke	*/
double visc(double r){
 
  	double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);
	
 	r_dze_i = 5.0;					
 	r_dze_o = 15.0;						/*	az atmenet 2*H	*/
  	Dr_dze_i = 2. * scale_height(r_dze_i);
  	Dr_dze_o = 2. * scale_height(r_dze_o);
  	a_mod = 0.01;

/*	alpha_r: a redukcio merteke	*/
 	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));	/* 	10 AU-nal a_mod-dal megvaltoztatja a viszkozitas merteket		*/
//	alpha_r = 1.0;

    	nu = alpha_visc * alpha_r * cs * H;
  
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


void Perem(double *sigmavec) {					/*	boundary condition for sigma		*/
  
  	sigmavec[0] = (sigmavec[1] / visc(RMIN)) * visc(RMIN - DD);
  	sigmavec[NGRID+1] = (sigmavec[NGRID] / visc(RMIN + (NGRID) * DD)) * visc(RMIN + (NGRID+1) * DD);
  
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

	return 1. / (2.0 * M_PI) * sigma / (HASP * r);

}


/* 	local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r) {

	return rho_mp(sigma,r) * c_sound(r) * c_sound(r);

}


void Perem_Press(double *pressvec) {					/*	boundary condition for pressure		*/
  
  	pressvec[0] = pressvec[1];
  	pressvec[NGRID+1] = pressvec[NGRID];
  
}


void Initial_Press(double *pressvec, double *sigmavec, double *rvec){		/*	initial profile of pressure		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		pressvec[i] = press(sigmavec[i],rvec[i]);		
  	}
  
  	Perem_Press(pressvec);
}


/*	a nyomas derivaltja	*/
void dpress(double *dp, double *p) {

	int i;

	for(i = 1; i <= NGRID; i++) {	
		dp[i] = (p[i+1] - p[i-1]) / (2.0 * DGRID);
	}


}	


void Perem_dPress(double *dpressvec) {					/*	boundary condition for pressure		*/
  
  	dpressvec[0] = dpressvec[1];
  	dpressvec[NGRID+1] = dpressvec[NGRID];
  
}


void Initial_dPress(double *dpressvec, double *pressvec){		/*	initial profile of pressure		*/

	dpress(dpressvec,pressvec);
   	Perem_dPress(dpressvec);

}


/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double r, double *drdt) {
     
     	double P, H, dPdr, St, csound;
    
	St = Stokes_Number(pradius,sigma);
	H = scale_height(r);  
	P = press(sigma,r);
	dPdr = dp;
	csound = c_sound(r); 

     	*drdt = St / (1 + St * St) * H / P * dPdr * csound;
     
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


//megiscsak vektorokkal kell majd csinalni? -- mit? :D -- ja! :D
/*	Runge-Kutta4 integrator	*/
void int_step(double prad, double *dpressvec, double *sigmavec, double *rvec, double step, double y, double *ynew) {
	double dy1,dy2,dy3,dy4;
	double ytemp;
	double sigma, dpress; 

	interpol(sigmavec,rvec,y,&sigma);
	interpol(dpressvec,rvec,y,&dpress);
	
	eqrhs(prad, dpress, sigma, y, &dy1);
	
	ytemp = y + 0.5 * step * dy1;
	eqrhs(prad, dpress, sigma, ytemp, &dy2);
		
	ytemp = y + 0.5 * step * dy2;
	eqrhs(prad, dpress, sigma, ytemp, &dy3);
	
	ytemp = y + step * dy3;
	eqrhs(prad, dpress, sigma, ytemp, &dy4);

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


int main() {
   
   	int i, L, linesout;
   	double y, y_out, t, t_integration, deltat, particle_radius, sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2];
	double max;
	char dens_name[64];
   
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	PARTICLE_NUMBER = reszecskek_szama(linesout);  		

/*	A korong parametereit beolvassa: a sigma profil kitevojet es a sigma0-t					*/
	disk_param_be(&SIGMA0, &SIGMAP_EXP);			

   	double radius[PARTICLE_NUMBER][2];
	double drdt[PARTICLE_NUMBER];


/*	A reszecskek adatainak beolvasasa file-bol								*/
	por_be(radius);				


 	for(i = 0; i <= NGRID+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = RMIN + (i-1) * DD;
	}

/*	Kezdeti profilok betoltese a megadott vektorokba	*/
	Initial_Profile(sigmavec,rvec);
	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	
  	t = 0.0;
   	t_integration = 500000.0 * 2.0 * M_PI; 					/*	numerikus integralas idotartama		*/
	deltat = time_step(rvec)/5.;						/*	idolepes	*/
	

   	L = 0;

   	fout = fopen("pormozgas.dat","w");

   	do {

		if(t == 0.0) {							/*	A t=0 idopillanatban elkesziti az alabbi file-okat:	*/
		
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

			for (i = 0; i < PARTICLE_NUMBER; i++) {			

				if (radius[i][0] > RMIN && radius[i][0] < RMAX) { 

					y = radius[i][0];
		     			particle_radius = radius[i][1]; 

		     			int_step(particle_radius,dpressvec,sigmavec,rvec,deltat,y,&y_out);

					drdt[i] = (fabs(y_out - y)/(deltat));
	
				} else {

					y_out = 0.0;
				}

				fprintf(fout2,"%lg %lg\n",radius[i][0], (radius[i][0] / drdt[i])/2.0/M_PI);

			}
	
			fclose(fout2);

		}

		max = find_max(radius);					/*	Megkeresi, hogy melyik a legtavolabbi reszecske a kozponti csillagtol	*/

		if(max >= RMIN) {							/*	Ha a legtavolabbi porreszke merete nagyobb RMIN-nel, akkor folytatodik az integralas, kulonben kilep	*/

  			for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
  				sigmavec[i] = (sigmavec[i] * visc(rvec[i]) + deltat * (Coeff_1(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - 2.0 * sigmavec[i] * visc(rvec[i]) + sigmavec[i-1] * visc(rvec[i-1])) / (DD * DD) + Coeff_2(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - sigmavec[i-1] * visc(rvec[i-1])) / (2.0 * DD))) / visc(rvec[i]);
				pressvec[i] = press(sigmavec[i],rvec[i]);
    
  			}

  			Perem(sigmavec);						/*	loading boundary condition in each timestep		*/
			dpress(dpressvec,pressvec);					/*	loading boundary condition in each timestep		*/
			Perem_Press(pressvec);						/*	loading boundary condition in each timestep		*/
			Perem_dPress(dpressvec);					/*	loading boundary condition in each timestep		*/
	
			t = t + deltat;							/*	Idoleptetes	*/
			L++;								/*	Idolepes szamlalo	*/
	
			for (i = 0; i < PARTICLE_NUMBER; i++) {	
		
				if (radius[i][0] >= RMIN && radius[i][0] <= RMAX) { 	/*	Ha a reszecske benne RMIN es RMAX kozott van, akkor megcsinalja az integralast	*/

					y = radius[i][0];
			     		particle_radius = radius[i][1]; 
//					particle_radius = 5.0 / AU2CM;

			     		int_step(particle_radius,dpressvec,sigmavec,rvec,deltat,y,&y_out);	/*	Integralas	*/

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
				fil=fopen(dens_name,"w");
		
		  		for(i = 1; i <= NGRID; i++) {							
					fprintf(fil,"%lg %lg %lg %lg %lg\n",rvec[i], sigmavec[i], pressvec[i], dpressvec[i],c_sound(rvec[i]));
				}
	
				fclose(fil);

				printf("t: %lg\n",t / (2.0 * M_PI));
		
			}


		} else {				/*	Ha a legmesszebbi reszecske tavolsaga mar nem nagyobb, vagy egyenlo, mint RMIN, akkor a program "figyelmezteto szoveg" mellett sikeresen kilep, nem fut "feleslegesen" tovabb.	*/
			printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",t/2./M_PI);
			getchar();
			exit(EXIT_SUCCESS);
		}

   	} while (t < t_integration);
  

/*	Az idoleptetes leteltevel a program sikeresen kilep	*/
	printf("\n\nA program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n",max); 
	getchar();


	return 0;

}
