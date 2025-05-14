#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TWOPI        	2.0 * M_PI
#define G            	0.01720209895  //Gaussian grav. constant
#define G2		G * G
#define STAR		1.0
#define HASP		0.05
#define SIGMAP_EXP  	-1.0
#define SIGMA0      	100 //g/cm2
#define alpha_visc 	0.001				/*	alpha parameter				*/

#define SDCONV      	1.12521e-7
#define AU2CM		1.496e13
#define PDENSITY    	4.0
#define SUN2GR	    	1.989e33

#define NGRID       	2000
#define RMIN		1.0
#define RMAX		100.0
#define DGRID 		(RMAX - RMIN) / (NGRID - 1)

#define filenev		"init_data.dat"

int NV, PARTICLE_NUMBER;
double DD = DGRID;	/*	valamiert ha az interpolalasnal csak siman DGRID-et hasznalok, akkor nem jol szamol!	*/

FILE *fin1, *fout, *fmo, *fil; 

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


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma) {		/*	in the Epstein drag regime	*/

	return PDENSITY * pradius * M_PI / (2.0 * sigma);

} 


/*	A feluletisuruseg erteke egyelore cgs-ben		*/
double sigmagas_mmsn_cgs(double r) {

	return SIGMA0*pow(r,SIGMAP_EXP);

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

	return v_kep(r) * scale_height(r);

}


double visc(double r){
 
  	double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);
	
 	r_dze_i = 5.0;					
 	r_dze_o = 20.0;						/*	az atmenet 2*h	*/
  	Dr_dze_i = 0.25;
  	Dr_dze_o = 1.0;
  	a_mod = 0.01;

  	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));	/* 	10 AU-nal a_mod-dal megvaltoztatja a viszkozitas merteket		*/
//	alpha_r = 1.0;

    	nu = alpha_visc * alpha_r * cs * H;
  
  	return nu;
  
}


double Coeff_1(double r){					/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
  
  	double A;
  	A = 3.0 * visc(r);
  	return A;
  
}


double Coeff_2(double r){					/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/		
  
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


/* 	local pressure of the gas	*/
double press(double sigma, double r) {

	return HASP * G2 * STAR * sigma / (sqrt(TWOPI * r * r));

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


/*	a nyomas analitikus derivaltja	*/
void dpress(double *dp, double *p) {

	int i;

	for(i = 1; i <= NGRID; i++) {	
		dp[i] = (p[i+1] - p[i-1]) / (2.0 * DGRID);
	}

//	return (SIGMAP_EXP - 2.0) * HASP * G2 * STAR * SIGMA0 * pow(r, (SIGMAP_EXP - 3.0)) / sqrt(TWOPI);

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
void int_step(double prad, double *dpressvec, double *sigmavec, double *rvec, double step, double y, double *ynew)
{
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


int main() {
   
   	int i, L, dummy, linesout, PARTICLE_NUMBER;
   	double distance, reprmass, y, y_out, t, t_integration, deltat, particle_radius, sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2];
	char dens_name[64];
   
   	PARTICLE_NUMBER = 0;
	linesout = 0;

	PARTICLE_NUMBER = reszecskek_szama(linesout);  		/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/

   	double radius[PARTICLE_NUMBER][2];

   	fout = fopen("pormozgas.dat","w");
   	fin1 = fopen(filenev,"r");
   	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/
            	if(fscanf(fin1,"%d %lg %lg %lg",&dummy,&distance,&particle_radius,&reprmass) == 4) {	/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 4, mert 4 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/
            		radius[i][0] = distance;
	   		radius[i][1] = particle_radius;
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  Nem sikerult a %i-ik sort beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n",dummy);
			getchar();
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(fin1);
	
	printf("\n\n *******   A file beolvasasa sikerult!   ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");	
	getchar();

	for(i = 0; i <= NGRID+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = RMIN + (i-1) * DD;
	}

	Initial_Profile(sigmavec,rvec);
	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	
  	t = 0.0;
   	t_integration = 5000.0 * 365.25; //numerikus integralas idotartama
//	deltat = time_step(rvec)/5.;
	deltat = 0.1;	

  	fmo = fopen("surface.dat","w");

 	for(i = 1; i <= NGRID; i++) {

   		fprintf(fmo, "%lg   %lg   %lg   %lg\n", rvec[i], sigmavec[i],pressvec[i],dpressvec[i]); 		/*	print sigma in a file	*/

 	}

	fclose(fmo);

   	L = 0;

   	do {

  		for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
  			sigmavec[i] = (sigmavec[i] * visc(rvec[i]) + deltat * (Coeff_1(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - 2.0 * sigmavec[i] * visc(rvec[i]) + sigmavec[i-1] * visc(rvec[i-1])) / (DD * DD) + Coeff_2(rvec[i]) * (sigmavec[i+1] * visc(rvec[i+1]) - sigmavec[i-1] * visc(rvec[i-1])) / (2.0 * DD))) / visc(rvec[i]);
			pressvec[i] = press(sigmavec[i],rvec[i]);
    
  		}

  		Perem(sigmavec);					/*	loading boundary condition in each timestep		*/
		dpress(dpressvec,pressvec);
		Perem_Press(pressvec);
		Perem_dPress(dpressvec);

		for (i = 0; i < PARTICLE_NUMBER; i++) {
		
			if (radius[i][0] > RMIN) { 
				y = radius[i][0];
		     		particle_radius = radius[i][1]; 
//				particle_radius = 1.0;

		     		int_step(particle_radius,dpressvec,sigmavec,rvec,deltat,y,&y_out);
		     		
				radius[i][0] = y_out;
			} else {
				radius[i][0] = 0.0;
			}
		     
		     	if (L%15000==0) {
//				printf("%lg %d %lg\n",t/365.25,i,radius[i][0]);
		         	fprintf(fout,"%lg %d %lg\n",t/365.25,i,radius[i][0]);
		         	fflush(fout);

		     	}

		}
		

		if (L%15000==0) {
			snprintf(dens_name,64,"dens_%lg.dat", t/365.25);
			fil=fopen(dens_name,"w");
	
	  		for(i = 1; i <= NGRID; i++) {							
				fprintf(fil,"%lg %lg %lg %lg\n",rvec[i], sigmavec[i], pressvec[i], dpressvec[i]);
			}
	
			fclose(fil);

			printf("t: %lg\n",t/365.25);
		
		}

		t = t + deltat;
		L++;

   	} while (t < t_integration);
   
	return 0;

}
