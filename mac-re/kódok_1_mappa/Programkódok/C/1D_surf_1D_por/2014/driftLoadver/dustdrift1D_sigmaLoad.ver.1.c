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

#define filenev		"dens_1D_100.dat"

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


/*	Statikus sigma profil beolvasasa	*/
void sigma_be(double *sigma, double *radius) {

	double rad, sig;
	int dummy = 0, i;
   	fin1 = fopen(filenev,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%lg %lg",&rad,&sig) == 2) {
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 2, mert 2 oszlopot szeretnénk beolvasni. 	*/

            		radius[i+1] = rad;
	   		sigma[i+1] = sig;		
			dummy++;

	    	} else {

/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/					
			printf("\n\n*******************     ERROR!     *********************\n\n  Nem sikerult a %i-ik sort beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n",dummy);
//			getchar();
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(fin1);
	
	printf("\n\n *******   A file beolvasasa sikerult!   ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");	

}


/*	local scale height	*/
double scale_height(double r) {

	return r * HASP;

}



/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r) {

	return sqrt(G2 * STAR / r / r / r);			/*	v_kepler in AU / (yr/2pi)	*/

}


/*	local sound speed		*/
double c_sound(double r) {

	return kep_freq(r) * scale_height(r);

}



/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {

	return 1. / (2.0 * M_PI) * sigma / (HASP * r);

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

//	Parabola(vec, NGRID - 2, NGRID - 1, NGRID, &a, &b, &c, DD, RMIN);
	Parabola(vec, PARTICLE_NUMBER - 2, PARTICLE_NUMBER - 1, PARTICLE_NUMBER, &a, &b, &c, DD, RMIN);
//	vec[NGRID+1] = a * (RMAX + DD) * (RMAX + DD) + b * (RMAX + DD) + c;
	vec[PARTICLE_NUMBER+1] = a * (RMAX + DD) * (RMAX + DD) + b * (RMAX + DD) + c;

	
}





/* 	local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r) {

	return rho_mp(sigma,r) * c_sound(r) * c_sound(r);

}


void load_Press(double *pressvec, double *sigmavec, double *rvec, double DD, double RMIN, double RMAX, int PARTICLE_NUMBER){		/*	initial profile of pressure		*/

  	int i;
  
//  	for(i = 1; i < NGRID; i++) {
	for(i=0; i < PARTICLE_NUMBER+2; i++) {
    		pressvec[i] = press(sigmavec[i],rvec[i]);		
  	}
  
//  	Perem(pressvec,DD,RMIN,RMAX);
}


/*	a nyomas derivaltja	*/
void dpress(double *dp, double *p, double DGRID, int PARTICLE_NUMBER) {

	int i;
	double ptemp, pvec[PARTICLE_NUMBER+1];

//	for(i = 1; i <= NGRID; i++) {
	for(i = 1; i <= PARTICLE_NUMBER; i++) {
		ptemp = (p[i+1] - p[i-1]) / (2.0 * DGRID);
		pvec[i] = ptemp;
	}
	
//	for(i = 1; i <= NGRID; i++) {
	for(i = 1; i <= PARTICLE_NUMBER; i++) {
		dp[i] = pvec[i];
	}


}	


void load_dPress(double *dpressvec, double *pressvec, double DD, double RMIN, double RMAX, int pnum){		/*	initial profile of pressure		*/

	dpress(dpressvec,pressvec,DD,pnum);
   	Perem(dpressvec,DD,RMIN,RMAX);

}



int main() {
   
   	int i, L, linesout, dim, j;
/*   	double y, y_out, t, t_integration, deltat, particle_radius, sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], u, u_bi, u_fi, sigma_temp[NGRID+2], sigma, uvec[NGRID+2];
	double max;
	char dens_name[64],size_name[64];
*/   
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	PARTICLE_NUMBER = reszecskek_szama(linesout); 
	double sigmavec[PARTICLE_NUMBER+2], rvec[PARTICLE_NUMBER+2], pressvec[PARTICLE_NUMBER+2], dpressvec[PARTICLE_NUMBER+2]; 
	sigma_be(sigmavec,rvec);

	RMIN = rvec[1];
	RMAX = rvec[PARTICLE_NUMBER-1];	
	double DD = (RMAX - RMIN) / (PARTICLE_NUMBER - 1);

	printf("RMIN: %lg, RMAX: %lg, DD: %lg\n", RMIN, RMAX, DD);

	Perem(rvec,DD,RMIN,RMAX);
	Perem(sigmavec,DD,RMIN,RMAX);
	
	load_Press(pressvec, sigmavec, rvec, DD, RMIN, RMAX, PARTICLE_NUMBER);
	load_dPress(dpressvec, pressvec, DD, RMIN, RMAX, PARTICLE_NUMBER);
	

	fout = fopen("out.dat","w");

	for(i=0; i<=PARTICLE_NUMBER; i++) fprintf(fout,"%lg %lg %lg %lg\n",rvec[i],sigmavec[i],pressvec[i], dpressvec[i]);

	fclose(fout);

	return 0;

}
