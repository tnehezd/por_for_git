#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TWOPI        	2.0 * M_PI
#define G            	0.01720209895  //Gaussian grav. constant
#define G2		G * G
#define STAR		1.0
#define ALPHA		0.996
#define SDCONV      	1.12521e-7
#define HASP		0.05
#define AU2CM		1.496e13
#define PDENSITY    	4.0
#define SUN2GR	    	1.989e33
#define NGRID       	5000
#define SIGMAP_EXP  	-1.5
#define SIGMA0      	1700 //g/cm2

#define filenev		"init_data.dat"

int NV, PARTICLE_NUMBER;
double RD[NGRID], SDISK[NGRID], PDISK[NGRID], DPDRDISK[NGRID], DD, RMIN, RMAX;

FILE *fin1, *fout; 

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


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius,double sigma) {		/*	in the Epstein drag regime	*/

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


/*	local sound speed		*/
double c_sound(double r) {

	return v_kep(r) * scale_height(r);

}


/* 	local pressure of the gas	*/
double press(double r) {

	double Sigma;
	Sigma = sigmagas_mmsn_cgs(r);
	return HASP * G2 * STAR * Sigma / (sqrt(2.0 * M_PI * r * r));

}


/*	a nyomas analitikus derivaltja	*/
double dpress(double r) {

	return (SIGMAP_EXP - 2.0) * HASP * G2 * STAR * SIGMA0 * pow(r, (SIGMAP_EXP - 3.0)) / sqrt(2.0 * M_PI);

}	

/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double r, double *drdt) {
     
     	double P, H, dPdr, Sigma, St, csound;
    
	Sigma = sigmagas_mmsn_cgs(r); 
	St = Stokes_Number(pradius,Sigma);
	H = scale_height(r);  
	P = press(r);
	dPdr = dpress(r);
	csound = c_sound(r); 

     	*drdt = St / (1 + St * St) * H / P * dPdr * csound;
     
}

//megiscsak vektorokkal kell majd csinalni? -- mit? :D -- ja! :D
void int_step(double prad, double step, double y, double *ynew)
{
	double dy1,dy2,dy3,dy4;
	double ytemp;
	
	eqrhs(prad, y, &dy1);
	
	ytemp = y + 0.5 * step * dy1;
	eqrhs(prad, ytemp, &dy2);
		
	ytemp = y + 0.5 * step * dy2;
	eqrhs(prad, ytemp, &dy3);
	
	ytemp = y + step * dy3;
	eqrhs(prad, ytemp, &dy4);
	
	*ynew = y + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
	
}


int main() {
   
   	int i, L, dummy, linesout, PARTICLE_NUMBER	;
   	double distance, reprmass, y, y_out, t, t_integration, deltat, particle_radius;
   
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

	
  	t = 0.0;
   	t_integration = 100000.0 * 365.25; //numerikus integralas idotartama
   	deltat = 0.1; //idolepes napban

   	L = 0;

   	do {
		
		for (i = 0; i < PARTICLE_NUMBER; i++) {
		
			y = radius[i][0];
		     	particle_radius = radius[i][1]; 
		     
		     	int_step(particle_radius,deltat,y,&y_out);
		     	radius[i][0] = y_out;
		     
		     	if (L%100000==0) {
				printf("%lg %d %lg\n",t/365.25,i,radius[i][0]);
		         	fprintf(fout,"%lg %d %lg\n",t/365.25,i,radius[i][0]);
		         	fflush(fout);
		     	}
		}
		
		t = t + deltat;
		L++;
   	} while (t < t_integration);
   
	return 0;

}
