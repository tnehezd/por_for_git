#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

/*	BUG: KISEBB IDOLEPESNEL	MAS A NOVEKEDESI GORBE
	     KETPOPULACIOS ESETBEN NEM EGYEZIK ME A GORBE A 15-OS VERZIOVAL	*/


/*------------------------------------------------------------------------------*/
/*					MAKROK					*/
/*------------------------------------------------------------------------------*/

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


/*	beolvasando file-ok nevei			 */
#define filenev1	"init_data.dat"
#define filenev2	"disk_param.dat"
#define filenev3	"time.dat"


/*------------------------------------------------------------------------------*/
/*				GLOBALIS VALTOZOK				*/
/*------------------------------------------------------------------------------*/


double RMIN, RMAX;					/*	korong kulso es belso hatara				*/
int PARTICLE_NUMBER, NGRID;				/*	a reszecskek szama, illetve a gridek szama		*/
double SIGMAP_EXP, SIGMA0;				/*	a felületi suruseg kitevoje, illetve erteke r=1AU-nal	*/
double FLIND;						/*	flaring index						*/			
double DD;						/*	grid felbontas						*/

double r_dze_i, r_dze_o; 				/*	a DEADZONE belso es kulso hatara			*/
double Dr_dze_i, Dr_dze_o;				/*	a DEADZONE hataranak atmenete				*/

double a_mod;						/*	a viszkozitas csokkentesenek merteke a DZ-ban		*/
double alpha_visc;					/*	a viszkozitas erteke					*/


double PDENSITY;					/*	a reszecskek atlagsurusege cgs-ben			*/
double PDENSITYDIMLESS;					/*	a reszecskek atlagsurusege dimenziotlan egysegekben = PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;	*/

double STAR;						/*	a kozponti csillag tomege naptomegben			*/
double TMAX;						/*	integralasi ido vege					*/
double WO;						/*	a file-okat adott idokozonkent irja ki a progam, ez ahhoz szukseges parameter	*/
double TCURR;						/*	a futas inditasa idopillanata --> kimenetnel erdekes	*/
double DT;						/*	idolepes beolvasasahoz szukseges valtozo		*/


double optev, optdr, optgr, opttwopop, optdze, optinp;	/*	parameterek, amelyekkel ki es be lehet kapcsolni a szamolasok egy reszet:
									optev: feluleti suruseg fejlodese (ha optev == 1)
									optdr: van-e drift a futasban
									optgr: van-e reszecskenovekedes a futasban
									opttwopop: 2 populaciora szamolunk-e
									optdze: hany deadzone edge van a futasban
									optinp: van-e bemeneti file a sigma-ra		*/
double fFrag, uFrag;					/*	a reszecske novekedeshez szukseges parameterek Birnstiel et al 2012 alapjan:
									fFrag: ?
									uFrag: fragmentacios sebesseg cgs-ben		*/

const char *inputsig;


FILE *fin1, *fin2, *fout, *foutmicr, *fmo, *fil, *fout2, *massfil, *fout3, *jelfut; 


/*	A futasok soran valtoztathato opciok letrehozasa egy strukturaval	*/
typedef struct options {

	// Option for dust drift
	double drift;
	// Option for dust growth
	double growth;
	// Option for solving the duffusion equation of the surface density
	double evol;
	// Option fot two population model
	double twopop;
	// Fragmentation barrier in cm/s
	double ufrag;
	// Fragmentation barrier constant
	double ffrag;
	// Number of grid cells
	int ngrid;
	// Input sigma file
	const char *input; 
	// Option for time stepping
	double tStep;

} options_t;


/*	A struktura elemeinek feltoltese alapertelmezett ertekekkel	*/
void create_default_options(options_t *opt) {

	opt->drift		 = 1.;		//Dust drift is included
	opt->growth		 = 1.;		//Particle growth is included
	opt->evol		 = 1.;		//The evolution of surface density is included
	opt->twopop		 = 1.;		//Two population simulation is included
	opt->ufrag		 = 1000.0;	//Fragmentation velocity in CGS -- Birnstiel et al 2012
	opt->ffrag		 = 0.37;	//?
	opt->ngrid		 = 2000;	//Number of the grid cells
	opt->input		 = "";
	opt->tStep		 = 0.;
	
}


/*	A valtoztathato parameterek beolvasa terminalbol	*/
int parse_options(int argc, const char **argv, options_t *opt){
	int i = 1;
	double temp = 1;

	while (i < argc) {
		char *p = malloc(sizeof(argv[i]));	/*	Ahhoz, hogy "szoveget" (char, mert a C nem tud stringet kezelni) tudjunk beolvasni, es ossze tudjuk vetni a lentebb megadott kapcsolokkal (pl. -drift), le kell foglalni a memoriateruletet p-nek. Ennek viszont elore nem tudjuk a meretet, ezert dinamikusan foglaljuk le azt	*/
   		strcpy(p, argv[i]);			/*	A p-be elmentjuk az adott argumentumot az strcpy paranccsal	*/	

		if(strcmp(p, "-drift") == 0) {		/*	Az strcmp parancs kepes osszevetni a karaktersorozatunkat egy masikkal: ha azok megegyeznek, akkor a kovetkezo beolvasott argumentumot elmenti a megadott struktura elembe	*/
			strcpy(p, argv[i]);
			i++;
			opt->drift = atof(argv[i]);
		}
		else if (strcmp(p, "-growth") == 0) {
			i++;
			opt->growth = atof(argv[i]);
		}
		else if (strcmp(p, "-evol") == 0) {
			i++;
			opt->evol = atof(argv[i]);
		}
		else if (strcmp(p, "-twopop") == 0) {
			i++;
			opt->twopop = atof(argv[i]);
		}
		else if (strcmp(p, "-tStep") == 0) {
			i++;
			opt->tStep = atof(argv[i]);
		}
		else if (strcmp(p, "-n") == 0) {
			i++;
			opt->ngrid = atoi(argv[i]);
		}
		else if (strcmp(p, "-i") == 0) {
			i++;
			opt->input = argv[i];
			temp = 0;
		}
		else {				/*	ha nem sikerul a parancssorbol a beolvasas, akkor a program kilep es feldobja, hogy mely kapcsolo nem volt jo, helyette mit probaljon meg	*/
			printf("\n\n**** Invalid switch on command-line: %s! ****\n\n\n",p);
			printf("**** Try following parameters: ****\n\n-drift\nif set to 0: dust drift is not included, if set to 1: dust drift is included (default)\n\n-growth\nif set to 0: partcle growth is not included, if set to 1: particle growth is included (default)\n\n-evol\nif set to 0: the evolution of the surface density is not included, if set to 1: evolution of the surface density is included (default)\n\n-n\nnumber of grid cells (2000 by default)\n");
			return 1;
		}
		i++;
		free(p);			/*	memoria foglalas miatt fel is kell szabaditanunk a memoria teruletet	*/
	}

	optinp = temp;
	
	return 0;
}


/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int reszecskek_szama(int numout, const char *filenev){

	char c;
	fin1 = fopen(filenev,"r+");		
	numout = 0;

/*	A porreszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal				*/

		while((c = fgetc(fin1)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 							*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert, ezzel beolvastuk, hogy hany soros a file				*/
	fclose(fin1);	

	return numout;

}


/*	Az idot tartalmazo file parametereinek beolvasasa	*/
void sigIn(double *sigvec, double *rvec) {

	double sig,r;
	int i;

	FILE *densin;
	densin = fopen(inputsig,"r");
	
	for(i = 0; i < NGRID; i++) {
           	if(fscanf(densin,"%lg  %lg",&r,&sig) == 2) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 3, mert 3 oszlopot szeretnénk beolvasni.	*/
			rvec[i+1] = r;				/*	r vektor	*/
			sigvec[i+1] = sig; 			/*	adott lepeskozonkent irja ki a file-t, megvan, hogy hany evente, tudjuk, hogy meddig fusson a kod, igy egy egyszeru osztassal meg lehet adni, hogy az mindig hanyadik idolepes	*/
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			exit(EXIT_FAILURE);
   	        }
	}

	RMIN = rvec[1];
	RMAX = rvec[NGRID];

	fclose(densin);
	
}


/*	Az idot tartalmazo file parametereinek beolvasasa	*/
void timePar(double *tMax, double *step, double *current) {

	double tmax,stepping,curr;

	fin2 = fopen(filenev3,"r");

           	if(fscanf(fin2,"%lg  %lg %lg",&tmax,&stepping,&curr) == 3) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 3, mert 3 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  tmax: %lg, a program %lg evenkent irja ki a file-okat\n\n\n",tmax,stepping);  
			*tMax = tmax;				/*	meddig fusson a program	*/
			stepping = tmax/stepping;
			*step = stepping; 			/*	adott lepeskozonkent irja ki a file-t, megvan, hogy hany evente, tudjuk, hogy meddig fusson a kod, igy egy egyszeru osztassal meg lehet adni, hogy az mindig hanyadik idolepes	*/
			*current = curr;			/*	beolvassa, hogy mennyi volt az ido eppen a futas inditasanak pillanataban -- ez a kiiratasnal lenyeges	*/
			
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

           	if(fscanf(fin2,"%lg  %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&rmin,&rmax,&dummy2,&exp,&sig0,&dummy,&rdzei,&rdzeo,&drdzei,&drdzeo,&amod,&rhop,&alpha,&mstar,&flind) == 15) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 15, mert 15 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  sigma0: %lg, sdexp: %lg\n\n\n", sig0,exp);  
			*sigma0 = sig0;						/*	sigma az r=1AU tavolsagon		*/
			*sdexp = exp; 						/*	a feluleti suruseg profiljanak kitevoje	*/
			*Rmin = rmin;						/*	a korong belso hatara			*/
			*Rmax = rmax;       					/*	a korong kulso hatara			*/
			*r_dzei = rdzei;					/*	a DZE belso hatara			*/
			*r_dzeo = rdzeo;					/*	a DZE kulso harara			*/
			*dr_dzei = drdzei;					/*	a DZE belso hatar atmenet vastagsaga	*/
			*dr_dzeo = drdzeo;					/*	a DZE kulso hatar atmenet vastagsaga	*/
			*alph_mod = amod;					/*	a viszkozitas redukcio merteke		*/
			*rho_p = rhop;						/*	a reszecskek atlagos surusege cgs-ben	*/
			*rho_p_dimless = rhop / SUN2GR * AU2CM * AU2CM * AU2CM;	/*	a reszecskek atlagos surusege (dimtlan)	*/
			*alphav = alpha;					/*	alpha viszkozitas merteke		*/
			*mStar = mstar;						/*	a kozponti csillag tomege		*/
			*gamma = flind;						/*	flaring index				*/

	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
//			getchar();
			exit(EXIT_FAILURE);
   	        }

	fclose(fin2);
	
}


/*	A porreszecskek adatainak beolvasasa	*/
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr) {

	int i, dummy;
	double distance, particle_radius, reprmass, reprmassmicr,radmicr;
	
   	fin1 = fopen(filenev1,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%d %lg %lg %lg %lg %lg",&dummy,&distance,&reprmass,&reprmassmicr,&particle_radius,&radmicr) == 6) {	
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 6, mert 6 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/

/*	A cm-es porreszecskek adatainak beolvasasa!		*/
           		radius[i][0] = distance;			/*	a reszecske tavolsaga AU-ban		*/
	   		radius[i][1] = particle_radius / AU2CM;		/* 	a részecske mérete AU-ban		*/
			mass[i] = reprmass;				/*	a porreszecske altal kepviselt reprezentativ tomeg dimenziotlan egysegekben					*/


/*	A mikronos reszecskek adatainak beolvasasa!		*/
            		radiusmicr[i][0] = distance;			/*	a mikronos reszecske tavolsaga AU-ban --> kezdetben ugyanolyan messze van az osszes, mint a cm-es reszecskek!!	*/
	   		radiusmicr[i][1] = radmicr / AU2CM;		/* 	a mikronos részecske mérete AU-ban	*/
			massmicr[i] = reprmassmicr;			/*	a porreszecske altal kepviselt reprezentativ tomeg dimenziotlan egysegekben					*/

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


/*	a sigmara kezdeti profil betoltese	*/
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


/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt) {
     
     	double P, H, dPdr, St, csound;
    
	St = Stokes_Number(pradius,sigma);
	H = scale_height(r);  
	P = press(sigma,r);
	dPdr = dp;
	csound = c_sound(r); 

     	*drdt = ug / (1. + St * St) + St / (1. + St * St) * H / P * dPdr * csound;	/* 	bearamlas sebessege: Birnstiel PHD	*/
     
}


/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt) {

	double rmid, rindex, coef1, temp;
	int index; 

     	rmid = pos - RMIN;
	rmid = rmid / rd;     						/* 	the integer part of this gives at which index is the body	*/
	index = (int) floor(rmid);					/* 	ez az rmid egesz resze	(kerekites 0.5-tol)			*/
	rindex = rvec[index];       					/* 	the corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

 	coef1 = (invec[index + 1] - invec[index]) / rd; 		/*	ez az alabbi ket sor a linearis interpolacio - remelem, jo!!!	*/
	temp = invec[index] + coef1 * (pos - rindex);          		/*	a beerkezo dimenzionak megfelelo mertekegysegben		*/

	if(opt == 1) if(temp < 0) temp = -1.*temp;

	*out = temp;

}

/*			BIRNSTIEL EL AL 2012			*/

//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret
/*double a_drift() {

	s_drift = f_drift * 2.0 / M_PI * Sigmad_cgs/rho_p * v_kep2 / c_s2 * fabs(1.0 / dlnPdlnr);

}
*/


// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
double a_turb(double sigma, double r, double rho_p) {

	double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

	u_frag = uFrag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
  	u_frag2 = u_frag * u_frag;
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r); // / CMPSECTOAUPYRP2PI;
	c_s2 = c_s * c_s;

	s_frag = fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(r)) * u_frag2 / c_s2;
	
	return s_frag;

}


// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
double a_df(double sigma, double r, double p, double dp, double rho_p) {

	double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

	u_frag = uFrag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r);
	c_s2 = c_s * c_s;
	dlnPdlnr = r / p * dp;
	vkep = v_kep(r);

	s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

	return s_df;

}


/*	a reszecskek novekedesenek idoskalaja	*/
double tauGr(double r, double eps) {

	double omega = kep_freq(r);
	double taugr = eps / omega;
	
	return taugr;

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


/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress, double dt) {

	double sturb = a_turb(sigma,y,pdens);			// cm-ben
	double sdf = a_df(sigma,y,p,dpress,pdens);		// cm-ben
	double smin = find_min(sturb,sdf,HUGE_VAL);		// cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja	
//	double eps = sigma / 100.;
	double eps = sigma / sigmad;
	double tau_gr = tauGr(y,eps);
	double rt;
	
	smin = smin / AU2CM;					// AU-ban

/*	kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
	if (prad < smin) rt = find_min(prad*exp(dt/tau_gr),smin,HUGE_VAL);
	if (prad >= smin) rt = smin;

	return rt;

}


/*	Runge-Kutta4 integrator	*/
// prad bemenet: AU-ban!
void int_step(double time, double prad, double *pressvec, double *dpressvec, double *sigmavec, double sigmad, double *rvec, double *ugvec, double step, double y, double *ynew, double *pradnew) {
	double dy1,dy2,dy3,dy4;
	double ytemp;
	double sigma, dpress, ugas; 
	double pdens, p;
	double pradtemp;

	int opt = 0;
	
/*	Mivel a kulongozo parametereket csak a megadott gridcella pontokban ismerjuk, de ez nem feltetlen egyezik meg a reszecskek pozicijaval, ezert minden fontos parametert interpolalunk a reszecskek tavolsagara	*/
	interpol(sigmavec,rvec,y,&sigma,DD,opt);
	interpol(dpressvec,rvec,y,&dpress,DD,opt);
	interpol(ugvec,rvec,y,&ugas,DD,opt);

	if(optgr == 1.) {		// ha van reszecskenovekedes

		if(time != 0.) {	// ha nem t0 idopontban vagyunk
	
			pradtemp = prad;
			interpol(pressvec,rvec,y,&p,DD,opt);
			pdens = PDENSITY; 
			pradtemp = getSize(prad,pdens,sigma,sigmad,y,p,dpress,step);	// itt szamolja a reszecskenovekedest
			prad = pradtemp;

		}
	
	}

	*pradnew = prad;

/*	Itt szamolja a reszecske poziciojat	*/
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
double find_max(double r[][2], int n) {

	int i;
	double maxim = -1.;

	for(i = 0; i < n; i++) {

		if (r[i][0] > maxim) {
			maxim = r[i][0];
		}
	
	}

	return maxim;
}


/*	counting the number of zero points of the pressure gradient function	*/
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


/*	mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke	*/
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
	
	if(((dp[i] * dp[i+1]) <= 0.) && (dp[i] > dp[i+1])) {		/*	Ha a ket pont szorzata negativ --> elojel valtas a dp-ben, nyomasi min/max. Maximum hely ott van, ahol pozitivbol negativba valt az ertek	*/
		r = find_r_zero(rvec[i],rvec[i+1],dp[i],dp[i+1]);	/*	Ha elojel valtas tortenik es nyomasi maximum van, akkor kiszamolja a ket pont kozott, hogy hol lenne a zerus hely pontosan	*/
	} else {
		r = 0.0;
	}

	return r;

}



/*	A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(double *rvec, double rin, double *ind_ii, double *ind_io, double rout, double *ind_oi, double *ind_oo) {

	int i;
	double rmid, rtemp;
	double roimH;
	double roipH;
	double roomH;
	double roopH;
	double riimH;
	double riipH;
	double riomH;
	double riopH;

	if(optdze == 0) {
	
		*ind_ii = 0;
		*ind_io = 0;
	
	}

	riimH = (rin - scale_height(rin)) - DD / 2.0;		/*	A nyomasi maximum az rout pontban van, ettol rout - 1/2H - DD / 2 es rout + 1*2H -DD / 2 kozott van a korgyuru belso hatara (azert DD/2, hogy biztosan 1 cellat tudjunk kijelolni, ne pedig egy tartomanyt)	*/
	riipH = (rin - scale_height(rin)) + DD / 2.0;		
	riomH = (rin + scale_height(rin)) - DD / 2.0;		/*	Az alabbi ketto pedig a kulso hatarat adja meg a korgyurunek	*/
	riopH = (rin + scale_height(rin)) + DD / 2.0;

	roimH = (rout - scale_height(rout)) - DD / 2.0;		/*	A nyomasi maximum az rout pontban van, ettol rout - 1/2H - DD / 2 es rout + 1*2H -DD / 2 kozott van a korgyuru belso hatara (azert DD/2, hogy biztosan 1 cellat tudjunk kijelolni, ne pedig egy tartomanyt)	*/
	roipH = (rout - scale_height(rout)) + DD / 2.0;		
	roomH = (rout + scale_height(rout)) - DD / 2.0;		/*	Az alabbi ketto pedig a kulso hatarat adja meg a korgyurunek	*/
	roopH = (rout + scale_height(rout)) + DD / 2.0;

	for(i = 1; i <= NGRID; i++) {

		if(optdze == 1) { 
/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
			if(rvec[i] > riimH && rvec[i] < riipH) {
			    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
				rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
				*ind_ii = rtemp;
			}
				
/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
			if(rvec[i] > riomH && rvec[i] < riopH) {
			    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
				rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
				*ind_io = rtemp;
			}
		}


/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
		if(rvec[i] > roimH && rvec[i] < roipH) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oi = rtemp;
		}
				
/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
		if(rvec[i] > roomH && rvec[i] < roopH) {
		    	rmid = (rvec[i] - RMIN)/ DD;     						/* 	The integer part of this gives at which index is the body	*/
			rtemp = (int) floor(rmid + 0.5);						/*	Rounding up (>.5) or down (<.5)					*/
			*ind_oo = rtemp;
		}

		if(rvec[i] > roopH) break;

	}

}


void GetMass(int n, double partmassind[][4],int indii, int indio, double tavi, double dzei, double *massiout,int indoi, int indoo, double tavo, double dzeo, double *massoout) {

	int i;
	double index_i;
	double massitemp = 0., massotemp = 0;

	if(optdze == 1) {
		for (i = 0; i < n; i++) {

/*	A partmassind tombben van eltarolva a reszecske indexe (ez az integralas soran allando, ezzel lehet abrazolni szinesben), az altala kepviselt reprezentativ tomeg (ez az integralas soran allando) es a grid cella indexe, amelyben eppen tartozkodik (ezt minden lepesben kiszamolja a program	*/
			index_i = partmassind[i][1];			/*	a gridcella indexe, amelyben a reszecske tartozkodik eppen	*/

/*	Ha a reszecse a kiszamolt korgyuruben tartozkodik, akkor a mass valtozohoz hozzaadjuk a tomeget --> igy tudjuk a tomegnovekedest kiszamolni	*/

			if(tavi != dzei) {

				if ((index_i >= (int)indii) && (index_i <= (int)indio) & partmassind[i][3] == 0) {
					partmassind[i][3] = 1;
					massitemp = massitemp + partmassind[i][0];
				}
			
			} else {
	
				if ((index_i >= (int)indii) && (index_i <= (int)indio)) {
					massitemp = partmassind[i][0] + massitemp;
				}
			}

			if(tavo != dzeo) {

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) & partmassind[i][3] == 0) {
					partmassind[i][3] = 1;
					massotemp = massotemp + partmassind[i][0];
				}

			} else {
	
				if ((index_i >= (int)indoi) && (index_i <= (int)indoo)) {
					massotemp = partmassind[i][0] + massotemp;
				}
			}

		}
	} else {

		for (i = 0; i < n; i++) {

/*	A partmassind tombben van eltarolva a reszecske indexe (ez az integralas soran allando, ezzel lehet abrazolni szinesben), az altala kepviselt reprezentativ tomeg (ez az integralas soran allando) es a grid cella indexe, amelyben eppen tartozkodik (ezt minden lepesben kiszamolja a program	*/
			index_i = partmassind[i][1];			/*	a gridcella indexe, amelyben a reszecske tartozkodik eppen	*/

/*	Ha a reszecse a kiszamolt korgyuruben tartozkodik, akkor a mass valtozohoz hozzaadjuk a tomeget --> igy tudjuk a tomegnovekedest kiszamolni	*/
		
			if(tavo != dzeo) {

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) & partmassind[i][3] == 0) {
					partmassind[i][3] = 1;
					massotemp = massotemp + partmassind[i][0];
				}

			} else {
	
				if ((index_i >= (int)indoi) && (index_i <= (int)indoo)) {
					massotemp = partmassind[i][0] + massotemp;
				}
			}
		}

	}

	*massiout = massitemp;
	*massoout = massotemp;
}



/*	Fuggveny a tomegfile kiiratasara	*/
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double t, double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout) {

	double ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;	

	tav = r_dze_o;	
	tav2 = r_dze_i;	

	int dim = find_num_zero(rvec,dpressvec);		// megnezi, hogy hany, nyomasi maximumbol szarmazo, 0 pont van a derivaltban
	double r_count[dim];					// a nullpontnak megfelelo elemu tombot letrehozza (pl. ha van dze_i es dze_o is, akkor kulon igy lehet elmenteni azok helyet)
	double temp_new = 0.;
	double temp = 0.;
	double rout = r_dze_o;
	double rin = r_dze_i;
	double rin_new, rout_new;

	int j, i;
	j=0;
	
	if(dim != 0) {						// ha van nullpont, akkor megkeresi, hogy hol
		for(i = 0; i < NGRID; i++) {
			temp_new = find_zero(i,rvec,dpressvec);	// ha van nullpont, akkor a temp_new valtozoba tarolja el --> ehhez vegig megy az egesz r-en 

			if(temp != temp_new && i > 3 && temp_new != 0.0) {
				r_count[j] = temp_new;		// a temp_new-t itt tarolja el, azaz a nullpontok szamanak megfeleloen itt tartolodnak el a nyomasi maximumok
				j++;
			}

/*	Ha csak kulso dze van, akkor a dze uj helyet itt tarolja el	*/
			if(optdze == 0) {
				if(temp_new > 0.) {			
					temp = temp_new;
					rout_new = temp;
				} 
			}
		}
	}


/*	Ha van belso dze is, akkor itt menti el a kuslo es belso dze helyet	*/		
	if(optdze == 1) {
		if(dim > 0) {
			if (dim == 1) {
				rin_new = r_count[0];
				rout_new = rout;
			} else {
				rin_new = r_count[0];
				rout_new = r_count[1];
			}
		} 
		if(dim == 0) {	// ha nincs nyomasi maximum meg, akkor a regi valtozok erteket (azaz a nyomasi dze_i es dze_o helyeket) menti el
			rin_new = rin;
			rout_new = rout;
		}
	}

	rin = rin_new;
	if(optdze == 0) rin = 0;	// ha nincs belso dze, akkor annak a helye 0 (ez vegulis ebben az esetben nem lenyeges)
	rout = rout_new;
	tav2 = rin;
	tav = rout;

/* 	MEG KELL OLDANI, HOGY AKKOR IS TUDJON TOMEGNOVEKEDEST SZAMOLNI, HA CSAK BELSO DZE VAN!	*/

	find_r_annulus(rvec,tav2,&ind_ii,&ind_io,tav,&ind_oi,&ind_oo);		/*	A belso es kuslo nyomasi maximum korul 2H tavolsagban keres korgyurut, a fuggveny visszaadja a cellak indexet	*/
	
	double masst0i = 0, massii = 0, massoi = 0;
	double masst0im = 0, massiim = 0,massoim = 0;
	double massis = 0, massos = 0;

	GetMass(PARTICLE_NUMBER,partmassind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massii,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massoi);
	
	if(opttwopop == 1) {
		GetMass(4*PARTICLE_NUMBER,partmasssecind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massis,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massos);
		GetMass(PARTICLE_NUMBER,partmassmicrind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massiim,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massoim);
	} 

	double massi, massim, masso, massom;

	if(tav2 != r_dze_i) {
		massi = massii + massbtempii + massis;
		massim = massiim + massmtempii;
	} else {
		massi = massii + massis;
		massim = massiim;
	}
	if(tav != r_dze_o) {
		masso = massoi + massbtempoi + massos;
		massom = massoim + massmtempoi;
	} else {
		masso = massoi + massos;
		massom = massoim;
	}

	*massbtempio = massi;
	*massbtempoo = masso;
	*massmtempio = massim;
	*massmtempoo = massom;

	*tavin = tav2;
	*tavout = tav;

	fprintf(massfil,"%lg %lg %lg %lg %lg\n",step,tav2,massi+massim,tav,masso+massom);
	fflush(massfil);

}



/*	Fuggveny a sigma, p, dp kiiratasara	*/
void Print_Sigma(char *dens_name, double *rvec, double *sigmavec, double *pressvec, double *dpressvec) {

	int i;
	fmo = fopen(dens_name,"w");				

 	for(i = 1; i <= NGRID; i++) {
   		fprintf(fmo, "%lg   %lg   %lg   %lg\n", rvec[i],sigmavec[i],pressvec[i],dpressvec[i]);
	}

	fclose(fmo);

}


/*	Fuggveny a por feluletisurusegenek kiiratasara	*/
void Print_Sigmad(char *dust_name, char *dust_name2, double min, double *r, double *rm, double *sigmad, double *sigmadm) {

	int i;
	double rtempvec[PARTICLE_NUMBER];
	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1);	/*	Mivel a jelen futas gridfelbontasa nem feltetlen egyezik meg a porreszecskeket generalo program gridfelbontasaval - ez a feluletisuruseg miatt lenyeges! - ezert itt szamolja ki a program	*/

	FILE *sid;	

	fil = fopen(dust_name,"w");
	if(opttwopop == 1) sid = fopen(dust_name2,"w");		/*	Ha 2pop --> megnyit egy kulon file-t a mikronos por feluletisurusegenek kiiratasara	*/

	for(i=0;i<PARTICLE_NUMBER;i++){

		if (r[i] >= RMIN) {			/*	a cm-es por feluletisurusege	*/
			fprintf(fil,"%lg  %lg \n",r[i],sigmad[i]);
		} 

		if(opttwopop == 1) {				/*	itt irja ki a mikronos por feluletisuruseget	*/

			if (rm[i] >= RMIN) {
				fprintf(sid,"%lg  %lg \n",rm[i],sigmadm[i]);
			}
		}
	}

	fclose(fil);
	if(opttwopop == 1) fclose(sid);

}


/*	Fuggveny a pormozgas kiiratasara	*/
void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2], double *rvec, double t){

	int i;
	if(optgr == 1) fout3 = fopen(size_name,"w");	/*	ha van pornovekedes	*/

	for(i=0;i<PARTICLE_NUMBER;i++){

/*	pormozgas.dat kiiratasa --> adott idokozonkent a cm-es reszecskek tavolsaga, indexe es az ido kiiaratasa egy file-ba	*/
		if(rad[i][0] >= RMIN) fprintf(fout,"%lg %d %lg\n",(double)step,i,rad[i][0]);
/*	ha a szimulacio 2 populacios, akkor a fentihez hasonlo file letrehozasa a mikronos reszecskekre	*/		
		if(opttwopop == 1.) if(radmicr[i][0] >= RMIN) fprintf(foutmicr,"%lg %d %lg\n",(double)step,i,radmicr[i][0]);
/*	ha van reszecske novekedes, akkor letrehoz egy file-t minden adott idolepesben, ami a centimeteres porreszecskek meretet tartalmazza	*/
		if(optgr == 1) {
			if (rad[i][0] >= RMIN) fprintf(fout3,"%lg  %lg  %lg \n",(double)step, rad[i][0], rad[i][1]*AU2CM);
		}
	}

	fflush(fout);
	if(opttwopop == 1.) fflush(foutmicr);
	if(optgr == 1) fclose(fout3);

}


/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sort(double *rv,int n) {

	double temp, temp2;
	int i, step;

	for(step = 1; step <= (n-1); step++) {

		for(i = 0; i <= (n-2); i++) {

			if(rv[i] > rv[i + 1]) {

				temp = rv[i];
				rv[i] = rv[i + 1];
				rv[i + 1] = temp;

			}
		}
	}
}


void histogram(double r, int *hist, double dd) {

	int index;
	double rmid, hist_i;
    	rmid = (r - RMIN) / dd;     						/* 	The integer part of this gives at which index is the body	*/
   	index = (int) floor(rmid);					/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
	hist_i = hist[index];
    	hist[index] = hist_i + 1;       					/* 	The corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

}


void loadSigDust(double radin[][2], double *massin, double out[][3], double dd, int n) {

	int i;

	for(i=0;i<n;i++){

/*	cm-es por feluletisurusegenek kiszamolasa	*/
/*	ha a reszecske tavolsaga nagyobb, mint RMIN, azaz a szamolas tartomanyan belul van, a feluletisuruseget az altala kepviselt tomegbol szamolja vissza a program	*/
		if((radin[i][0] >= RMIN)) {
			out[i][0] = massin[i] / (2. * (radin[i][0]-dd/2.) * M_PI * dd);	// sigma = m /(2 * r * pi * dr) --> itt a dr az a tavolsag, ami porreszecske "generalo" programban az eredeti gridfelbontas
			out[i][1] = radin[i][0];					// elmenti a reszecske tavolsagat is

  			double rmid = (radin[i][0] - RMIN) / dd;     						/* 	The integer part of this gives at which index is the body			*/
			int rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
			out[i][2] = (double) rindex;

/*	ha a reszecske RMIN-en belul van, akkor az o "tavolsagaban" a feluletisuruseg 0		*/	
		} else {
			out[i][0] = 0;
			out[i][1] = 0;					// r = 0, mert ha RMIN-en belulre kerult a reszecske, akkor a program automatikusan kinullazza a reszecske tavolsagat. Itt tehat a sigdtemp[i][1] = 0 lesz!
			out[i][2] = 0;
		}
	}
}


void contract(double in[][3], double dd, int n) {

	int i;
	int j;
	int k;
	double sig = 0, radout[n], sigout[n];

	for(i = 0; i < n; i++) {
		radout[i] = 0;
		sigout[i] = 0;
		in[i][2] = 0;
	}

	i = 0;
	j = 0;
	k = 0;

	do {
		if(in[i][1] != in[i+1][1]) {
			radout[j] = in[i][1];
			sigout[j] = in[i][0];
			sig = 0;
			k = 0;
			j++;
			i++;
		} else {
			do {
				radout[j] = in[i][1];
				sig = sig + in[i+k][0];
				sigout[j] = sig;
				k++;
			} while (in[i][1] == in[i+k][1]);
			i = i+k;
			k =0;
			j++;
		}

	} while (i < n);

	for(i = 0; i < n; i++) {
	
		in[i][0] = sigout[i];
		in[i][1] = radout[i];
  		double rmid = (in[i][1] - RMIN) / dd;     						/* 	The integer part of this gives at which index is the body			*/
		int rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
		in[i][2] = (double)rindex;
	}

}


/*	Fuggveny a por feluletisurusegenek kiszamolasara	*/
void Get_Sigmad(double L, double max, double min, double rad[][2], double radmicr[][2], double radsec[][2], double *sigma_d, double *sigma_dm, double *sigma_ds, double *massvec, double *massmicrvec, double *masssecvec, double *rd, double *rmic, double *rs) {

	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1);	/*	Mivel a jelen futas gridfelbontasa nem feltetlen egyezik meg a porreszecskeket generalo program gridfelbontasaval - ez a feluletisuruseg miatt lenyeges! - ezert itt szamolja ki a program	*/
	int i,j,k;
	double sigdtemp[PARTICLE_NUMBER][3], sigdmicrtemp[PARTICLE_NUMBER][3], sigdsectemp[4*PARTICLE_NUMBER][3], rtempvec[PARTICLE_NUMBER][2], sigtempvec[PARTICLE_NUMBER];
	double radtemp[PARTICLE_NUMBER], sdtemp[PARTICLE_NUMBER], radmicrtemp[PARTICLE_NUMBER], sdmicrtemp[PARTICLE_NUMBER], radsectemp[4*PARTICLE_NUMBER], sdsectemp[4*PARTICLE_NUMBER];
	double intsig[PARTICLE_NUMBER], intsigmicr[PARTICLE_NUMBER];
	double rtemp[PARTICLE_NUMBER];
	int rtempi[PARTICLE_NUMBER];

	for(i=0;i<PARTICLE_NUMBER;i++){

		rtempvec[i][0] = RMIN + i * dd;
		rtempvec[i][0] = rtempvec[i][0] + dd / 2.0;
		rtempvec[i][1] = i;
		radtemp[i] = 0;
		radmicrtemp[i] = 0;
		sdtemp[i] = 0;
		sdmicrtemp[i] = 0;
		intsig[i] = 0;
		intsigmicr[i] = 0;
		rd[i] = 0;
		rmic[i] = 0;
		rtemp[i] = 0;
		rtempi[i] = 0;
		sigtempvec[i] = 0;
	}

	for(i=0;i<4*PARTICLE_NUMBER;i++){
		radsectemp[i] = 0;
		sdsectemp[i] = 0;
		sigdsectemp[i][0] = 0;
		sigdsectemp[i][1] = 0;
		sigdsectemp[i][2] = 0;
		rs[i] = 0;
	}

	loadSigDust(rad,massvec,sigdtemp,dd,PARTICLE_NUMBER);
	if(opttwopop == 1) {
		loadSigDust(radmicr,massmicrvec,sigdmicrtemp,dd,PARTICLE_NUMBER);
		loadSigDust(radsec,masssecvec,sigdsectemp,dd,4*PARTICLE_NUMBER);
	}

	contract(sigdtemp,dd,PARTICLE_NUMBER);
	if(opttwopop == 1) {
		contract(sigdmicrtemp,dd,PARTICLE_NUMBER);
		contract(sigdsectemp,dd,4*PARTICLE_NUMBER);
	}

	
/*	a feluletisurusegek interpolalasa 2 populacio eseten	*/
	if(opttwopop == 1) {
		for(i=0;i<PARTICLE_NUMBER;i++) {
			for(j=0; j < PARTICLE_NUMBER; j++) {
				if((int)rtempvec[i][1] <= (int)sigdtemp[j][2] && (int)sigdtemp[j][2] < (int)rtempvec[i+1][1]) {
					intsig[i] = intsig[i] + sigdtemp[j][0];
				}
				if((int)rtempvec[i][1] <= (int)sigdmicrtemp[j][2] && (int)sigdmicrtemp[j][2] < (int)rtempvec[i+1][1]) {
					intsigmicr[i] = intsigmicr[i] + sigdmicrtemp[j][0];
				}
				if(((int)rtempvec[i][1] <= (int)sigdsectemp[j][2] && (int)sigdsectemp[j][2] < (int)rtempvec[i+1][1]) && sigdsectemp[j][1] != 0) {
					intsig[i] = intsig[i] + sigdsectemp[j][0];
				}
				if(((int)rtempvec[i][1] <= (int)sigdsectemp[j+PARTICLE_NUMBER][2] && (int)sigdsectemp[j+PARTICLE_NUMBER][2] < (int)rtempvec[i+1][1]) && sigdsectemp[j][1] != 0) {
					intsig[i] = intsig[i] + sigdsectemp[j+PARTICLE_NUMBER][0];
				}
				if(((int)rtempvec[i][1] <= (int)sigdsectemp[j+2*PARTICLE_NUMBER][2] && (int)sigdsectemp[j+2*PARTICLE_NUMBER][2] < (int)rtempvec[i+1][1]) && sigdsectemp[j][1] != 0) {
					intsig[i] = intsig[i] + sigdsectemp[j+2*PARTICLE_NUMBER][0];
				}
				if(((int)rtempvec[i][1] <= (int)sigdsectemp[j+3*PARTICLE_NUMBER][2] && (int)sigdsectemp[j+3*PARTICLE_NUMBER][2] < (int)rtempvec[i+1][1]) && sigdsectemp[j][1] != 0) {
					intsig[i] = intsig[i] + sigdsectemp[j+3*PARTICLE_NUMBER][0];
				}
			}
		}
	}

	if(opttwopop == 1) {
		j = 0;
		for(i=0; i < PARTICLE_NUMBER; i++){
			sigtempvec[j] = intsig[i] + intsigmicr[i];
			rtemp[j] = rtempvec[i][0];
			rtempi[j] = rtempvec[i][1];
			j++;
		}
	}

	if(opttwopop == 1) {

		double y, dust;
		double step1, step2;
		int k = -1;
		int l = -1;
	
		for(i = 0; i < PARTICLE_NUMBER; i++) {
			double simt = 0;
			double r = 0;
			for(j=0; j<PARTICLE_NUMBER;j++) {
	
				if((int)sigdtemp[j][2] >= rtempi[i] && (int)sigdtemp[j][2] < rtempi[i+1]){
					double y = sigdtemp[j][1];
					double rd = rtemp[i+1]-rtemp[i];
					interpol(sigtempvec,rtemp,y,&simt,dd,opttwopop);
					r = sigdtemp[j][1];
					k++;
				} 
			}

			if(simt != 0 && r != 0) {
				sigma_d[k] = simt;	/*	ide rakja az interpolalt feluletisuruseget	*/
				rd[k] = r;
			}		
				
			simt = 0;
			r = 0;
			for(j=0; j<PARTICLE_NUMBER;j++) {
	
/*	Az elozohoz hasonloan interpolalja a mikronos feluletisuruseget a cm-es porszemcsek tavolsagara		*/	
				if((int)sigdmicrtemp[j][2] >= rtempi[i] && (int)sigdmicrtemp[j][2] < rtempi[i+1]){
					double y = sigdmicrtemp[j][1];
					double rd = rtemp[i+1]-rtemp[i];
					interpol(sigtempvec,rtemp,y,&simt,dd,opttwopop);
					r = sigdmicrtemp[j][1];
					l++;
				} 
			}

			if(simt != 0) {
				sigma_dm[l] = simt;	/*	ide rakja az interpolalt feluletisuruseget	*/
				rmic[l] = r;
			}
		}

		k = -1;
		l = -1;

		for(i = 0; i < PARTICLE_NUMBER-1; i++) {
			double simt = 0;
			double r = 0;
			for(j=0; j<4*PARTICLE_NUMBER;j++) {
	
				if((int)sigdsectemp[j][2] >= rtempi[i] && (int)sigdsectemp[j][2] < rtempi[i+1]){
					double y = sigdsectemp[j][1];
					double rd = rtemp[i+1]-rtemp[i];
					interpol(sigtempvec,rtemp,y,&simt,dd,opttwopop);
					r = sigdsectemp[j][1];
					k++;
				} 
			}

			if(simt != 0 && r != 0) {
				sigma_ds[k] = simt;	/*	ide rakja az interpolalt feluletisuruseget	*/
				rs[k] = r;
			}
		}

	} else {

/*	ha nem 2 pop a futas, akkor nincs interpolacio, szimplan a lenti tomb feltoltes	*/
		for(i=0; i < PARTICLE_NUMBER; i++) {
			rd[i] = sigdtemp[i][1];
			rmic[i] = 0;
			sigma_d[i] = sigdtemp[i][0];
			sigma_dm[i] = 0;
		}

		for(i=0; i < 4*PARTICLE_NUMBER; i++) {
			rs[i] = 0;
			sigma_ds[i] = 0;
		}
	}
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


/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(char *nev, int opt, double radius[][2], double *pressvec, double *dpressvec, double *sigmavec, double *sigmad, double *rvec, double *ugvec, double deltat, double t, int n) {		/*	a bemenet dimenziótlan egységekben van!!!	*/

	int i;
	double y, y_out, prad_new, particle_radius;
	char scout[1024];

/*	A futas kezdeten kiirja egy file-ba a drift idoskalajat!	*/
	if(t==0) {
		sprintf(scout,"%s/timescale.dat",nev);
		fout2 = fopen(scout,"w");
	}

	for (i = 0; i < n; i++) {	

		double drdt[n];

		if (radius[i][0] > RMIN && radius[i][0] < RMAX) { 

			y = radius[i][0];					// a reszecske tavolsaga
     			particle_radius = radius[i][1]; 			// a reszecske merete	

/*	a reszecske tavolsaganak kiszamolasa	*/	
     			int_step(t,particle_radius,pressvec,dpressvec,sigmavec,sigmad[i],rvec,ugvec,deltat,y,&y_out,&prad_new);

/*	itt az opt a 2 populacios modellre utal, ha az 1, akkor a mikronos reszecskekre szamol - ehelyett egyebkent valoszinu lehetne az opttwopop-ot is hasznalni	*/
			if(t == 0) {
				if(opt == 0) {
					drdt[i] = (fabs(y_out - y)/(deltat));	/*	csak a cm-es esetben irja ki a drift idoskalajat	*/
						
/*	timescale.dat, ez tartalmazza a reszecskek bearamlasanak idoskalajat (evben)				*/
/*	A timescale.dat-hoz a t_drift = r / v(r,peak) kepletet a kovetkezo modon oldja meg a program:		*/
/*	Kiszamolja, hogy a kovetkezo idopillanatban hol lenne minden reszecske, majd az uj koordinatabol (x(1))	*/
/*	a regit (x(0)) levonva, majd az egeszet elosztva dt-vel megkapjuk a bearamlas sebesseget (drdt)		*/
/*	ez utan kiszamolhato az eredeti t_drift keplet, ezt iratom ki a file-ba evben (ezert osztok le 2PI-vel	*/

					fprintf(fout2,"%lg %lg\n",radius[i][0], (radius[i][0] / drdt[i])/2.0/M_PI);
				}

			}
	
/*	cm-es porokra az uj meret es pozicio elmentese				*/
			if (opt != 1){
				radius[i][1] = prad_new;
		     		radius[i][0] = y_out;
			} 

/*	mikronos esetben a az uj pozicio elmentese, azonban nincs novekedes	*/	
			if (opt == 1){
		     		radius[i][0] = y_out;
			}

/*	ha a reszecske RMIN-en belulre kerul, akkor a tavolsaga 0 lesz	*/		
		} else {
			radius[i][0] = 0.0;
		}
	}

	if(t==0) fclose(fout2);
}


/*	fuggveny a reszecskek tomegenek, indexenek es r_indexenek (a reszecske tavolsagaban levo gridcella indexe) frissitesere	*/
void Count_Mass(double radin[][2], double partmassindin[][4], double *massvecin, double t, int n) {

	int i, rindex;
	double rmid;	

	for (i = 0; i < n; i++) {	

  		rmid = (radin[i][0] - RMIN) / DD;     						/* 	The integer part of this gives at which index is the body			*/
		rindex = (int) floor(rmid+0.5);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
		if(rmid < 0) rindex = 0;
		if(isnan(rmid)) rindex = 0;	/*	neha az indexre - ha az mar RMIN-en belul van, nan-t ad, ezert abban az esetben is 0 lesz a rindex	 */

		if(n == PARTICLE_NUMBER) partmassindin[i][0] = massvecin[i];							/*	mass of the particles				*/
		partmassindin[i][1] = rindex;							/*	initial distance of the particles				*/
		if(t == 0) {
			partmassindin[i][2] = partmassindin[i][0];							/*	initial distance of the particles				*/
			partmassindin[i][3] = 0;
		}
 	}
}


/*	Fuggveny az adott futashoz mappa letrehozasara, igy egy adott futas mindig kulon mappaba kerul es nem kavarodnak ossze az adatok	*/
void Mk_Dir(char *nev) {

	int i, step, mappa;
	char comm[2048];

	step = 0;

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


void secondaryGrowth(double rad[][2], double radmicr[][2], double radsec[][2], double partmicind[][4], double partsecind[][4], double *massvec, double *massmicvec, double *masssecvec, double L) {

	int i, j, k=0;
	int step = 0;
	double temp[PARTICLE_NUMBER], temp2[4*PARTICLE_NUMBER], tempmic[PARTICLE_NUMBER], temp2ch1[PARTICLE_NUMBER], temp2ch2[PARTICLE_NUMBER], temp2ch3[PARTICLE_NUMBER], temp2ch4[PARTICLE_NUMBER];
	int hist[PARTICLE_NUMBER], histmic[PARTICLE_NUMBER];

/*	A porszemcse generalo program gridfelbontasa -- nem feltetlen egyezik meg a jelenlegi gridfelbontassal!	*/
	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1), rtempvec[PARTICLE_NUMBER];
	double temptemp[4*PARTICLE_NUMBER], masstemp[4*PARTICLE_NUMBER][2], sizetemp[4*PARTICLE_NUMBER],mtp[4*PARTICLE_NUMBER], indtemp[4*PARTICLE_NUMBER][2];

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {	
		temp2[i] = radsec[i][0];			// masodlagosan novesztett porszemcsek tavolsaga
		sizetemp[i] = radsec[i][1];
		temptemp[i] = 0;				// atmeneti vektor, amely elmeni, hogy hol noveszt uj reszecsket a program
		masstemp[i][0] = 0;				// atemeneti tomb, amely a novesztett reszecske tomeget tarolja el
		masstemp[i][1] = 0;				// atmeneti tomb, amely a novesztett reszecske tavolsagat tarolja el
		mtp[i] = masssecvec[i];
		if(temp2[i] != 0) step++;			// megszamolja, hogy hany novesztett reszecske van a fuggvenybe valo belepeskor
		indtemp[i][0] = partsecind[i][3];
	}

	for(i=0; i < PARTICLE_NUMBER; i++) {

		temp[i] = rad[i][0];				// cm-es porszemcsek tavolsaga
		tempmic[i] = radmicr[i][0];			// mikoronos porszemcsek tavolsaga
		rtempvec[i] = RMIN + i * dd;			// atmeneti r vektor
		rtempvec[i] = rtempvec[i] + dd / 2.0;		// ez adja a reszecskek eredeti helyet (a generalo programban ide helyeztuk el a reszecskeket	
		hist[i] = 0;					// cm-es porszemcsek hisztogramja
		histmic[i] = 0;					// mikoronos porszemcsek hisztogramja

		temp2ch1[i] = temp2[i];
		temp2ch2[i] = temp2[i+PARTICLE_NUMBER];
		temp2ch3[i] = temp2[i+2*PARTICLE_NUMBER];
		temp2ch4[i] = temp2[i+3*PARTICLE_NUMBER];
	}



/*	Megvizsgalja, hogy az adott tavolsagon hany reszecske talalhato	*/
/*	Valamiert csak ugy ad jo kepet, ha az if-es felteleket beirom.	*/
	for(i=0;i<PARTICLE_NUMBER;i++) {
		if(temp[i] != 0) histogram(temp[i],hist,dd);
		if(tempmic[i] != 0) histogram(tempmic[i],histmic,dd);
		if(temp2ch1[i] != 0) histogram(temp2ch1[i],hist,dd);
		if(temp2ch2[i] != 0) histogram(temp2ch2[i],hist,dd);
		if(temp2ch3[i] != 0) histogram(temp2ch3[i],hist,dd);
		if(temp2ch4[i] != 0) histogram(temp2ch4[i],hist,dd);
	}

	j=0;
	for(i=0;i<PARTICLE_NUMBER;i++) {
/*	Megvizsgalja, hogy kifogyott-e az adott tavolsagrol a cm-es por. Ha igen, viszont az adott tavolsagon van mikronos, akkor novel egy ujabb reszecsket, ha a mikronos por altal kepviselt tomeg az eredeti helyen levo tomeg szazadanal nagyobb	*/
		if(hist[i] == 0 && histmic[i] != 0) {
			if(rtempvec[i] != 0) {
				double sigcurr, siginit;
				sigcurr = partmicind[i][0]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
				siginit = partmicind[i][2]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
				if(sigcurr >= siginit / 100.) {
/*	Uj reszecske novelese --> az adott helyen levo mikronos reszecske tomegenek 75%-at kapja az uj reszecske, ezzel egyutt a mikronos reszecske altal kepviselt reszecske tomege 75%-kal csokken	*/
					massmicvec[i] = partmicind[i][0] * (1.-.75);	// a mikronos reszecske tomegenek csokkentese
					masstemp[j][1] = partmicind[i][0] * .75;	// novesztett reszecske altal kepviselt tomeg
					masstemp[j][0] = rtempvec[i];			// novesztett reszecske altal kepviselt tavolsag
					temptemp[j] = -1*rtempvec[i];			// atmeneti vektor, amely elmenti, hogy mely helyen noveszt uj reszecsket a program
					j++;
				}
			}
		}
	}

/*	A novesztett reszecske tavolsagat tartalmazo atmeneti vektor sorbarendezese	*/
	sort(temptemp,4*PARTICLE_NUMBER);	

	for(i=0; i < 4*PARTICLE_NUMBER-step; i++) {
		if((temptemp[i] != 0)) {
			radsec[i][0] = -1. * temptemp[i];
			temptemp[i] = 0;
		} else {
			radsec[i][0] = 0;
			radsec[i][1] = 0;
		}
	}

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {
		temptemp[i] = radsec[i][0];
		radsec[i][0] = 0;
		radsec[i][1] = 0;
	}

	sort(temptemp,4*PARTICLE_NUMBER);
	double index = find_max(indtemp,4*PARTICLE_NUMBER);

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {
		if(temptemp[i] >= RMIN) {
			radsec[i][0] = temptemp[i];
			for(j = 0; j < 4*PARTICLE_NUMBER; j++) {
				if(masstemp[j][0] == radsec[i][0]) {
					masssecvec[i] = masstemp[j][1];
					partsecind[i][0] = masssecvec[i];
				}
				if(temp2[j] == radsec[i][0]) {
					masssecvec[i] = mtp[j];
					radsec[i][1] = sizetemp[j];
					partsecind[i][3] = indtemp[j][0];
					partsecind[i][0] = masssecvec[i];
				} else {
					radsec[i][1] = 1e-6/AU2CM;
				}
			}
		}
	}
}


/*	Itt vegzi el az integralast, ha szukseg van ra	*/
void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec) {

   	int linesout;
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A reszecskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe, amennyiben a futas soran szamol a program driftet	*/
	if(optdr == 1.) {
		PARTICLE_NUMBER = reszecskek_szama(linesout,filenev1);  		
	} else {
		PARTICLE_NUMBER = 0;
	}

	char mv[1024], dens_name[1024],size_name[1024], porout[1024], poroutmicr[1024], massout[1024], dust_name[1024], dust_name2[1024];
   	double radius[PARTICLE_NUMBER][2],radiusmicr[PARTICLE_NUMBER][2],radiussec[4*PARTICLE_NUMBER][2],radius_rec[PARTICLE_NUMBER][2];
	double massvec[PARTICLE_NUMBER], massmicrvec[PARTICLE_NUMBER], masssecvec[4*PARTICLE_NUMBER];
	double max, min, max2, min2;
	int i;

/*	A reszecskek adatainak beolvasasa file-bol, ha az optdr erteke 1, egyebkent nem szamol a program driftet		*/
	if(optdr == 1.) {

		por_be(radius,radiusmicr,massvec,massmicrvec);			/*	porreszecskek adatainak beolvasasa	*/ 
		int dummy;
		sprintf(mv,"cp %s %s/",filenev1,nev);				/*	az adott helyre, ahova eppen menti a futast eredmenyeit a program, a porreszecskek adatait tartalmazo file atmasolasa, igy kesobb is visszanezheto, hogy mik voltak a kezdeti adatok	*/
		dummy = system(mv);						/*	itt masolja at a file-t			*/	
/*	az aktualis mappaban a pormozgas.dat file letrehozasa: ebbe kerul be a porreszecske tavolsaga es indexe, valamint az adott idolepes	*/
		snprintf(porout,1024,"%s/pormozgas.dat",nev);			
/*	ha 2 populacios a futas, akkor a mikronos pornak is letrehoz egy pormozgas file-t, ebbe kerul be a tavolsag, index es az ido */
		if(opttwopop == 1.) snprintf(poroutmicr,1024,"%s/pormozgasmic.dat",nev);
/*	tomegnovekedesi file letrehozasa az aktualis mappaba - ez lehet, hogy egy kulon opcio lesz a kimeneti adatok meretenek csokkentesere	*/
		snprintf(massout,1024,"%s/mass.dat",nev);
   		fout = fopen(porout,"w");
		if(opttwopop == 1.) foutmicr = fopen(poroutmicr,"w");
   		massfil = fopen(massout,"w");

	}

  	double t = 0.0;
   	double t_integration = TMAX * 2.0 * M_PI; 						/*	numerikus integralas idotartama		*/
	double deltat = time_step(rvec)/5.;							/*	idolepes	*/

	if(DT <= deltat && DT != 0) deltat = DT;	

	double partmassind[PARTICLE_NUMBER][4], partmassmicrind[PARTICLE_NUMBER][4], partmasssecind[4*PARTICLE_NUMBER][4];
	double L = 0.;
	double sigmad[PARTICLE_NUMBER], sigmadm[PARTICLE_NUMBER], sigmads[4*PARTICLE_NUMBER], rdvec[PARTICLE_NUMBER], rmicvec[PARTICLE_NUMBER], rsvec[4*PARTICLE_NUMBER];
	int num = PARTICLE_NUMBER;
	
	double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0; // a kulso es belso dze-n felgyulemlett por mennyisege -- bemeneti adat (Print_Mass fuggvenybe)
	double masstempiout = 0, masstempoout = 0, massmtempiout = 0, massmtempoout = 0; // a kulso es belso dze-n felgyulemlett por mennyisege -- kimeneti adat (Print_Mass fuggvenybol)
	double tavin,tavout;

	for(i = 0; i < 4*PARTICLE_NUMBER; i++) {
		radiussec[i][0] = 0;
		radiussec[i][1] = 0;
		partmasssecind[i][0] = 0;
		partmasssecind[i][1] = 0;
		masssecvec[i] = 0;
		sigmads[i] = 0;
		rsvec[i] = 0;
	}

	if(opttwopop == 0) {
		for(i = 0; i < PARTICLE_NUMBER; i++) {
			radiusmicr[i][0] = 0;
			radiusmicr[i][1] = 0;
			partmassmicrind[i][0] = 0;
			partmassmicrind[i][1] = 0;
			massmicrvec[i] = 0;
		}
	}

   	do {

/*	Ha van drift:	*/	

		if(optdr == 1.) {

			if(opttwopop == 1) secondaryGrowth(radius,radiusmicr,radiussec,partmassmicrind,partmasssecind,massvec,massmicrvec,masssecvec, L);
/*	A minimum kereseshez letrehozza a cm-es reszecskek tavolsaganak reciprokat	*/
			for (i=0; i < PARTICLE_NUMBER; i++) {
				if (radius[i][0] > 0. && radius[i][0] > RMIN) {
					radius_rec[i][0] = 1. / radius[i][0];
				} else {
					radius_rec[i][0] = 0.;
				}
			}
		
			max = find_max(radius,PARTICLE_NUMBER);					/*	Megkeresi, hogy melyik a legtavolabbi cm-es reszecske a kozponti csillagtol	*/
			min = find_max(radius_rec,PARTICLE_NUMBER);				/*	Megkeresi a tavolsag reciprokanak maximumat, azaz a legkisebb tavolsagra levo cm-es reszecsket	*/
			min = 1. / min;

			double mint, maxt;

/*	ha 2 populacios a szimulacio, a fentihez hasonloan megkeresi a legnagyobb es a legkisebb tavolsagra levo mikronos reszecsket	*/
			if(opttwopop == 1) {
				for (i=0; i < PARTICLE_NUMBER; i++) {
					if (radiusmicr[i][0] > 0. && radiusmicr[i][0] > RMIN) {
						radius_rec[i][0] = 1. / radiusmicr[i][0];
					} else {
						radius_rec[i][0] = 0.;
					}
				}

				max2 = find_max(radiusmicr,PARTICLE_NUMBER);					/*	Megkeresi, hogy melyik a legtavolabbi reszecske a kozponti csillagtol	*/
				min2 = find_max(radius_rec,PARTICLE_NUMBER);
				min2 = 1. / min2;

/*	megnezi, hogy mely reszecske van a legkozelebb, illetve legtavolabb (mikronos, vagy cm-es)	*/
				mint = find_min(min,min2,HUGE_VAL);
				maxt = find_min(1. / max, 1./max2,HUGE_VAL);
				maxt =  1./ maxt;
			} else {
				mint = min;
				maxt = max;
			}

/*	Ha a legtavolabbi reszecske tavolsaga nagyobb, mint RMIN (azaz meg a szimulacio tartomanyaban van), es a min es a max nem egyenlo (azaz nem gyult pl. ossze az osszes porreszecske 1 helyen), akkor a program tovabb szamol, egyebkent a futas leall	-- ez persze 1 dze eseten mukodik, meg kell oldani, hogy 2 dze eseten is lealljon akkor, ha az osszes por osszegyult a nyomasi maximumokban - meg kell persze csinalni, hogy ez is opcionalis legyen	*/
			if(maxt >= RMIN && mint != maxt) {

				double time = t / 2.0 / M_PI;

				if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){

/*	Az adatok kiirasahoz szukseges file-ok neveinek elmentese	*/
					if (t==0) {
						snprintf(dens_name,1024,"%s/surface.dat",nev);
					} else {
						if(optev == 1) {
							snprintf(dens_name,1024,"%s/dens.%d.dat",nev,(int)L);
						}
					}

					snprintf(dust_name,1024,"%s/dust.%i.dat",nev,(int)L);
					snprintf(dust_name2,1024,"%s/dustmic.%i.dat",nev,(int)L);
					snprintf(size_name,1024,"%s/size.%d.dat",nev,(int)L);
	
					if(t==0) {
/*	A reszecskek tomeget tartalmazo tomb inicializalasa		*/
						Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER);
						if(opttwopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER);
						if(opttwopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER);

/*	Ha van tomegnovekedes, akkor a por feluletisurusegenek kiszamolasa itt tortenik	*/
						if(optgr == 1.) {			
							Get_Sigmad(L,max,min,radius,radiusmicr,radiussec,sigmad,sigmadm,sigmads,massvec,massmicrvec,masssecvec,rdvec,rmicvec,rsvec);
						}
					}
		
/*	A sigma, p, dp kiirasa egy file-ba	*/
					if(optev == 1 || time == 0) {
						Print_Sigma(dens_name, rvec, sigmavec, pressvec, dpressvec);
					}

/*	Ha szamol a futas driftet, itt irja ki a reszecskek tavolsagat es meretet	*/
					if(optdr == 1) {
						Print_Pormozg_Size(size_name,(int)L,radius,radiusmicr,rvec,t);
					}

/*	A tomegnovekedesi file-ba az adatok kiirasa	*/
					masstempiout = 0; 
					massmtempiout = 0;
					masstempoout = 0;
					massmtempoout = 0;

					Print_Mass(L,rvec,partmassind,partmassmicrind,partmasssecind,t,dpressvec,masstempiin,masstempoin,massmtempiin,massmtempoin,&masstempiout,&masstempoout,&massmtempiout,&massmtempoout,&tavin,&tavout);

					if(r_dze_i != tavin) {
						masstempiin = masstempiout; 
						massmtempiin = massmtempiout;
					}
					if(r_dze_o != tavout) {
						masstempoin = masstempoout;
						massmtempoin = massmtempoout;
					}

/*	Ha van pornovekedes, kiirja a por felultisuruseget egy file-ba --> a pornovekedeshez szukseges egyaltalan ezt kiszamolni!	*/
					if(optgr == 1.) {
						Print_Sigmad(dust_name,dust_name2,mint,rdvec,rmicvec,sigmad,sigmadm);
					}

					L = L+(double)(TMAX/WO);

				}

/*	Ha az optev erteke 1, akkor megoldja minden lepesben a sigmara vonatkozo diffuzios egyenletet	*/
				if(optev == 1.) {
					Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
				} 

/*	A reszecskek tomeget tartalmazo tomb adatainak frissitese	*/
				Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER);
				if(opttwopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER);
				if(opttwopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER);

/*	Ha van reszecskenovekedes, akkor kiszamolja a por feluletisuruseget	*/
				if(optgr == 1.) {
					Get_Sigmad(L,max,min,radius,radiusmicr,radiussec,sigmad,sigmadm,sigmads,massvec,massmicrvec,masssecvec,rdvec,rmicvec,rsvec);
				}

				int optsize = 0;		// ezt ki lehetne siman valtani opttwopop-pal!
/*	A cm-es reszecskek eseten az optsize erteke 0	*/
/*	Itt szamolja ki a program a cm-es reszecskek uj tavolsagat (es meretet, ha kell)	*/
				Get_Radius(nev,optsize,radius,pressvec,dpressvec,sigmavec,sigmad,rvec,ugvec,deltat,t,PARTICLE_NUMBER);
 
/*	Ha a futas 2 populacios, akkor az optsize erteke 1	*/
/*	Itt szamolja ki a program a mikoronos reszecskek uj tavolsagat	*/
				if(opttwopop == 1.) {
					optsize = 1;
					Get_Radius(nev,optsize,radiusmicr,pressvec,dpressvec,sigmavec,sigmad,rvec,ugvec,deltat,t,PARTICLE_NUMBER); 
					optsize = 2;
					Get_Radius(nev,optsize,radiussec,pressvec,dpressvec,sigmavec,sigmad,rvec,ugvec,deltat,t,4*PARTICLE_NUMBER); 
				}

				t = t + deltat;						/*	Idoleptetes	*/

			} else {				/*	Ha a legmesszebbi reszecske tavolsaga mar nem nagyobb, vagy egyenlo, mint RMIN, vagy a legkisebb tavolsagra levo reszecske tavolsaga, akkor a program "figyelmezteto szoveg" mellett sikeresen kilep, nem fut "feleslegesen" tovabb.	*/
				printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",L);
				exit(EXIT_SUCCESS);
			}	
	
		} else {	/*	Ez az az eset, ha a program nem szamol driftet, azaz csak a gaz feluletisurusegenek fejlodesere vagyunk kivancsiak	*/

			double time = t / 2.0 / M_PI;

			if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){	

				if (t==0) {
					snprintf(dens_name,1024,"%s/surface.dat",nev);
				} else {
					snprintf(dens_name,1024,"%s/dens.%d.dat",nev,(int)L);
				}
	
				Print_Sigma(dens_name, rvec, sigmavec, pressvec, dpressvec);

				L = L+(double)(TMAX/WO);
			}

			Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
			t = t + deltat;						/*	Idoleptetes		*/
		}

   	} while (t <= t_integration);
  

/*	Az idoleptetes leteltevel a program sikeresen kilep	*/
	printf("\n\nA program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n",max); 


}


/*	Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok	*/
void infoCurrent(char *nev) {

	char out[1024];

/*	A TCURR tartalmazza a a szimulacio inditasanak pillanataban az idot: honap, nap, ora, perc, mp formatumban	 */
	sprintf(out,"run_%i.dat",(int)TCURR);
	jelfut = fopen(out,"w");
	
/*	Abba a mappaba hoz letre egy file-t, ahonnan inditottuk a programot. A file-ban leolvashatjuk, hogy a szimulacio kimenete mely mappaban talalhato, illetve a szimulaciorol nehany informacio -- ezt ki fogom meg egeszitani	*/
	fprintf(jelfut,"A jelenlegi futás a %s mappaban taláható!\n",nev);
	fprintf(jelfut,"\n\nA korong paraméterei:\nRMIN: %lg, RMAX: %lg\nSIGMA0: %lg, SIGMA_EXP: %lg, flaring index: %lg\nALPHA_VISC: %lg, ALPHA_MOD: %lg\nR_DZE_I: %lg, R_DZE_O: %lg, DR_DZEI: %lg, DR_DZE_O: %lg  (*** R_DZE_I/O = 0, akkor azt a DZE-t nem szimulálja a futás! ***)\n\n\n",RMIN,RMAX,SIGMA0,SIGMAP_EXP,FLIND,alpha_visc,a_mod,r_dze_i,r_dze_o,Dr_dze_i,Dr_dze_o);
	fprintf(jelfut,"A központi csillag tömege: %lg M_Sun\n",STAR);
	fclose(jelfut);

}


int main(int argc, const char **argv) {
   
	options_t def;				/*	A letrehozott struktura elemeire .def-el lehet majd hivatkozni	*/
	create_default_options(&def);		/*	A struktura elemeinek feltoltese alapertelmezett parameterekkel	*/
/*	A terminalbol keresi a kapcsolok segitsegevel a struktura elemeire vonatkozo adatokat. A beolvasas kimenetet eltarolja egy integerbe	*/
	int retCode = parse_options(argc, argv, &def);
/*	Ha a parameterek beolvasasa sikertelen volt, akkor a program kilep	*/	
	if (0 != retCode) {
		exit(retCode);
	}

/*	A globalis valtozok feltoltese a parameterek ertekevel	*/
	optev = def.evol;
	optdr = def.drift;
	optgr = def.growth;
	opttwopop = def.twopop;
	fFrag = def.ffrag;
	uFrag = def.ufrag;
	inputsig = def.input;
	DT = def.tStep;

	if(optinp == 0) {

		int lout, nout;
		lout = 0, nout = 0;

/*	A sigmat tartalmazo file sorainak szama elmentve egy integerbe	*/
		nout = reszecskek_szama(lout,inputsig); 
		NGRID = nout;

	} else {

		NGRID = def.ngrid;

	}

   	double sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], ugvec[NGRID+2];
	char dens_name[1024], nev[1024], mv[1024];

/*	A korong parametereinek beolvassa				*/
	disk_param_be(&SIGMA0, &SIGMAP_EXP, &RMIN, &RMAX, &r_dze_i, &r_dze_o, &Dr_dze_i, &Dr_dze_o, &a_mod, &PDENSITY, &PDENSITYDIMLESS, &alpha_visc,&STAR,&FLIND);
	DD = (RMAX - RMIN) / (NGRID - 1);						/*	rácsfelbontás	*/

/*	Ha nincs bemeneti sigma file	*/
	if(optinp == 0) {

		sigIn(sigmavec,rvec);
		Perem(rvec);
		Perem(sigmavec);
		
	} else {

/*	Kezdeti profilok betoltese a megadott vektorokba	*/
		load_R(rvec);
		Initial_Profile(sigmavec,rvec);

	}

	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	Initial_Ugas(sigmavec,rvec,ugvec);

	timePar(&TMAX,&WO,&TCURR);

/*	az optdze globalis valtozo ertekenek megadasa: ha van belso nyomasi maximum is, akkor az erteke 1, egyebkent nulla. Ennek a tomegnovekedes szamolasanak szempontjabol van lenyege	*/ 
	optdze = 0;
	if(r_dze_i != 0) optdze = 1.;

/*	Mappa letrehozasa az adatok eltarolasahoz	*/
	Mk_Dir(nev);								

/*	Az aktualis mappaba a kezdeti adatokat tartalmazo file-ok atmasolasa - kesobb ezek az adatok visszanezhetok igy		*/
	int dummy;
	sprintf(mv,"cp %s %s/",filenev2,nev);
	dummy = system(mv);	
	sprintf(mv,"cp %s %s/",filenev3,nev);
	dummy = system(mv);

/*	Sigma file beolvasasa, ha szukseges	*/
	if(optinp == 0) {
		sprintf(mv,"cp %s %s/",inputsig,nev);
		dummy = system(mv);
	}

/*	Abban a mappaban, ahol a futast inditottuk, egy file letrehozasa, amely az aktualis futasrol infokat ir ki (hol talalhatoak a kimeneti file-ok, es milyen parameterei vannak pl. a korongnak az adott futas eseten	*/
	infoCurrent(nev);

/*	Ha nincs sigma fejlodes es drift, akkor a kezdeti profilt kiirja egy file-ba es kilep "figyelmeztetes" mellett	*/
	if(optev == 0. && optdr == 0.) {
		printf("A megadott opciok szerint nem szamol sem sigmat, sem driftet, ezert a progam kilep!\n\nA kezdeti file-ok a %s mappaban talalhatoak!\n",nev);
/*	t=0-ban kiirja a sigma-t, a nyomast es a nyomasderivaltjat	*/
		snprintf(dens_name,1024,"%s/surface.dat",nev);
		Print_Sigma(dens_name,rvec,sigmavec,pressvec,dpressvec);
	} else {
		tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec);
	}

	return 0;

}
