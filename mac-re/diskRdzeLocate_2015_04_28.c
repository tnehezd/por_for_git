//peldaprogram arra, hogyan lehet kiszamolni adott parameterekkel, hogy hova esik a dze kulso pereme

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>


#define SDCONV			1.12521e-7
#define ICEFACTOR       	3.0
#define SNOWLINE        	2.7
#define G			1.0
#define G2			G*G
#define AUPDAY2CMPSEC   	1.7314568e8    	//conversion between the velocities 
#define CMPSECTOAUPYRP2PI	3.35725e-07	//conversion between the velocities
#define GRPCM32MSUNAU3		1.68329e6
#define SDCRIT		200.*SDCONV		/*	kritikus feluleti suruseg, a dze kulso hatara	*/



double ALPHA, STAR, DISKMASS, HASP,  FLIND, DUST2GASR;
double RMIN, RMAX;
double NGRID;
double r_dze_i;
double Dr_dze_i;
double r_dze_o;
double Dr_dze_o;
double a_mod;
double SIGMAP_EXP;
long double SIGMA0;
long double output;

double SDRCRIT;						/*	az a tavolsag, ahol a feluleti suruseg a kritikus ala csokken, es ujra teljesen ionizaltta valik a korong	*/

FILE *fin, *fout, *fout2;

typedef struct options {
	//! Number of particles
	double	n;
	//! Inner edge of the gaseous disk
	double	ri;
	//! Outer edge of the gaseous disk
	double	ro;
	//! Surface density of the gaseous disk
	long double	sigma0;
	//! Surface density of the gaseous disk in gr/cm/cm
	long double	sigma0cgs;
	//! The index of the profile of the initial surface density of the gaseous disk
	double	index;
	//! The inner deadzone edge
	double 	rdze_i;
	//! The outer deadzone edge
	double 	rdze_o;
	//! The width of the viscosity jump at the inner deadzone edge
	double 	drdze_i;
	//! The width of the viscosity jump at the outer deadzone edge
	double 	drdze_o;
	//! The alpha parameter (Shakura & Sunyaev (1973))
	double	alphaParam;
	//! Viscosity reduction parameter
	double 	amod;
	//! The aspect ratio of the disk (h = H/r)
	double	h;
	//! Flaring index
	double	flind;
	//! The mass of the central star (solar mass)
	double	m0;
	//! The mass of the disk (solar mass)
	double	md;
	//!Option for output file to overwrite or add output
	double output;

} options_t;


void create_default_options(options_t *opt) {



	opt->n				 = 1000;	
	opt->ri				 = 0.1;		
	opt->ro				 = 5.0;
	opt->sigma0			 = 0.0001;
	opt->sigma0cgs			 = 1.0e-4/SDCONV;
	opt->index			 = 0.5;
	opt->rdze_i			 = 0.0;		// rdze_i: nincs belso dze
	opt->rdze_o			 = 0.0;		// rdze_o: nincs kulso dze
	opt->drdze_i			 = 0.0;		
	opt->drdze_o			 = 0.0;
	opt->alphaParam		 	 = 1.0e-2;
	opt->amod			 = 0.01;
	opt->h				 = 5.0e-2;
	opt->flind			 = 0.0;		// flat korong!
	opt->m0				 = 1.0;		
	opt->md				 = 0.01;
	opt->output			 = 0.;
	
}


int parse_options(int argc, const char **argv, options_t *opt, double *par, double *dpar){
	int i = 1;
	double temp1 = 0., temp2 = 0., temp3 = 0.;

	while (i < argc) {
		char *p = malloc(sizeof(argv[i]));
   		strcpy(p, argv[i]);

		// Number of intervals
		if(strcmp(p, "-n") == 0) {
			i++;
			opt->n = atof(argv[i]);
		}
		else if (strcmp(p, "-ri") == 0) {
			i++;
			opt->ri = atof(argv[i]);
		}
		else if (strcmp(p, "-ro") == 0) {
			i++;
			opt->ro = atof(argv[i]);
		}
		else if (strcmp(p, "-sigma0") == 0) {
			temp1 = 1.;
			i++;
			opt->sigma0 = atof(argv[i]);
		}
		else if (strcmp(p, "-sigma0cgs") == 0) {
			temp1 = 1.;
			i++;
			opt->sigma0 = atof(argv[i])/SDCONV;
			opt->sigma0cgs = atof(argv[i]);
		}
		else if (strcmp(p, "-index") == 0) {
			temp3 = 1.0;
			i++;
			opt->index = atof(argv[i]);
		}
		else if (strcmp(p, "-rdzei") == 0) {
			i++;
			opt->rdze_i = atof(argv[i]);
		}
		else if (strcmp(p, "-rdzeo") == 0) {
			i++;
			opt->rdze_o = atof(argv[i]);
		}
		else if (strcmp(p, "-drdzei") == 0) {
			i++;
			opt->drdze_i = atof(argv[i]);
		}
		else if (strcmp(p, "-drdzeo") == 0) {
			i++;
			opt->drdze_o = atof(argv[i]);
		}
		else if (strcmp(p, "-alpha") == 0) {
			i++;
			opt->alphaParam = atof(argv[i]);
		}	
		else if (strcmp(p, "-amod") == 0) {
			i++;
			opt->amod = atof(argv[i]);
		}
		else if (strcmp(p, "-h") == 0) {
			i++;
			opt->h = atof(argv[i]);
		}

		else if (strcmp(p, "-flind") == 0) {
			i++;
			opt->flind = atof(argv[i]);
		}
		else if (strcmp(p, "-m0") == 0) {
			i++;
			opt->m0 = atof(argv[i]);
		}
		else if (strcmp(p, "-md") == 0) {
			temp2= 1.;
			i++;
			opt->md = atof(argv[i]);
		}
		else if (strcmp(p, "-merge") == 0) {
			i++;
			opt->output = atof(argv[i]);
		}
		else {
			printf("\n\n**** Invalid switch on command-line: %s! ****\n\n\n",p);
			printf("**** Try following parameters: ****\n\n-n: number of grids (particles)\n-ri\ninner edge of disk\n-ro\nother edge of disk\n-sigma0\nsigma at r=1AU\n-sigma0cgs\nsigma at r=1AU in cgs\n-index\nsigma exponent\n-rdzei\ninner DZE\n-rdzeo\nother DZE\n-drdzei\nthe width (in H, eg. 2H, 1.5H, etc) of the viscosity jump at the inner deadzone edge (if not given, but rdzei is given, set to 2*H)\n-drdzeo\nthe width (in H, eg. 2H, 1.5H, etc) of the viscosity jump at the outer deadzone edge (if not given, but rdzeo is given, set to 2*H)\n-alpha\nalpha viscosity\n-amod: viscosity reduction\n-h\naspect ratio\n-flind\nflaring index\n-m0\nmass of the central star (in solar masses)\n-md\nmass of the disk\n-eps\ndust to gas ratio\n\n\n");
			return 1;
		}
		i++;
		free(p);
	}
	
	*par = temp1;
	*dpar = temp2;
	
	return 0;
}



/*	sigma kiszámolása r helyen M_NAP / AU / AU-ban	*/
double sigma_gas(double r) {

	return SIGMA0 * pow(r,SIGMAP_EXP);
	
}



/*	sigma0 kiszámolása (illetve megadása) M_NAP / AU / AU-ban	*/
double sigma_null() {
	
	double alpha2 = SIGMAP_EXP + 2;
	return (alpha2 / (2.0 * M_PI)) * DISKMASS / (pow(RMAX,alpha2) - pow(RMIN,alpha2));
	
}


/*	mivel a sig csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne SDCRIT az erteke	*/
/*	solving a*x + b = y (here a = r, y = dp)	*/
double find_r(double r1, double r2, double dp1, double dp2) {

	double a, b, r;
	a = (dp2 - dp1) / (r2 - r1);
	b = dp1 - a * r1;
	r = (SDCRIT - b) / a;

	return r;

}



int main(int argc, const char **argv) {


	options_t def;
	create_default_options(&def);
	double par, par2, par3, par4, par5;
	int retCode = parse_options(argc, argv, &def, &par, &par2);

	if (0 != retCode) {
		exit(retCode);
	}

	NGRID = def.n;

	SIGMAP_EXP = def.index;
	SIGMAP_EXP = SIGMAP_EXP*(-1.0);

	ALPHA = def.alphaParam;
	STAR = def.m0;
	DISKMASS = def.md; 
	HASP = def.h;
	FLIND = def.flind;
	RMIN = def.ri;
	RMAX = def.ro;
	a_mod = def.amod;
	SIGMA0 = def.sigma0;
	par3 = def.output;

	if(par2 == 1) SIGMA0 = sigma_null();

	char filnev[1024];
	sprintf(filnev,"rdzrloc_%lg.dat",SIGMAP_EXP);

	if(par3 == 0) fout = fopen(filnev,"w");
	if(par3 == 1) fout = fopen(filnev,"a");

	double rcritin = 0, rcritout = 0, scritin = 0, scritout = 0;
	double sigmavec[(int)NGRID], r, DD;  
	int i;

	DD = (RMAX-RMIN)/(NGRID - 1);

  	for(i = 0; i < NGRID; i++) {
		r = RMIN + i * DD;
    		sigmavec[i] = SIGMA0 * pow(r,SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/

		if(sigmavec[i]<=SDCRIT && rcritout == 0) {
			rcritout = r;	
			rcritin = r-DD;
			scritout = sigmavec[i];
			scritin = sigmavec[i-1];
		}
  	}

	SDRCRIT = find_r(rcritin,rcritout,scritin,scritout);		/*	kiszamolja azt a helyet, ahol a feluletisuruseg eppen a kritikus erteket veszi fel	*/

	fprintf(fout,"%lg %lg %lg %lg %Lg %lg %lg %lg\n",RMIN,RMAX,SDRCRIT,DISKMASS,SIGMA0,SDCRIT,STAR,SIGMAP_EXP);
	fclose(fout);
	

	return 0;

}
