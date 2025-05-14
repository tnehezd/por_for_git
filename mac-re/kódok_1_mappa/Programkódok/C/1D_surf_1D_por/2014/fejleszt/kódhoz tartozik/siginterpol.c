#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

const char *inputsig;
const char *outputsig;

double RMIN, RMAX, DD;
int NGRID, optout;

/*	A futasok soran valtoztathato opciok letrehozasa egy strukturaval	*/
typedef struct options {

	// Number of grid cells
	int ngrid;
	// RMIN
	double rmin;
	// RMAX
	double rmax;
	// Input sigma file
	const char *input; 
	// Output sigma file
	const char *output;
	// Option for output
	int opto;

} options_t;


/*	A struktura elemeinek feltoltese alapertelmezett ertekekkel	*/
void create_default_options(options_t *opt) {

	opt->ngrid		 = 0;
	opt->rmin   		 = 0.;
	opt->rmax 		 = 0.;
	opt->input		 = "";
	opt->output		 = "";
	opt->opto		 = 0;
	
}


/*	A valtoztathato parameterek beolvasa terminalbol	*/
int parse_options(int argc, const char **argv, options_t *opt){
	int i = 1;

	while (i < argc) {
		char *p = malloc(sizeof(argv[i]));	/*	Ahhoz, hogy "szoveget" (char, mert a C nem tud stringet kezelni) tudjunk beolvasni, es ossze tudjuk vetni a lentebb megadott kapcsolokkal (pl. -drift), le kell foglalni a memoriateruletet p-nek. Ennek viszont elore nem tudjuk a meretet, ezert dinamikusan foglaljuk le azt	*/
   		strcpy(p, argv[i]);			/*	A p-be elmentjuk az adott argumentumot az strcpy paranccsal	*/	

		if (strcmp(p, "-o") == 0) {
			i++;
			opt->output = argv[i];
			opt->opto = 1;
		}
		else if (strcmp(p, "-i") == 0) {
			i++;
			opt->input = argv[i];
		}
		else if (strcmp(p, "-ri") == 0) {
			i++;
			opt->rmin = atof(argv[i]);
		}
		else if (strcmp(p, "-ro") == 0) {
			i++;
			opt->rmax = atof(argv[i]);
		}
		else if (strcmp(p, "-n") == 0) {
			i++;
			opt->ngrid = atoi(argv[i]);
		}
		else {				/*	ha nem sikerul a parancssorbol a beolvasas, akkor a program kilep es feldobja, hogy mely kapcsolo nem volt jo, helyette mit probaljon meg	*/
			printf("\n\n**** Invalid switch on command-line: %s! ****\n\n\n",p);
			return 1;
		}
		i++;
		free(p);			/*	memoria foglalas miatt fel is kell szabaditanunk a memoria teruletet	*/
	}

	return 0;
}




/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int reszecskek_szama(int numout, const char *filenev){

	char c;
	FILE *fin1;
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
void sigIn(double *sigvec, double *rvec, int n) {

	double sig,r;
	int i;

	FILE *densin;
	densin = fopen(inputsig,"r");
	
	for(i = 0; i < n; i++) {
           	if(fscanf(densin,"%lg  %lg",&r,&sig) == 2) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 3, mert 3 oszlopot szeretnénk beolvasni.	*/
			rvec[i] = r;				/*	r vektor	*/
			sigvec[i] = sig; 			/*	adott lepeskozonkent irja ki a file-t, megvan, hogy hany evente, tudjuk, hogy meddig fusson a kod, igy egy egyszeru osztassal meg lehet adni, hogy az mindig hanyadik idolepes	*/
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(densin);
	
}

/*	r vektor (gridcellák) inicializálása	*/
void load_R(double *rint) {
	
	int i;
 	for(i = 0; i < NGRID; i++) {						/*	load an array of radii	*/
 		rint[i] = RMIN + i * DD;
	}

}



/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out, double n) {

	double rmid, rindex, coef1, temp, rd;
	int index; 

	double a = log10(RMIN);
	double b = log10(RMAX);
	double dd = (b-a) / (n-1);
	double ytemp = log10(pos);
	double stemp = log10(pos);

     	rmid = (ytemp - a)/dd;
	index = (int) floor(rmid);					/* 	ez az rmid egesz resze	(kerekites 0.5-tol)			*/

	double f;

//y = 10^[log(y1) + (log(x) - log(x1)) * (log(y2) - log(y1)) / (log(x2) - log(x1))]
 

	f = ytemp - log10(rvec[index]);
	f = f * (log10(invec[index+1]) - log10(invec[index]));
	f = f / (log10(rvec[index+1]) - log10(rvec[index]));
	f = f + log10(invec[index]);

	f = pow(10,f);

	*out = f;

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

	Parabola(vec, NGRID - 4, NGRID - 3, NGRID - 2, &a, &b, &c, DD);
	vec[NGRID-1] = a * (RMAX) * (RMAX) + b * (RMAX) + c;

}


int main(int argc, const char **argv) {

	FILE *outf;

	options_t def;				/*	A letrehozott struktura elemeire .def-el lehet majd hivatkozni	*/
	create_default_options(&def);		/*	A struktura elemeinek feltoltese alapertelmezett parameterekkel	*/
/*	A terminalbol keresi a kapcsolok segitsegevel a struktura elemeire vonatkozo adatokat. A beolvasas kimenetet eltarolja egy integerbe	*/
	int retCode = parse_options(argc, argv, &def);
/*	Ha a parameterek beolvasasa sikertelen volt, akkor a program kilep	*/	
	if (0 != retCode) {
		exit(retCode);
	}

/*	A globalis valtozok feltoltese a parameterek ertekevel	*/
	inputsig = def.input;
	outputsig = def.output;
	RMIN = def.rmin;
	RMAX = def.rmax;
	NGRID = def.ngrid;
	optout = def.opto;
	
	int lout, nout;
	lout = 0, nout = 0;

/*	A sigmat tartalmazo file sorainak szama elmentve egy integerbe	*/
	nout = reszecskek_szama(lout,inputsig); 

   	double sigmavec[nout], rvec[nout];
	sigIn(sigmavec,rvec,nout);

	if(NGRID == 0) NGRID = nout;
	if(RMIN == 0) RMIN = rvec[0];
	if(RMAX == 0) RMAX = rvec[nout-1];

	DD = (RMAX-RMIN) / (NGRID - 1);

   	double sigint[NGRID], rint[NGRID];

	printf("GRID: %i  rmin: %lg  rmax: %lg  N: %i DD: %lg   opt: %i\n",nout,RMIN,RMAX,NGRID,DD,optout);

	load_R(rint);

	int i;
	
	double y, interp;

	for(i=0; i<NGRID-1; i++) {
		y = rint[i];
		interpol(sigmavec,rvec,y,&interp,nout);
		sigint[i] = interp;
	}

	Perem(sigint);

	if(optout == 1) {
		outf = fopen(outputsig,"w");
		for(i=0; i < NGRID; i++) fprintf(outf,"%lg	%lg\n",rint[i],sigint[i]);
		fclose(outf);
	}
	

	return 0;

}
