/*	A porszemcse ulepedese koagulacioval laminaris akkrecios korongban	*/

#include<stdio.h>
#include<math.h>

#define	ASP_RATIO	0.05		/*	aspect ratio					*/
#define SIGMA0		0.0001		/*	feluleti suruseg 1 AU-nal (M_NAP / AU / AU)	*/
#define SD_EXP		-0.5		/*	feluleti suruseg kitevoje			*/
#define	M_STAR		1.0		/*	csillag tomege naptomegben			*/
#define	G		0.01720209895	/*	Gauss-allando					*/
#define	DUST2GAS	0.01		/*	por gaz arany					*/
#define PDENSITY	3.0		/*	a porszemcse surusege (g / cm / cm / cm)	*/
#define	PRADIUS		1.0		/*	a porszemcse sugara cm-ben			*/
#define AU2CM		1.496e13	/*	1 AU --> 0.496e13 cm konverzio			*/
#define SUN2GR		1.9892e33	/*	M_nap --> gramm	konverzio			*/

/*	feluleti suruseg	*/
double sigma(double r) {

	return SIGMA0 * pow(r,SD_EXP);				/*	diszkprofil			*/

}


/*	termalis sebesseg	*/
double v_thermal(double r) {

	double c_s;
	c_s = ASP_RATIO * G * G * sqrt( M_STAR / r);		/* 	H * omega = h * v_kep --> AU / D_KOZEP	*/
	return 2.0 * sqrt(2.0 / M_PI) * c_s;			/*	sqrt(8.0 / PI) 			*/

}


/*	gaz surusege a csillagtol, valamint a korong fosikjatol mert tavolsag fuggvenyeben (Gauss-profil z-ben)	*/
double gas_density(double r, double z) {

	double density, H;
	H = ASP_RATIO * r;					/*	nyomasi skalamagassag (AU)	*/
	density = sigma(r) / sqrt(2.0 * M_PI) / H * exp(-z * z / (2.0 * H * H));	
	return density;						/* 	M_STAR / AU / AU / AU		*/

}

/*	kepler korfrekvencia	*/
double omega(double r) {
	
	double r3;
	r3 = r * r * r;
	return G * sqrt(M_STAR / (r3));

}

/*	diffegyenlet jobb oldal fuggvenye, az y tombben vannak az m, z, az ynew-ban az mnew,ynew elemek	
y[0] = m, y[1] = z; y_new[0] = mnew, y_new[1] = znew	*/
void eq_rhs(double r, double *y, double *ynew) {

	double omega2, density, vtherm, conv_fact, m, z, m_new, z_new;

	m = y[0];				/*	por reszecske tomege			*/
	z = y[1];				/*	por reszecske fosiktol valo tavolsaga	*/
	m_new = ynew[0];
	z_new = ynew[1];

	density = gas_density(r,z);
	omega2 = omega(r) * omega(r);
	vtherm = v_thermal(r);
	conv_fact = SUN2GR / (AU2CM * AU2CM * AU2CM);

	m_new = 0.75 * omega2 * DUST2GAS / vtherm * z * m;
	z_new = -PDENSITY / (density * conv_fact) * (PRADIUS / AU2CM) / vtherm * omega2 * z;	/*	PRADIUS: cm --> AU	*/

}

/*	ez az allando lepeskozu negyedrendu Runge-Kutta modszer	*/
void time_stepping(double r, double dt, double *y, double *ynew) {

	int i;
	double dy1[2], dy2[2], dy3[2], dy4[2];
	double ytemp[2];

	eq_rhs(r,y,dy1);

	for (i = 0; i < 2; i++) ytemp[i] = y[i] + 0.5 * dt * dy1[i];
	eq_rhs(r,ytemp,dy2);

	for (i = 0; i < 2; i++) ytemp[i] = y[i] + 0.5 * dt * dy2[i];
	eq_rhs(r,ytemp,dy3);

	for (i = 0; i < 2; i++) ytemp[i] = y[i] + dt * dy3[i];
	eq_rhs(r,ytemp,dy4);

	for (i = 0; i < 2; i++) ynew[i] = y[i] + dt * (dy1[i] + 2.0 * dy2[i] + 2.0 * dy3[i] + dy4[i]) / 6.0;

}


int main() {
  
	double time, time_integration, dt, z0, m0, r, prad3, z, m, Y[2], Y_new[2];
	FILE *fout;

	fout = fopen("por.dat","w");


	time = 0.0;
	dt = 0.05;						/*	idolepes				*/
	time_integration = 10.0 * 365.25;			/*	numerikus integralas idotartama		*/

	prad3 = PRADIUS * PRADIUS * PRADIUS;
	r = 1.0;						/*	a porszemcse csillagtol valo tavolsaga AU-ban	*/
	z0 = 5.0 * ASP_RATIO * r;					/*	a porszemcsenek a korong fosikjatol mert tavolsaga a nyomasi skalamagassag egysegeben	*/
	m0 = PDENSITY * 4.0 * M_PI / 3.0 * prad3 / SUN2GR;	/*	a reszecske tomege gramm --> M_NAP			*/
	printf("Kezdeti magassag AU-ban: %lg, kezdeti tomeg gr-ban: %lg \n", z0, m0 * SUN2GR);
	printf("nyomj egy entert a folytatashoz! \n");
	getchar();
	
	Y[0] = m0;
	Y[1] = z0;


	do {
	
		time_stepping(r, dt, Y, Y_new);			/*	Runge-Kutta modszer	*/

		Y[0] = Y_new[0];
		Y[1] = Y_new[1];

		time = time + dt;
		printf("%lg\n",time);

/*	Kiiratas a file-ba, itt az egesz idointervallum alatt 1000 adatot ir ki. Elso oszlop az ido (evben), a masodik a tomeg (grammban), a harmadik a porszemcse tavolsaga a fosiktol (a nyomasi skalamagassag egysegeben)	*/	
		if(fmod(time,(int)(time_integration/1000.0)) < dt || time == 0.0) {
			fprintf(fout,"%lg \t %lg \t %lg \n",time/365.25 ,Y[0] * SUN2GR, Y[1]);
			fflush(fout);
		}	
		
	} while (time < time_integration);  
	
	fclose(fout);

	return 0;
  
}
