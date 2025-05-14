#include<stdio.h>
#include<math.h>
#include<stdlib.h>


#define M_star 		1.0				/*	classical T Tauri mass				*/
#define G 		2.96002e-04 			/*	Newton-grav const in AU^3 / M_Sun / day^2	*/
#define k		sqrt(G * M_star)		/*	Gaussian graviation constant			*/
#define ngrid 		1000				/*	a gridpontok szama				*/
#define rmin 		0.1				/*	0.1 AU-tol szamol				*/
#define rmax 		10.0				/*	10 AU-ig szamol					*/
#define dgrid 		(rmax-rmin)/(ngrid-1)		/*	dr 						*/
#define mu		(k * k * M_star)
#define timemax		700

FILE *der, *filgnu, *filplot;

double kepler_vel(double r) {

	return sqrt(G * M_star / r);		/*	ez a körsebesség 	*/

}


double r_calc(double *rvec) {			/*	az adott koordináták mellett az r távolság	*/
	
	double r;
	r = sqrt(rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]);
	return r;	

}


void init_pos(double *rvec) {				/*	kezdeti pozicio		*/

	rvec[0] = 1.0;
	rvec[1] = 0.0;
	rvec[2] = 0.0;

}


void init_vel(double *velvec, double *rvec) {		/*	kezdeti sebesség	*/

	double r;

	init_pos(rvec);

	r = r_calc(rvec);
	velvec[0] = 0.0;
	velvec[1] = kepler_vel(r);
	velvec[2] = 0.0;

}


double calc_per(double *rvec) {				/*	a periodus kiszamolasa		*/

	double p, r;

	init_pos(rvec);

	r = r_calc(rvec);
	p = sqrt(r * r * r * 4.0 * M_PI * M_PI / (G * M_star));		/*	Kepler III.	*/
	return p;

}


double calc_acc(double coord, double r) {				/*	gyorsulas kiszamolasa	*/

	return (-G * M_star * coord / (r * r * r));
}


void calc_coeff1(double time, double *velvec, double *coeff) {		/*	RK4 az elso diffegyenlethez a koefficiensek kiszamolasa	*/

	coeff[0] = velvec[0];
	coeff[1] = velvec[1];
	coeff[2] = velvec[2];
	
}

void calc_coeff2(double time, double *rvec, double *coeff) {		/*	RK4 a masodik diffegyenlethez a koefficiensek kiszamolasa	*/

	double r;
	r = r_calc(rvec);

	coeff[0] = calc_acc(rvec[0],r);
	coeff[1] = calc_acc(rvec[1],r);
	coeff[2] = calc_acc(rvec[2],r);
	
}


void Gnuplot(const char * gnucommand){				/*	gnuplot parancs meghivasa	*/
  
	char syscommand[1024];

	printf("... A gnuplot megnyitasa es az animaciok elkeszitese ...\n\n\n");
 
	sprintf(syscommand, "echo \"%s\" | gnuplot", gnucommand);
	system(syscommand);

}


int main(){
  
  	double rvec[3], velvec[3], rvecnew[3], velvecnew[3], accvec[3], x, y, z, vx, vy, vz, per, dt, time, r;
	double an1[3], bn1[3], cn1[3], dn1[3], an2[3], bn2[3], cn2[3], dn2[3];
	int i,fcount;
	char der_name[64],animplot[50], gifnev[50], timechar[50], filekinev2[100];

	time = 0.0; 
	fcount = 0;	

	sprintf(animplot,"anim.gnuplot");
	sprintf(timechar,"t=%%i");
	sprintf(filekinev2, "der_%%i.dat");
	sprintf(gifnev,"rungekutta.gif");

	 
	init_pos(rvec);
	init_vel(velvec,rvec);
	per = calc_per(rvec);
	dt = calc_per(rvec) / 100.0;

	x = rvec[0];
	y = rvec[1];
	z = rvec[2];

	vx = velvec[0];
	vy = velvec[1];
	vz = velvec[2];


	for (i = 0; i < 3; i++) {
		printf("%lg\n",rvec[i]);
		accvec[i] = calc_acc(rvec[i],r_calc(rvec));
	}

//	printf("%lg %lg %lg\n", accvec[0],accvec[1],accvec[2]);

	do {	

		for (i = 0; i < 3; i++) {
	
			calc_coeff1(time,velvec,an1);			/*	elso koefficiensek kiszamolosa	*/
			calc_coeff2(time,rvec,an2);

			rvecnew[i] = rvec[i] + dt / 2.0 * an1[i];
			velvecnew[i] = velvec[i] + dt / 2.0 * an2[i];

			calc_coeff1(time + dt / 2.0,velvecnew,bn1);	/*	masodik koefficiensek kiszamolsa	*/
			calc_coeff2(time + dt / 2.0,rvecnew,bn2);


			rvecnew[i] = rvec[i] + dt / 2.0 * bn1[i];
			velvecnew[i] = velvec[i] + dt / 2.0 * bn2[i];

			calc_coeff1(time + dt / 2.0,velvecnew,cn1);	/*	harmadik koefficiensek kiszamolsa	*/
			calc_coeff2(time + dt / 2.0,rvecnew,cn2);

			rvecnew[i] = rvec[i] + dt * cn1[i];
			velvecnew[i] = velvec[i] + dt * cn2[i];

			calc_coeff1(time + dt,velvecnew,dn1);		/*	negyedik koefficiensek kiszamolsa	*/
			calc_coeff2(time + dt,rvecnew,dn2);
			rvec[i] = rvec[i] + dt / 6.0 * (an1[i] + 2.0 * bn1[i] + 2.0 * cn1[i] + dn1[i]);
			velvec[i] = velvec[i] + dt / 6.0 * (an2[i] + 2.0 * bn2[i] + 2.0 * cn2[i] + dn2[i]);
		
		}

		if(fmod(time, (int)(timemax/200.0)) < dt || time == 0) {

			snprintf(der_name,64,"der_%i.dat",fcount);
			der=fopen(der_name,"w");
			fprintf(der,"%lg \t %lg \t %lg \t %lg \t %lg \t %lg \n",rvec[0],rvec[1],rvec[2], velvec[0], velvec[1], velvec[2]);
			filgnu=fopen(animplot,"w");
			fprintf(filgnu,"filename(n) = sprintf('%s', n) \n set xlabel 'Distance [AU]' \n set xrange[-2:2] \n set yrange[-2:2] \n set ylabel 'Distance [AU]' \n set title 'Central force problem' \n plot filename(i) using 1:2 title sprintf('%s * 3.652 nap',i) with points pt 7 ps 0.5  \n i=i+1 \n if (i < n) reread \n ",filekinev2,timechar); 

			fclose(filgnu);

			filplot=fopen("plot.rungekutta","w");

/*	Animacio letrehozasa az elkeszult file-okbok. A szamozott file-okat egy for ciklussal nyitja meg a gnuplot.	*/		
			fprintf(filplot,"set term gif giant animate size 450, 450 \n set output '%s' \n n=192 \n i=0 \n load '%s' \n", gifnev, animplot); 

			fclose(filplot);
			
		}
		fclose(der);	

	fcount++;						/*	az adatok file-ba valo kiiratasahoz a futoindex	*/
	time = time + dt;
	} while (time < timemax);


	printf("%lg %lg %lg %lg %lg %lg %lg\n",x,y,z,vx,vy,vz,per);
  

	Gnuplot("load 'plot.rungekutta'");	


	return 0;
  
}


