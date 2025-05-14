//peldaprogram arra, hogyan kellene kiszamolni a reprezentativ reszecskek tomeget
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define SDCONV		1.12521e-7
#define SIGMAP_EXP  	-1.0
#define ICEFACTOR    	3.0
#define DUST2GASR    	0.01
#define SNOWLINE     	2.7
#define RHOP         	4.0 		//density of the dust particle in g/cm3

#define DISKMASS     	0.01

#define RMIN         	1.5
#define RMAX        	99.

FILE *fout, *fout2;

double my_rand() {
	return (double) rand()/RAND_MAX;
}


double sigma_null() {

       double alpha2 = SIGMAP_EXP + 2;
       return (alpha2/2.0*M_PI)*DISKMASS/(pow(RMAX,alpha2) - pow(RMIN,alpha2));

}


//ez adja meg a por feluleti suruseget MMSN-ben (Hayashi 1980 nyoman)
//ide mas fuggvenyt is lehet tenni
//a lenyeg, hogy ebben az esetben 2.7 AU-nal, ahol a hohatar van, a
//feluleti suruseg ugrik egyet

double sigmap_mmsn(double r) {

	double sigma0, sigmap;
	
	if (r <= 2.7) {
		sigma0 = 10.0;
		sigmap = sigma0*pow(r,SIGMAP_EXP)*SDCONV;	
	} else {
		sigma0 = 30.0;  //itt eredetileg ennel is nagyobb ugras van, az eredeti szorzotenyezo 4
		sigmap = sigma0*pow(r,SIGMAP_EXP)*SDCONV;
    	}

	return sigmap;
}

double sigmagas_mmsn_cgs(double r) {
	
	double sigma0 = 100.;
	return sigma0*pow(r,SIGMAP_EXP);
	
	
}


	
double sigmap_general(double r) {


	double sigma0, sig_dust;
	
	sigma0 = sigma_null();

    	if (r <= SNOWLINE) {
	       	sig_dust = DUST2GASR*sigma0*pow(r,SIGMAP_EXP);
    	} else {
        	sig_dust = DUST2GASR*sigma0*ICEFACTOR*pow(r,SIGMAP_EXP);
	}

        return sig_dust;
}



int main() {  

  	int i, NGRID; //ez egyben a reprezentativ reszecskek szamat is jelenti
  	double rmin, rmax, r, delta_r, reval, reppmass, pradius_cm;

  	fout = fopen("init_data.dat","w");

  	NGRID = 1000;
  	rmin = RMIN;   //ez itt hulyen van megirva, atirhatod, ha szeretned :-)
  	rmax = RMAX;
  
  	delta_r = (rmax - rmin)/((double) NGRID);
  
  
	printf("a parameterek osszegzese\n");
  	printf("a korong tomege Naptomegben %lg\n",DISKMASS);
  	printf("korong belso hatara %lg CsE kulso hatara %lg CsE\n",RMIN,RMAX);
  	printf("a feluleti suruseg hatvanyfuggveny alakjanak kitevoje %lg\n",SIGMAP_EXP);
  	printf("a hohatar helyzete %lg CsE \n", SNOWLINE);
  	printf("a por feluleti suruseg ugrasanak tenyezoje a hohataron tul %lg\n",ICEFACTOR);
  	printf("a gaz feluleti surusege a csillagtol 1 CsE tavolsagban %lg\n",sigma_null());
  	printf("a por gaz arany %lg\n", DUST2GASR);
  	printf("a reprezentativ reszecskek szama %d\n",NGRID);
  	printf("nyomj ENTER-t a folytatashoz!\n");
  	getchar();
  
	for (i=0; i<NGRID; i++) {

		r = rmin + i*delta_r;
       
       		reval = r + delta_r/2.0;
       		reppmass = 2.0*M_PI*r*delta_r*sigmap_general(reval); 
       		pradius_cm = 2.0*sigmagas_mmsn_cgs(reval)/(M_PI*RHOP);
       
//       	printf("a reprezentativ reszecske sorszama, csillagtol valo tavolsaga es tomege %d %lg %lg\n",i, reval, reppmass);
       		fprintf(fout,"%d %lg %lg %lg\n",i,reval,pradius_cm,reppmass);
//		getchar();
  
  	}
	
	fclose(fout);

  	printf("A reszecskeket tartalmazo file elkeszult (init_data.dat), a korong parametereit tartalmazo file kiirasa! \n\n A file tartalma: \n i (a reszecske sorszama), r (a reszecske tavolsaga), prad (a reszecske merete cm-ben), reppmass (a reszecske reprezentativ tomege)\n");
  	printf("A folytatashoz nyomj ENTER-t!\n");
  	getchar();


  	fout2 = fopen("disk_param.dat","w"); 
       	fprintf(fout2,"%lg %lg %d %lg %lg\n",RMIN, RMAX, NGRID, SIGMAP_EXP, 100.*SDCONV);

  	fclose(fout2);   

  	printf("A korong parametereit tartalmazo file elkeszult (disk_param.dat) \n\n A file tartalma: \n RMIN, RMAX, NGRID, a profil kitevoje (SIGMAP_EXP), és sigma0 (M_Sun / AU / AU)\n\n");    

	return 0;

}
