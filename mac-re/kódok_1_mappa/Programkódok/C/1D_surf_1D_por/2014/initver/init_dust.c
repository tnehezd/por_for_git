//peldaprogram arra, hogyan kellene kiszamolni a reprezentativ reszecskek tomeget
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define SDCONV		 1.12521e-7
#define SIGMAP_EXP  -1.5
#define ICEFACTOR    3.0
#define DUST2GASR    0.01
#define SNOWLINE     2.7
#define RHOP         4.0 //density of the dust particle in g/cm3

#define DISKMASS     0.01

#define RMIN         5.0
#define RMAX        100.0

FILE *fout;

double my_rand()
{
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
	double sigma0;
	
	if (r <= 2.7) {
	    sigma0 = 10.0;
	    return sigma0*pow(r,SIGMAP_EXP)*SDCONV;	
	}
    if (r >	2.7) {
        sigma0 = 30.0;  //itt eredetileg ennel is nagyobb ugras van, az eredeti szorzotenyezo 4
        return sigma0*pow(r,SIGMAP_EXP)*SDCONV;
    }
}

double sigmagas_mmsn_cgs(double r) {
	
	double sigma0 = 1700.0;
	return sigma0*pow(r,SIGMAP_EXP);
	
	
}


	
double sigmap_general(double r) {
    double sigma0 = sigma_null();
    if (r <= SNOWLINE) {
        return DUST2GASR*sigma0*pow(r,SIGMAP_EXP);
    }
    
    if (r > SNOWLINE) {
        sigma0 = sigma0*ICEFACTOR;
        return DUST2GASR*sigma0*pow(r,SIGMAP_EXP);
    }
    
}



main() {  

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
       
       printf("a reprezentativ reszecske sorszama, csillagtol valo tavolsaga es tomege %d %lg %lg\n",i, reval, reppmass);
       fprintf(fout,"%d %lg %lg %lg\n",i,reval,pradius_cm,reppmass);
       //getchar();
  
  }

}
