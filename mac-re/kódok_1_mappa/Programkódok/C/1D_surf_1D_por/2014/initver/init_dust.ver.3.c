//peldaprogram arra, hogyan kellene kiszamolni a reprezentativ reszecskek tomeget
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define SDCONV		1.12521e-7
#define SIGMAP_EXP     	-1.5
#define ICEFACTOR       3.0
#define DUST2GASR       0.01
#define SNOWLINE        2.7
#define ALPHA           0.01
#define STAR            1.0
#define G               0.01720209895  //Gaussian grav. constant
#define DISKMASS        0.01           //Solar mass units
#define HASP		0.05
#define AUPDAY2CMPSEC   1.7314568e8    //conversion between the velocities 

#define RMIN         	1.5
#define RMAX        	99.


FILE *fin, *fout, *fout2;


/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r) {

   	double r_dze_i = 2.7;
   	double Dr_dze_i= 0.27;
   	double r_dze_o = 15.0;
   	double Dr_dze_o= 1.5;
  	double a_mod = 0.01;
   
   	double alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
    
   	return alpha_r * ALPHA;

}


/*	sigma0 kiszámolása (illetve megadása) M_NAP / AU / AU-ban	*/
double sigma_null() {
	
//	double alpha2 = SIGMAP_EXP + 2;
	return 0.0002; //(alpha2/2.0*M_PI) *  DISKMASS/pow(RMAX,alpha2);
	
}


/*	sigma kiszámolása r helyen	*/
double sigma_gas(double r) {
	
	return sigma_null() * pow(r,SIGMAP_EXP);
	
}


/*	a por feluletisurusegenek kiszamitasa, a hohatar figyelembevetelevel	*/
double sigma_dust(double r) {

	double sigma0, sig_dust;
	
	sigma0 = sigma_null();

    	if (r <= SNOWLINE) {
	       	sig_dust = DUST2GASR * sigma0 * pow(r,SIGMAP_EXP);
    	} else {
        	sig_dust = DUST2GASR * sigma0 * ICEFACTOR * pow(r,SIGMAP_EXP);
	}

        return sig_dust;

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


int main() {  

  	int i, NGRID; 		//ez egyben a reprezentativ reszecskek szamat is jelenti
  	double r, DD, reval, reval2, reppmass;
  	double f_drift, f_frag, u_frag, u_frag2; //H
	double G2, v_kep, c_s, v_kep2, c_s2, Sigma, Sigma_cgs, P, dPdr, dlnPdlnr, rho_p, s_drift, s_frag, s_df;
  	double s_max;
  
  	fout = fopen("init_data2.dat","w");
  	fin = fopen("denspress.dat","r");

  	NGRID = 1000;			/*	a reszecskek szama	*/
    
	DD = (RMAX-RMIN)/NGRID;
   
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
  
  	f_drift = 0.55;
  	f_frag  = 0.37;
  
  	rho_p = 1.6;    	// g/cm^3
  	u_frag = 1000.0; 	// cm/s 
	//persze ez a reszecske osszeteteletol is fugg, valahogyan ezt is figyelembe lehetne/kellene venni
  
  	u_frag = u_frag / AUPDAY2CMPSEC; //u_frag CsE/nap mertekegysegben
  	u_frag2 = u_frag * u_frag;
  
  	for (i=0; i<NGRID; i++) {
       
       		r = RMIN + i * DD;
       		reval = r + DD / 2.0;
	   	reval2 = reval * reval;
      
       		reppmass = 2.0 * M_PI * r * DD * sigma_dust(reval); //ez itt a reprezentativ reszecske tomege
       
       		G2 = G * G;
//     		H = reval*HASP;
       		v_kep = sqrt(G2 * STAR/reval);  
                                       
	   	c_s = v_kep*HASP;
        
       		v_kep2 = v_kep * v_kep;
       		c_s2 = c_s * c_s;
              
	   	Sigma = sigma_gas(reval);
	   	P = HASP * G2 * STAR * Sigma / (sqrt(2.0 * M_PI) * reval2);
	   	dPdr = (SIGMAP_EXP - 2) * pow(reval,(SIGMAP_EXP - 3.0)) * HASP * G2 * STAR * sigma_null() / sqrt(2.0 * M_PI);
	   	dlnPdlnr = reval/P * dPdr;
	  
	  	Sigma_cgs = sigma_gas(reval) / SDCONV;
//	  	printf("Sigma_cgs %lg %lg\n",Sigma_cgs,sigma_null());

   		//reprezentativ reszecske kezdeti meretenek meghatarozasa
       		// 1. radialis drift altal meghatarozott maximalis meret
       		s_drift = f_drift * 2.0 / M_PI * Sigma_cgs/(100.0 * rho_p) * v_kep2 / c_s2 * fabs(1.0 / dlnPdlnr);
       
       		// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret
       		s_frag = f_frag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(reval)) * u_frag2 / c_s2;
       		//s_frag = f_frag * 2.0/(3.0*M_PI) * Sigma_cgs/(rho_p*ALPHA) * u_frag2/c_s2;

       		// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret
       		s_df = u_frag * v_kep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

       		s_max = find_min(s_drift,s_frag,s_df);

      		printf("s_drift: %lg, s_frag: %lg, s_df: %lg, smin: %lg\n", s_drift,s_frag,s_df, s_max);
       
 //      	printf("a reprezentativ reszecske sorszama, csillagtol valo tavolsaga es tomege %d %lg %lg %lg %lg %lg\n", i, reval, reppmass,s_drift,s_frag,s_df);
       		fprintf(fout,"%d %lg %lg %lg %lg %lg %lg\n",i,reval,reppmass,s_max,s_drift,s_frag,s_df);
       		//getchar();
  
  	}


  	printf("A reszecskeket tartalmazo file elkeszult (init_data.dat), a korong parametereit tartalmazo file kiirasa! \n\n A file tartalma: \n i (a reszecske sorszama), r (a reszecske tavolsaga), prad (a reszecske merete cm-ben), reppmass (a reszecske reprezentativ tomege)\n");
  	printf("A folytatashoz nyomj ENTER-t!\n");
  	getchar();


  	fout2 = fopen("disk_param.dat","w"); 
       	fprintf(fout2,"%lg %lg %d %lg %lg\n",RMIN, RMAX, NGRID, SIGMAP_EXP, sigma_null() * SDCONV);

  	fclose(fout2);   

  	printf("A korong parametereit tartalmazo file elkeszult (disk_param.dat) \n\n A file tartalma: \n RMIN, RMAX, NGRID, a profil kitevoje (SIGMAP_EXP), és sigma0 (M_Sun / AU / AU)\n\n");    


  	return 0;

}
