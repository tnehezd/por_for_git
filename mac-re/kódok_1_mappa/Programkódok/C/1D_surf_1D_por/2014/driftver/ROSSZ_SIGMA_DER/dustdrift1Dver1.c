#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TWOPI       2.0*M_PI
#define G           0.01720209895  //Gaussian grav. constant
#define STAR        1.0
#define ALPHA		0.996
#define SDCONV      1.12521e-7
#define HASP		0.05
#define AU2CM		1.496e13
#define PDENSITY    4.0
#define SUN2GR	    1.989e33
#define NGRID       5000
#define SIGMAP_EXP  -1.5
#define SIGMA0      1700 //g/cm2

int NV, PARTICLE_NUMBER;

double RD[NGRID], SDISK[NGRID], PDISK[NGRID], DPDRDISK[NGRID], DD, RMIN, RMAX;

FILE *fin1, *fin2, *fout; 


double Stokes_Number(double pradius, double r, double sigma) {
	//in the Epstein drag regime
	return PDENSITY*pradius*M_PI/(2.0*sigma);
	
} 

double sigmagas_mmsn_cgs(double r) {
	
	
	return SIGMA0*pow(r,SIGMAP_EXP);
	
}

/*    dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs    */

void eqrhs(double pradius, double r, double *drdt) {
     
     double rmid, rindex, coef1, coef2, P, coef3, H, dPdr, Sigma, St, csound, vkep, G2;
     int index;
     
         
       
     G2 = G*G;
     
     H = r*HASP;
     
     vkep = sqrt(G2*STAR/r);                                   
	 
	 csound = vkep*HASP;   
	
	 Sigma = sigmagas_mmsn_cgs(r);
	
	 P = HASP*G2*STAR*Sigma/(sqrt(2.0*M_PI*r*r));
	 dPdr = (SIGMAP_EXP-2)*HASP*G2*STAR*SIGMA0*pow(r,(SIGMAP_EXP-3.0))/sqrt(2.0*M_PI);
	
	 St = Stokes_Number(pradius,r,Sigma);	
     
     *drdt = St/(1+St*St)*H/P*dPdr*csound;
     

}

//megiscsak vektorokkal kell majd csinalni?
void int_step(double prad, double step, double y, double *ynew)
{
	int n;
	double dy1,dy2,dy3,dy4;
	double ytemp;
	
	eqrhs(prad, y, &dy1);
	
	ytemp = y + 0.5*step*dy1;
	eqrhs(prad, ytemp, &dy2);
		
	ytemp=y+0.5*step*dy2;
	eqrhs(prad, ytemp, &dy3);
	
	ytemp=y+step*dy3;
	eqrhs(prad, ytemp, &dy4);
	
	*ynew=y+step*(dy1+2.0*dy2+2.0*dy3+dy4)/6.0;
	
}


main() {
   
   fout = fopen("pormozgas.dat","w");
   fin1 = fopen("init_data.dat","r");
//   fin2 = fopen("denspress.dat","r");
   
   int i, L, dummy;
   double distance, reprmass, r, sd, pressure, dpdr, StN, y, y_out, t, t_integration, deltat, particle_radius;
   
   PARTICLE_NUMBER = 1000;
   
   double radius[PARTICLE_NUMBER][2];
   
   
    
   /*for (i = 0; i < NGRID; i++) {
        fscanf(fin2,"%lg %lg %lg %lg",&r,&sd,&pressure,&dpdr);
        RD[i] = r;
        SDISK[i] = sd;
        PDISK[i] = pressure;
        DPDRDISK[i] = dpdr;
	   
	}*/

	
	for (i=0; i<PARTICLE_NUMBER; i++) {
            fscanf(fin1,"%d %lg %lg %lg",&dummy,&distance,&particle_radius,&reprmass);
            radius[i][0] = distance;
	    radius[i][1] = particle_radius;
	}
	

	
   t = 0.0;
   t_integration = 100000.0*365.25; //numerikus integralas idotartama
   deltat = 0.1; //idolepes napban

   L = 0;

   do {
		
		for (i=0; i<PARTICLE_NUMBER; i++) {
		
		     y   = radius[i][0];
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
   


}

