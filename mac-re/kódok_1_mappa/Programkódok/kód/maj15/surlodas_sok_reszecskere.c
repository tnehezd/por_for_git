#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define G 		0.01720209895			/*	Gaussian-grav const in astron. units		*/
#define G2		G*G
#define GCGS		6.674e-8			/*	Gaussian-grav const in CGS			*/
#define	GCGS2		GCGS*GCGS
#define AU2CM		1.496e13			/*	1 AU --> 1 cm					*/
#define asp_ratio 	0.05				/*	aspect ratio					*/

#define MMOL		3.9e-24
#define SIGMOL		2.0e-15
#define PDENSITY    	3.0
#define CBODY		1.0				/*	central body mass in solar units		*/
#define GCM2MSAU	1.125211e-7			/*	g/cm2 --> M_Sun/AU2				*/
#define	DDENSASTR2CGS	7.95890e-17			/*	deriv p: M_Sun/day/day/AU/AU/AU --> g/s/s/m/m/m	*/
#define M_SOL		1.9891e33			/*	mass of Sun in gramms				*/

#define filenev		"dens_559000.dat"
#define filenev2	"distr.dat"


FILE *filebe, *fildrag, *outcgs;

int Sorokszama(int numout){

	char c;
	filebe = fopen(filenev,"r+");		
	numout = 0;

/*	A feluletisuruseg profilt tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal	*/

		while((c = fgetc(filebe)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 					*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert								*/
	fclose(filebe);	
	return numout;

}


int Reszecskekszama(int numout){

	char c;
	filebe = fopen(filenev2,"r+");		
	numout = 0;

/*	A random generalt reszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal	*/

		while((c = fgetc(filebe)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 						*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert									*/
	fclose(filebe);				
	return numout;

}


/*	A feluletisuruseg profilt tartalmazo file tartalmanak beolvasasa: r a sugar, sigma a feluletisuruseg, dp a nyomas derivaltja	*/
void load_r_sigma_dp(int num, double *rvec, double *sigmavec, double *dpressvec, double *DD, double *RMIN, double *RMAX) {
		
	int i;
	double temp[num];					/*	homerseklet tarolasa, egyelore nem hasznaljuk	*/

	filebe = fopen(filenev,"r+");

	for (i = 1; i <= num; i++){ 
/*	file-bol beolvasas: r,sigma,homerseklet,dp	*/
		fscanf(filebe, "%lf %lf %lf %lf",&rvec[i],&sigmavec[i],&temp[i],&dpressvec[i]); 
//		printf("%i %lg %lg %lg\n",i,rvec[i],sigmavec[i],dpressvec[i]);
	}
	fclose(filebe);

	rvec[0] = rvec[1];
	rvec[num] = rvec[num-1];

	*DD = rvec[1] - rvec[0];		/*	lepeskoz	*/
	*RMIN = rvec[0];			/*	"rmin"		*/
	*RMAX = rvec[num-1];			/*	"rmax"		*/
}


/*	a reszecske adatainak betoltese: rad a sugara (cm-ben!!!), x,y,z a hely-, vx,vy,vz a sebessegvektor koordinatai (AU es AU/day)	*/
void load_rad_xyz_vxvyvz(int partnum, double *rad, double *x, double *y, double *z, double *vx, double *vy, double *vz) {
		
	int i;

	filebe = fopen(filenev2,"r+");

	for (i = 0;i < partnum; i++){ 
/*	file-bol beolvasas: r,sigma,homerseklet,dp	*/
		fscanf(filebe, "%lf %lf %lf %lf %lf %lf %lf",&rad[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]); 
	}
	fclose(filebe);
}


/*	a beoltott adatok egimechanikai mertekegysegekbol cgs-be valo atvaltasa (surusegprofilra)	*/
void r_sigma_dp_to_cgs(double r, double sigma, double dpress, double DD, double RMIN, double RMAX, double *rvec_cgs, double *sigmavec_cgs, double *dpressvec_cgs, double *DD_cgs, double *RMIN_cgs, double *RMAX_cgs) {

	rvec_cgs = r * AU2CM;				/*	AU --> cm					*/	
	sigmavec_cgs = sigma / GCM2MSAU;		/*	M_Sun/AU2 --> g/cm2				*/
	dpressvec_cgs = dpress * DDENSASTR2CGS;	/*	M_Sun/Day/Day/AU/AU/AU --> g/s/s/cm/cm/cm	*/

	*DD_cgs = DD * AU2CM;						/*	AU --> cm	*/
	*RMIN_cgs = RMIN * AU2CM;					/*	AU --> cm	*/
	*RMAX_cgs = RMAX * AU2CM;					/*	AU --> cm	*/

}


/*	a beoltott adatok egimechanikai mertekegysegekbol cgs-be valo atvaltasa (reszecskekre)	*/
void xyz_vxvyvz_to_cgs(double x, double y, double z, double vx, double vy, double vz, double *x_cgs, double *y_cgs, double *z_cgs, double *vx_cgs, double *vy_cgs, double *vz_cgs) {

	x_cgs = x * AU2CM;				/*	AU --> cm	*/	
	y_cgs = y * AU2CM;				/*	AU --> cm	*/
	z_cgs = z * AU2CM;				/*	AU --> cm	*/
	vx_cgs = vx * AU2CM / 86400.0;			/*	AU/day --> cm/s	*/	
	vy_cgs = vy * AU2CM / 86400.0;			/*	AU/day --> cm/s	*/
	vz_cgs = vz * AU2CM / 86400.0;			/*	AU/day --> cm/s	*/

}


double rho_midplane(double localsigmacgs, double rcgs) {

	return (localsigmacgs / (sqrt(2.0 * M_PI) * asp_ratio * rcgs));              		/*	this is rho in the midplane of the disk in cgs	*/

}


double rho_gas(double localsigmacgs, double r, double z) {

	double rho0;
	rho0 = rho_midplane(localsigmacgs,r);
	return (rho0 * exp(-z * z / (2.0 * asp_ratio * asp_ratio * r * r)));           		/*	in g/cm^3 -- itt bevezettuk a z-tol valo fuggest, ugyanis a por-reszecske elhagyhatja a korongot	*/

}


/*	mean free-path in cm	*/
double lambda(double localsigma, double r, double z) {

	return (MMOL / SIGMOL) / rho_gas(localsigma,r,z);

}



double Knudsen(double body_radius, double localsigma, double r, double z) {			/*	Knudsen number, if lamda and body_radius are in cm --> dimensionless		*/

	return (0.5 * lambda(localsigma,r,z) / body_radius);

}


double v_kepler_cgs(double r) {

	double vkep;
	vkep = sqrt(GCGS2 * CBODY * M_SOL / (r));							/*	v_kepler in AU/day	*/
	return vkep;											/*	v_kepler in cm/s	*/

}





double speed_of_sound(double r) {
	
	return v_kepler_cgs(r) * asp_ratio;
	
}


double Mach_number(double r, double delv) {							/*	dimensionless, both delv and csound are in AU/day		*/

	return (delv / speed_of_sound(r));

}


double Reynolds_number(double r, double delv, double body_radius, double localsigma, double z) {	/*	dimensionless, both body_radius and lambda are in cm	*/

	return (6.0 * (body_radius) * delv / (sqrt(8.0 / M_PI) * speed_of_sound(r) * lambda(localsigma,r,z)));

}


double Epstein_coeff(double r, double delv) {

	return (2.0 * sqrt(1.0 + 128.0 / (9.0 * M_PI * Mach_number(r,delv) * Mach_number(r,delv))));

}


double Stokes_coeff(double body_radius, double delv, double r, double localsigma, double z) {

	double cstk,ren;
	
	ren = Reynolds_number(r,delv,body_radius,localsigma,z);  

     	if (ren <= 500.0) {
         	cstk = 24.0 / ren + 3.6 * pow(ren,-0.313);     
     	} else 

   
     	if (ren <= 1500.0) {
         	cstk = 9.5e-5 * pow(ren,1.397);       
     	} else {
         	cstk = 2.61;
     	}

	return cstk;	

}
   

double Drag_Coeff(double body_radius, double delv, double r, double localsigma, double z) {         

	double knud,ceps,cstk;
	knud = Knudsen(body_radius,localsigma,r,z);
	ceps = Epstein_coeff(r,delv);
	cstk = Stokes_coeff(body_radius,delv,r,localsigma,z);
	return (9.0 * knud * knud * ceps + cstk) / ((3.0 * knud + 1.0) * (3.0 * knud + 1.0));
}


void interpol(double RMIN, double DD, double r, double *sigmavec, double *dpressvec, double *sigma_interp, double *dpress_interp) {

	double coef1, coef2, rmid;
	int index, rindex;

    	rmid = (r - RMIN)/DD;     						/* 	the integer part of this gives at which index is the body	*/
   	index = (int) floor(rmid);						/* 	ez az rmid egesz resze						*/
    	rindex = rvec[index];       						/* 	the corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

   	coef1 = (sigmavec[index + 1] - sigmavec[index]) / DD;      		/*	ez az alabbi ket sor a linearis interpolacio - remelem, jo!!!	*/
     	*sigma_interp = sigmavec[index] + coef1 * (r - rindex);          	/*	g/cm^2	*/

     	coef2 = (dpressvec[index + 1] - dpressvec[index]) / DD; 		/*	ugyanaz az interpolacio a nyomas derivaltjara			*/
     	*dpress_interp  = dpressvec[index] + coef2 * (r - rindex);

}



/*Olyan fuggveny kell, mely kiszamolja a gaz surlodasi erot. Ehhez kell a gaz sebessege, a reszecske sebesseget szamolja a kettest-problemas program, a ketto kulonbsege pedig a reszecske gazhoz viszonyitott relativ sebessege. 

A gaz sebessegehez kell nekunk a nyomas tavolsag szerinti derivaltja a jolismert egyenlet szerint:

(lasd papir )
*/


void get_drag_acceleration(double body_radius, double x, double y, double z, double vx, double vy, double vz, double *adrag, double *dvx, double *dvy, double *dvz, double sigma_i, double sigma i_1, double dpress_i, double dpress_i_1, double DD, double RMIN)  {
/*	body_radius: a reszecske sugara; body_index:ez a reszecskek tomeget/sugarat tartalmazo tombbeli index; x,y,z,vx,vy,vz: a reszecske hely es sebesseg koordinatai a csillaghoz viszonyitva; adrag: a gaz surlodasbol szarmazo gyorsulas; dvx,dvy,dvz: a reszecskenek a gazhoz viszonyitott sebessege	*/
                           
   	double r, r2, r_cgs, r2_cgs, sigma_cgs, dpress_cgs, rho0_cgs, rhogas_cgs, ugas, ux, uy, uz, ldvx, ldvy, ldvz, delv, vkep,rhogas,cdrag,sigma_interp,dpress_interp;                    
     	int index;
   
     	r2 = x * x + y * y + z * z;				/*	r^2	*/
 	r = sqrt(r2);

	r_sigma_dp_to_cgs(r,sigma_i,dpress_i,DD,RMIN,RMAX,&r_cgs,&sigma_cgs,&dpress_cgs,*DD_cgs,*RMIN_cgs,*RMAX_cgs);
	xyz_vxvyvz_to_cgs(x,y,z,vx,vy,vz,x_cgs,y_cgs,z_cgs,vx_cgs,vy_cgs,vz_cgs);

     	r2_cgs = x_cgs * x_cgs + y_cgs * y_cgs + z_cgs * z_cgs;				/*	r^2	*/
 	r_cgs = sqrt(r2_cgs);

	interpol(RMIN,DD,r,sigmavec,dpressvec,&sigma_interp,&dpress_interp,);
	
	rhogas = rho_gas(sigma_interp,r,z);

/*	a gaz sebessege a nyomas jelenleteben a rendszerben	*/
     	ugas = sqrt(v_kepler_cgs(r) * v_kepler_cgs(r) + r * dpress_interp / rhogas);  		/* 	second term in cm^2/sec^2, ugas in cm/sec			*/
     
   	ugas = ugas*(86400.0/AU2CM);             		/* 	conversion back to AU/day units					*/
     
 	ux =-ugas*y/r; // a gaz sebessege 
     	uy = ugas*x/r;
     	uz = 0.0;
       
     	ldvx  = vx - ux; // a reszecske es a gaz relativ sebessege 
     	ldvy  = vy - uy;
     	ldvz  = vz - uz;
     
     	*dvx = ldvx;     // ez kell a mozgasegyenletekhez, ha mar kiszamoltuk...
     	*dvy = ldvy;
     	*dvz = ldvz;
            
     	delv = sqrt(ldvx*ldvx + ldvy*ldvy + ldvz*ldvz);
	cdrag = Drag_Coeff(body_radius,delv,r,sigma_interp,z);
   	*adrag = -3.0 * rhogas * cdrag * delv / (8.0 * PDENSITY * body_radius); //RHO bulk density of the planetesimal population   

} 

/*

 for(i=1;i<=imax;i++)
  {
   j = i-1;
   r=rmin+(i-1)*dr;
   Sigma[j]= u[i]/Func_nu( r )/SDCONV;
   Density[j] = 1.0/(sqrt(2.0*M_PI))*Sigma[j]/(Asp_Ratio*r*AU2CM);
   csound = Asp_Ratio*k_gauss*sqrt(CBODY/r);
   csound_cgs = csound*AU2CM/86400.0;
   Pressure[j] = Density[j]*csound_cgs*csound_cgs;
  }

*/


int main() {

	double t, deltat, t_integration, DD, RMIN, RMAX, DD_cgs, RMIN_cgs, RMAX_cgs,rho,adrag_i,dvx_i,dvy_i,dvz_i;
	int i,num,linesout,count,body_index,partnum;
	linesout = 0;
	num = 0;
	count = 0;
	partnum = 0;

	num = Sorokszama(linesout);		/*	A feluletisuruseg profilt tartalmazo file sorainak szama elmentve egy integerbe	*/
	linesout = 0;
	partnum = Reszecskekszama(linesout); 	/*	A resezcskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe	*/
	
	double body_radius[partnum],x[partnum],y[partnum],z[partnum],vx[partnum],vy[partnum],vz[partnum];	/*	a reszecskek adatait tarolo vektorok: body_radius a reszecske sugara, x,y,z a reszecsek hely-, vx,vy,vz a reszecskek sebessegvektorainak koordinatai	*/
	double sigmavec[num+2],rvec[num+2],dpressvec[num+2];				/*	sigmavec: feluletisuruseget, rvec: az ehhez tartalmazo x pontokat, dpressvec: a nyomas derivaltat tartalmazo vektorok	*/
	double adrag[partnum], dvx[partnum], dvy[partnum], dvz[partnum];	/*	a reszecskekre hato surlodasbol szarmazo gyorsulas, es a relativ sebessegkomponenseket tarolo vektorok	*/	
	double sigma_i, sigma_i_1, dpress_i, dpress_i_1;

	load_r_sigma_dp(num,rvec,sigmavec,dpressvec,&DD,&RMIN,&RMAX);		/*	a feluletisuruseg es a nyomasderivalt profil betoltese es eltarolasa	*/
	load_rad_xyz_vxvyvz(partnum,body_radius,x,y,z,vx,vy,vz);		/*	a porreszecskek adatainak betoltese es eltarolasa	*/

	sigma

	for(i = 0; i < partnum; i++) {
		sigma_i = sigmavec[i];
		sigma_i_1 = sigmavec[i+1];
		dpress_i = dpress[i];
		dpress_i_1 = dpress[i+1];
		get_drag_acceleration(body_radius[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],&adrag_i,&dvx_i,&dvy_i,&dvz_i,sigma_i,sigma_i_1,dpress_i,dpress_i_1,DD,RMIN);
		adrag[i] = adrag_i;
		printf("%lg %lg %lg\n",adrag_i,dvx_i,dvy_i,dvz_i);
		dvx[i] = dvx_i;
		dvy[i] = dvy_i;
		dvz[i] = dvz_i;
	}

	return 0;

}
