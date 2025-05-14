
// A derivs-be be kell adni egy csomo uj parametert!! pl. dt, r, stb. !!!
// Hogy legyen ez? Mas rk legyen? Mas derivs? Sajat?


#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>

#define rmin 5.0			/* ... r-tartomany ettol ... */
#define rmax 200.0			/* ... eddig tart, es    ... */
#define imax 2500			/* ... ennyi racspontra van felosztva (vigyazat, rmin-dr maradjon pozitiv!) ... */	
#define G    0.01720209895		//  Gaussian constant of gravity
#define Asp_Ratio 0.05			    //  h = H(r)/r
#define ALPHA 0.01			        //  alpha viscosity parameter
#define M_DISK 0.01
#define SIGMA_0 1.45316e-5
#define R_S    10.0
#define SDCONV 1.12521e-7
#define AU2CM  1.496e13
#define Central_Mass 1.0
#define temp_exp 1.0
#define Planet_Mass 6.0e-5		// 3.0e-6 is ~ the Earth mass in Solar Mass units
#define Mass_Par  Planet_Mass/Central_Mass
#define dr (rmax-rmin)/(imax-1)		/* terbeli felosztas lepeskoze */
#define RHO_P  1.6
#define F_FRAG 0.37
#define U_FRAG 1000.0 //10 m/s
#define AUPDAY2CMPSEC   1.7314568e8

FILE *f_sigma, *f_pos;



double Initial_rp()
{
 return(30.0);
}


double Initial_sigma(int i)
{
  double r;
  r=rmin+(i-1)*dr;
  //u[i] = M_DISK*exp(-r/R_S)/(2.0*M_PI*R_S*r); 
  return(0.0002*(pow(r,-1.5)));
}


double Func_nu(double r)
{
	double nu_r, r_dze, Dr_dze, alpha_1, alpha_2, alpha_r;
	//itt adjuk meg ezeket az ertekeket ezen eljarason belul. Nem elegans, de ez van :-) 
	
	double r_dze_i = 2.7;
	double Dr_dze_i= 0.27;
	double r_dze_o = 25.0;
	double Dr_dze_o= 2.5;
	double a_mod = 0.01;
	
	alpha_r = 1.0 - 0.5*(1.0 - a_mod)*(tanh((r-r_dze_i)/Dr_dze_i)+tanh((r_dze_o - r)/Dr_dze_o));
	
	nu_r = ALPHA*alpha_r*Asp_Ratio*Asp_Ratio*G*sqrt(Central_Mass*r);
	
	return nu_r;	
}

double alpha_turb(double r) {
	
	double r_dze_i = 2.7;
	double Dr_dze_i= 0.27;
	double r_dze_o = 25.0;
	double Dr_dze_o= 2.5;
	double a_mod = 0.01;
	
	double alpha_r = 1.0 - 0.5*(1.0 - a_mod)*(tanh((r-r_dze_i)/Dr_dze_i)+tanh((r_dze_o - r)/Dr_dze_o));
	
	return alpha_r * ALPHA;
}



double Func_A(double r)
{
 return(3.0*Func_nu(r));
}


double Func_B(double r)
{
 return(9.0*Func_nu(r)/(2.0*r));
}

double Func_C(double r)
{
 return(0.0);
}


void Parabola(double *u,int i1,int i2,int i3,double *a,double *b,double *c)
{
 double x1,x2,x3,y1,y2,y3;
 double av,bv,cv;
 
 x1=rmin+(i1-1)*dr;
 x2=rmin+(i2-1)*dr;
 x3=rmin+(i3-1)*dr;
 y1=u[i1];
 y2=u[i2];
 y3=u[i3];

 av=((y1-y3)/(x1-x3)-(y1-y2)/(x1-x2))/(x3-x2);
 bv=(y1-y2)/(x1-x2)-av*(x1+x2);
 cv=y1-av*x1*x1-bv*x1;
 
 *a=av; *b=bv; *c=cv;
}


void Perem(double *u)		/* Kifolyas -- ez ad erteket a ghost cellaknak */
{
 double a,b,c;
 //eredeti parabolas extrapolalos
 // bal oldal
 Parabola(u,1,2,3,&a,&b,&c);
 u[0]=a*(rmin-dr)*(rmin-dr)+b*(rmin-dr)+c;

 // jobb oldal
 Parabola(u,imax-2,imax-1,imax,&a,&b,&c);
 u[imax+1]=a*(rmax+dr)*(rmax+dr)+b*(rmax+dr)+c;

 //Rezso-fele:
 //u[0] = (u[1] / Func_nu(rmin))*Func_nu(rmin-dr);
 //u[imax+1] = (u[imax] / Func_nu(rmin+(imax)*dr)) * Func_nu(rmin+(imax+1)*dr);
}


double Derivs_rp(double r, double p_radius, double Sigma, double P, double dPdr)
{
 double G2, H, vkep, csound, Sigma_cgs, St;

 
	G2 = G*G;
	
	H = r*Asp_Ratio;
	
	vkep = sqrt(G2*Central_Mass/r);                                   
	
	csound = vkep*Asp_Ratio; 
	Sigma_cgs = Sigma/SDCONV;
	//Stokes szam
	St = RHO_P*p_radius*M_PI/(2.0*Sigma_cgs);
    return St/(1+St*St)*H/P*dPdr*csound;
}


void Compute_Sigma_P_dP(double time,double step,double position,double *sigma,double *Sigma,double *Pressure, double *dPdr) {
	int i,index;
	double r;
	double localtime,localstep;
	double u[imax+2], pressure[imax+2];
	double small, sigma_r, pressure_r, dsigma_left, dsigma_right, DSigma, ize;
	double dpressure_left, dpressure_right, DPressure; 
	
	localtime=time;
	localstep=step;
	
	for(i=0;i<=imax+1;i++) 
		u[i]=sigma[i]*Func_nu(rmin+(i-1)*dr);
	
	while(localtime<time+step) {
		
		for(i=1;i<=imax;i++) {
			r=rmin+(i-1)*dr;
			u[i]=u[i]+localstep*(Func_A(r)*(u[i+1]-2*u[i]+u[i-1])/(1.0*dr*dr)+Func_B(r)*(u[i+1]-u[i-1])/(2.0*dr)+Func_C(r));	// EULER STEP !!!
		}
		
		Perem(u);
		
		if((localtime<localtime+step)&&(localtime+localstep>localtime+step)) 
			localstep=time+step-localtime;
		
		localtime=localtime+localstep;
	}
	
	for(i=0;i<=imax+1;i++) {
		r=rmin+(i-1)*dr;
		sigma[i]= u[i]/Func_nu(r);
		pressure[i] = Asp_Ratio*G*G*Central_Mass*sigma[i]/(sqrt(2.0*M_PI)*r*r);
	}
	
		
	index=(int)(floor((position-rmin)/(1.0*dr))+1);
	small=position-rmin-(index-1)*(1.0*dr);
	sigma_r=(1-small)*sigma[index]+small*sigma[index+1];
	pressure_r=(1-small)*pressure[index]+small*pressure[index+1];
	*Sigma=sigma_r;
	*Pressure=pressure_r;
 	
	//na es itt kellene kiszamolnom a nyomas derivaltjat
	
	dpressure_left=(pressure[index+1]-pressure[index-1])/(2.0*dr);
	dpressure_right=(pressure[index+2]-pressure[index])/(2.0*dr);
	DPressure=(1-small)*dpressure_left+small*dpressure_right;
	*dPdr = DPressure;
	
}



void Compute_rp(double time, double step, double position, double p_radius, double *position_new, double *sigma)
{
 int i;
 double dr1, dr2, dr3, dr4, rtemp;
 double Sigma, Pressure, dPdr;

 
 Compute_Sigma_P_dP(time,0.0,position,sigma,&Sigma,&Pressure,&dPdr);	
	
 dr1=Derivs_rp(position,p_radius,Sigma,Pressure,dPdr);

 rtemp=position+0.5*step*dr1;
 Compute_Sigma_P_dP(time,0.5*step,rtemp,sigma,&Sigma,&Pressure,&dPdr);	
 dr2=Derivs_rp(rtemp,p_radius,Sigma,Pressure,dPdr);	

 rtemp=position+0.5*step*dr2;
 Compute_Sigma_P_dP(time+0.5*step,0.0,rtemp,sigma,&Sigma,&Pressure,&dPdr);	
 dr3=Derivs_rp(rtemp,p_radius,Sigma,Pressure,dPdr);
	
 rtemp=position+step*dr3;
 Compute_Sigma_P_dP(time+0.5*step,0.5*step,rtemp,sigma,&Sigma,&Pressure,&dPdr);	
 dr4=Derivs_rp(rtemp,p_radius,Sigma,Pressure,dPdr);
 

 *position_new=position+step*(dr1+2.0*dr2+2.0*dr3+dr4)/6.0;
}


double Radius(double time, double timestep, double r, double radius0, double *sigma) {
	
	double small, sigma_gas, sigma_dust, omega_kep, v_kep, csound, c_s2, Sigma_cgs, epsilon,
	       s_frag, s_df, s_min, new_radius, tau_gr, pressure, dPdr, dlnPdlnr, u_frag, u_frag2, rt;
	
	u_frag = U_FRAG;
    u_frag = u_frag/AUPDAY2CMPSEC;
	u_frag2 = u_frag*u_frag;
	
	omega_kep = G*sqrt(Central_Mass/(r*r*r));
	v_kep = r*omega_kep;
	csound = Asp_Ratio*v_kep;
	c_s2 = csound*csound;
	
	Compute_Sigma_P_dP(time,0.0,r,sigma,&sigma_gas,&pressure,&dPdr);
	
	Sigma_cgs = sigma_gas/SDCONV;
	
	dlnPdlnr = r/pressure * dPdr;
	
	s_frag = F_FRAG * 2.0/(3.0*M_PI) * Sigma_cgs/(RHO_P*alpha_turb(r)) * u_frag2/c_s2;
	
	s_df = u_frag*v_kep/fabs(dlnPdlnr*c_s2*0.5) * 2.0*Sigma_cgs/(M_PI*RHO_P);
	
	s_min = fmin(s_frag,s_df);
	
	//epsilon = sigma_gas/sigma_dust;
	epsilon = 100.0;
	
	tau_gr = epsilon/omega_kep; //time unit conversion !!!
	
	if (radius0 < s_min) rt = fmin(radius0*exp(timestep/tau_gr),s_min);
	
	if (radius0 >= s_min) rt = s_min;
	
	
	return rt;
	
}



int main(int argn,char *arg[])
{
 int i;
 double r, rp, radius, radius_old, rp_new, rdensmax;
 double sigma[imax+2];
 double t,t_old,dt,tmax;
 double dt_max,A_max;
 FILE *f_rp;
 double Proba_Sigma,Proba_Sd_exp;
 
 f_pos = fopen("pozici.dat","w");

 if(argn==1) printf("***   file_rp(time:rp)   file_sigma(r:sigma)   tmax   dt\n");
 else
 {

 // Input arguments
 f_rp=fopen(arg[1],"w");	// output file for r_planet
 f_sigma=fopen(arg[2],"w");	// output file for sigma
 tmax=atof(arg[3]);		    // time length of integration
 dt=atof(arg[4]);		    // time step

 // Time step check for the diffusion term
 A_max = -10000.0;
 for(i=1;i<=imax;i++)
 {
  r=rmin+(i-1)*dr;
  if (Func_A(r) > A_max) A_max = Func_A(r);
 }
   
 dt_max = dr*dr/(2.0*A_max);
 if(dt>dt_max)
 {
  printf("**  spatial resolution = %lg  **  corresponding maximum stepsize = %lg\n",dr,dt_max);
  exit(1);
 }


 // Initial values
 t=0.0;
 rp=Initial_rp();					                // csak egy szam t=0-ban :)
 rp_new=0.0;
 for(i=0;i<=imax+1;i++) sigma[i]=Initial_sigma(i);	// igazi racspontokra
 fprintf(f_rp,"%lg %.15E\n",t,rp);

 radius_old = 0.0001;

 // Time loop
 while(t<tmax)
 {
		 	 
  radius = Radius(t,dt,rp,radius_old,sigma);
	 
  Compute_rp(t,dt,rp,radius,&rp_new,sigma);
  rp=rp_new;
  t=t+dt;
  radius_old = radius;	 
  fprintf(f_rp,"%lg %.15E %.15E\n",t/365.25,rp,radius);
  //fprintf(f_pos,"%lg %lg\n",rp,0.00001);
  fflush(f_rp);
 }
 for(i=1;i<=imax;i++) {
     fprintf(f_sigma,"%lg %lg\n",rmin+(i-1)*dr,sigma[i]);
 }
 
 fclose(f_rp); fclose(f_sigma);
 }
}
