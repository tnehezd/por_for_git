#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define ngrid 1000
#define rmin 0.1
#define rmax 15.0
#define dgrid (rmax-rmin)/(ngrid-1)	//dr

#define sigma0 0.0001			//(M_Nap/AU^2) ---> midplane ro kiszámolható (jegyzetből)
#define sdexp -0.5 			//surface density profile exponent 
#define alpha_visc 0.01
#define asp_ratio 0.05

#define G 1.0 				/* ekkor r=1-nél a periódus 2*pi */
#define rho 0.01683361792   /* ---> rho0 kiszámolva, innen majd rho is lesz */
#define sb_const 1.26903855E19        /* W/AU/AU/K^4 */
#define T_bg 10.0			/* background temperature 10K */
#define M_star	1.0

#define mean_molweight (1.156406492*pow(10.0,-33.0))
#define k_boltzmann (1.3806504*pow(10.0,-23.0))
#define m_prot (8.40969858*pow(10.0,-58.0))




FILE *fmo, *fil;
 
double viscosity(double r){
 
  double nu, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, alpha_r;

  r_dze_i=0.0;
  r_dze_o=10.0;
  Dr_dze_i=0.0;
  Dr_dze_o=0.5;
  a_mod=0.01;
  
  alpha_r=1.0-0.5*(1.0-a_mod)*(tanh((r-r_dze_i)/Dr_dze_i)+tanh((r_dze_o-r)/Dr_dze_o));
  
  nu=alpha_visc*alpha_r*asp_ratio*asp_ratio*G*sqrt(r);
  
  return nu;
  
}


double Coeff_1(double r){
  
  double A;
  
  A=3.0*viscosity(r);
  
  return A;
  
}

double Coeff_2(double r){
  
  double B;
  
  B=9.0*viscosity(r)/(2*r);
  
  return B;
    
}


void Perem(double *u){
  
  u[0]=(u[1]/viscosity(rmin))*viscosity(rmin-dgrid);
  u[ngrid+1]=(u[ngrid]/viscosity(rmin+(ngrid)*dgrid))*viscosity(rmin+(ngrid+1)*dgrid);
  
}


void Initial_Profile(double *u){
 // u a nu*sigma  
 //*u-val mondom meg, hogy tömb, a for-ban már u[i]-val töltöm fel

  int i;
  double r;
  
  Perem(u);

  for(i=1;i<=ngrid;i++)
  {
    
    r=rmin+(i-1)*dgrid;					//mindig az adott r-el tölti fel az u tömböt
    u[i]=sigma0*viscosity(r)*pow(r,sdexp);
  }
  
}

double Coeff_m(){
 
  double m_sigma;
  
  m_sigma=2*sb_const;
  
  return m_sigma;
  
}


double kep_freq(double r){
  
  double omega;
   
  omega=sqrt((G*M_star)/pow(r,3));
  
  return omega;
  
}




int main(){

  int i;
  double dt, tmax, t, r, u[ngrid+2];
  double temp, tempreg, a,b,c, k0, h,k,l, eps, t_midplane[ngrid+2], kappa, optd, opteff;
  
  t=0.0;
  dt=0.1;
  tmax=100.0;
  
  fil=fopen("temperature_2.dat","w");
  fmo=fopen("surface_2.dat","w");

  Initial_Profile(u);

 
   eps=1.0e-6;

   temp=50.0;  

  
  
  do {
  
  for(i=1;i<=ngrid;i++) {
    
    r=rmin+(i-1)*dgrid;
    u[i]=u[i]+dt*(Coeff_1(r)*(u[i+1]-2*u[i]+u[i-1])/(dgrid*dgrid)+Coeff_2(r)*(u[i+1]-u[i-1])/(2.0*dgrid));
    
  }

  Perem(u);
  t=t+dt;

  printf("számolok, t = %.4lf\n",t);  

  } while(t<tmax);

    
 for(i=1;i<=ngrid;i++) {
  
   //tmax-ban veszi az u(i)-t, és azt iratjuk ki
   
   r=rmin+(i-1)*dgrid;
   fprintf(fmo, "%lg   %lg\n", r, u[i]/viscosity(r)); 		//csak a szigma kell
   
 }
  
 

   for(i=1;i<ngrid;i++){
     
   r=rmin+(i-1)*dgrid;     
 
   l=4.6E3*pow(rho,1.0/15.0);
   h=1.1E4*pow(rho,1/21);
   k=3.0E4*pow(rho,4/75);
     
     do{
     
       tempreg=temp;
       
       if(temp<=170.){ 
	 k0=2.250419625E3;
	 a=0.0;
	 b=2.0;	 	 
      }
      
      
      
      else if(temp<=210.0){
	k0=2.250419637E23;
	a=0.0; 
	b=-7.0;	 
      }
      
      else if(temp<=l){
	k0=5.626049092E29;
	a=0.0; 
	b=1.0;	
      }
      
      else if(temp<=3000.0){
	k0=2.250419637E41;
	a=2.0/3.0; 
	b=-9.0;	
      }
      
      else if(temp<=h){
	k0=2.250419637e-1;
	a=2.0/3.0; 
	b=3.0;	
      }
      
      else if(temp<=k){
	k0=1.125209818e-29;
	a=1.0/3.0; 
	b=10.0;	
      } 
      
      else if(temp>k){
	k0=1.687814744E28;
	a=1.0; 
	b=-5.0/2.0;	
      }
    
    temp=pow((3.0*k0*pow(rho,a)*pow(tempreg,b)*sigma0/2.0/8.0+sqrt(3)/4.0+1.0/(4.0*k0*pow(rho,a)*pow(tempreg,b)*sigma0/2.0)*(9.0*sigma0*viscosity(r)*kep_freq(r)*kep_freq(r)/4.0)+2.0*sigma0*T_bg*T_bg*T_bg*T_bg)/(2.0*sigma0),1.0/4.0);	
	
     }while(fabs(tempreg-temp)>eps);   
    
     t_midplane[i]=temp;
    
     printf("r: %.07lf t_midplane: %0.7lf\n",r, t_midplane[i]);


   } 
           
  
  do{
    
    for(i=0;i<ngrid;i++){
  
      t_midplane[i]=t_midplane[i]+dt*(Coeff_1(r)*(t_midplane[i+1]-2*t_midplane[i]+t_midplane[i-1])/(dgrid*dgrid)+Coeff_2(r)*(t_midplane[i+1]-t_midplane[i-1])/(2.0*dgrid));
      
    } 
   
   t=t+dt;
      
  }while(t<tmax);


  for(i=0;i<ngrid;i++){      
    r=rmin+(i-1)*dgrid;     
    fprintf(fil, "%.4lg   %.8lg\n", r,t_midplane[i]);
  }
 
    
}
