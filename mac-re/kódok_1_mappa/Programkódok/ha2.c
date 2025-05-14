#include<stdio.h>
#include<math.h>

#define rho 1.501780435e-3
#define sigma0 0.0001			
#define ngrid 100
#define rmin 0.1
#define rmax 15.0
#define dgrid (rmax-rmin)/(ngrid-1)	//dr

#define sdexp -0.5 			//surface density profile exponent 
#define alpha_visc 0.01
#define asp_ratio 0.05
#define G 1.0 				//ekkor r=1-nél a periódus 2*pi 
#define sb_const 2.533732811E14        // W/AU/AU/K^4
#define T_bg 10.0			//background temperature 10K
#define M_star	1.0

FILE *fil;


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

int main(){

 double fx, temp, tempreg, a,b,c, k0, h,k,l, eps, r, t_midplane[ngrid+2], t, dt, tmax, kappa, optd, opteff;
 int i,j;
 
  t=0.0;
  dt=0.1;
  tmax=10.0;
 
  
  fil=fopen("temperature.dat","w");


/*  for(i=1;i<=ngrid;i++) {
  
   r=rmin+(i-1)*dgrid;
   
 }
*/
   l=4.6E3*pow(rho,1/15);
   h=1.1E4*pow(rho,1/21);
   k=3.0E4*pow(rho,4/75);
 
   eps=1e-4;

   temp=1000.0;  
  
   for(i=1;i<ngrid;i++){
     
   r=rmin+(i-1)*dgrid;     
     
     do{
     
       tempreg=temp;
       
       if(temp<=170.){ 
	 k0=2.250419625E3;
	 a=0.0;
	 b=2.0;	 	 
      }
      
      
      
      else if(temp>170.0, temp<=210.0){
	k0=2.250419637E23;
	a=0.0; 
	b=-7.0;	 
      }
      
      else if(temp>210.0, temp<=l){
	k0=5.626049092E29;
	a=0.0; 
	b=1.0;	
      }
      
      else if(temp>l, temp<=3000.0){
	k0=2.250419637E41;
	a=2.0/3.0; 
	b=-9.0;	
      }
      
      else if(temp>3000.0, temp<=h){
	k0=2.250419637e-1;
	a=2.0/3.0; 
	b=3.0;	
      }
      
      else if(temp>h, temp<=k){
	k0=1.125209818e-29;
	a=1.0/3.0; 
	b=10.0;	
      } 
      
      else if(temp>k){
	k0=1.687814744E28;
	a=1.0; 
	b=-5.0/2.0;	
      }
    
/*    kappa=k0*pow(rho,a)*pow(tempreg,b);
    optd=kappa*sigma0/2.0;
    opteff=(3.0*optd/8.0)+(sqrt(3.0)/4.0)+(1/(4.0*optd));
    
    temp=pow((opteff*(9.0/4.0*sigma0*viscosity(r)*kep_freq(r)*kep_freq(r))+Coeff_m()*T_bg*T_bg*T_bg*T_bg)/Coeff_m(),1.0/4.0);
*/    
    temp=pow((3.0*k0*pow(rho,a)*pow(tempreg,b)*sigma0/2.0/8.0+sqrt(3)/4.0+1.0/(4.0*k0*pow(rho,a)*pow(tempreg,b)*sigma0/2.0)*(9.0*sigma0*viscosity(r)*kep_freq(r)*kep_freq(r)/4.0)+2.0*sigma0*T_bg*T_bg*T_bg*T_bg)/(2.0*sigma0),1.0/4.0);	
	
//    printf("tempreg: %.8lf, temp: %.8lf, hiba: %.8lf\n",tempreg, temp, fabs(tempreg-temp)); 
//    getchar();
//     temp=temp_temp;
     }while(fabs(tempreg-temp)>1.0e-6);   
    
    t_midplane[i]=temp;
    
     printf("r: %.07lf t_midplane: %0.7lf\n",r, t_midplane[i]);


   } 
           
  
  do{
    
    for(i=0;i<ngrid;i++){
  
      t_midplane[i]=t_midplane[i]+dt*(Coeff_1(r)*(t_midplane[i+1]-2*t_midplane[i]+t_midplane[i-1])/(dgrid*dgrid)+Coeff_2(r)*(t_midplane[i+1]-t_midplane[i-1])/(2.0*dgrid));
      
//      printf("%.20lf\n", t_midplane[i]);
    } 
   
   t=t+dt;
  
//   printf("time: %lf\n", t);
      
  }while(t<tmax);


  for(i=0;i<ngrid;i++)      fprintf(fil, "%0.7lf\n",t_midplane[i]);


/*
 ffx(temp);
 ggx(temp,k0,a,b,r);
*/

} 	
