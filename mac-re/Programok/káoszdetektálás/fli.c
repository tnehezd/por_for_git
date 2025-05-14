/*  	KISZAMOLJA EGY ALTALUNK MEGADOTT PALYA FLI-JET  	*/
/* 	Standard Map -- 2013. 01. 31. 		*/

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>


/************************************************************************************************/


double my_rand() {
  	return (double) rand()/RAND_MAX;

}


int main(int argn,char *arg[])	{				/*  a parancssorbol keri be az adatokat  */

 	double K, log_wlength;
 	long double x,p,zero,w_1,w_2,w_old_1,w_old_2,wlength,lce, small;
 	int n,nmax,init,initmax;
 	FILE *fmap,*fmap0;

 	if(argn==1) printf("***   K   nmax   0<x<1   0<p<1   orbitfile  vecfile \n");	/*  K: nemlin param  */
 	else {
 	 	K=atof(arg[1]);
  		nmax=atoi(arg[2]);						/*  iteracio szama; atoi --> ascii-bol integer  */
  		x=atof(arg[3]);						/* atof --> asciibol float  */
  		p=atof(arg[4]);
 
  		fmap=fopen(arg[5],"w");
  		fmap0=fopen(arg[6],"w");
  
	 	w_old_1=1.0;
  		w_old_2=0.0;  

  		zero = 0.0;

  		if(K<0) {
   			printf(" ===>   K should be positive.\n");
   			exit;
  		}

   
	   	x=x*2*M_PI;
   		p=p*2*M_PI;
   
	   	if(x<zero) x=x+2*M_PI;
   		if(p<zero) p=p+2*M_PI;
   		if(x>2*M_PI) x=x-2*M_PI;
   		if(p>2*M_PI) p=p-2*M_PI;

   		lce=0.0;
		small=-HUGE_VAL;

   		for(n=1;n<=nmax;n++) {						/*  hanyszor iteraljuk az egyes palyakat  */
   
    			p=fmod(p-K*sin(x+M_PI),2*M_PI);
    			x=fmod(x+p+M_PI,2*M_PI);

    			w_1=w_old_1+K*cos(x)*w_old_2; 
    			w_2=w_old_1+w_old_2+K*cos(x)*w_old_2;

    			wlength=sqrt(w_1*w_1+w_2*w_2);
			   
    			lce=lce+log(wlength);		

			log_wlength=log(wlength);

			if(log_wlength>small){
				small=log_wlength;
			}  
   
    			printf("lepes: %i, elteres vektor hossza: %Lg,   LCE: %Lg,  small: %Lg\n",n,wlength,lce/n,small);
    			fprintf(fmap0,"%i %Lg  %Lg %Lg\n",n,wlength,lce/n, small);
    			fflush(fmap0);
/*    getchar();	*/
 
    			w_old_1=w_1; /* /wlength;*/
    			w_old_2=w_2; /* /wlength; */
   
    			if(x<zero) x=x+2*M_PI;
    			if(p<zero) p=p+2*M_PI;
    			if(x>2*M_PI) x=x-2*M_PI;
    			if(p>2*M_PI) p=p-2*M_PI;

  
    			fprintf(fmap,"%Lg %Lg\n",x/(2*M_PI),p/(2*M_PI)); 
    			fflush;
   		}

 	} 

 	fclose(fmap), fclose(fmap0);

}

/************************************************************************************************/
