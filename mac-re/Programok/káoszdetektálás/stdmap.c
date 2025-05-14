
/* Standard Map -- 2012.06.14. */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>


/************************************************************************************************/


double my_rand()
 {
  return (double) rand()/RAND_MAX;
 }


int main(int argn,char *arg[])					/*  a parancssorbol keri be az adatokat  */
{
 double K;
 long double x,p,zero;
 int n,nmax,init,initmax;
 FILE *fmap,*fmap0;

 if(argn==1) printf("***   K   nmax   initmax   file \n");	/*  K: nemlin param  */
 else
 {

  K=atof(arg[1]);
  nmax=atoi(arg[2]);						/*  iteracio szama  */
  initmax=atoi(arg[3]);						/*  hany kezdeti feltetel van, ebbol random modon valogat  */
  fmap=fopen(arg[4],"w");
  fmap0=fopen("stdmap0.d","w");

  zero = 0.0;

  if(K<0)
  {
   printf(" ===>   K should be positive.\n");
   exit;
  }

  for(init=1;init<=initmax;init++)				/*  kezdeti feltetelek generalasa  */
  {

/*   x=fmod(rand(),2*M_PI);					  x-et random vesszuk mod 2*pi */
/*   p=fmod(rand(),2*M_PI);							         */

   x=2*M_PI*my_rand();
   p=2*M_PI*my_rand();

   if(x<zero) x=x+2*M_PI;
   if(p<zero) p=p+2*M_PI;
   if(x>2*M_PI) x=x-2*M_PI;
   if(p>2*M_PI) p=p-2*M_PI;

   fprintf(fmap0,"%Lg %Lg\n",x/(2*M_PI),p/(2*M_PI)); 		/*  kezdeti feltetelek vannak benne  */

   for(n=1;n<=nmax;n++)						/*  hanyszor iteraljuk az egyes palyakat  */
   {
    p=fmod(p-K*sin(x+M_PI),2*M_PI);
    x=fmod(x+p+M_PI,2*M_PI);

    if(x<zero) x=x+2*M_PI;
    if(p<zero) p=p+2*M_PI;
    if(x>2*M_PI) x=x-2*M_PI;
    if(p>2*M_PI) p=p-2*M_PI;

  
    fprintf(fmap,"%Lg %Lg\n",x/(2*M_PI),p/(2*M_PI)); 
    fflush;
   }
  }

 } 

 fclose(fmap0); fclose(fmap);

}

/************************************************************************************************/
