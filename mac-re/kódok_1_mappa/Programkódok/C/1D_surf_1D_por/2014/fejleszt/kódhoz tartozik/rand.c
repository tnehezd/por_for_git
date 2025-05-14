#include<stdlib.h>
#include<stdio.h>
#include<math.h>


#define filenev1 "por_r.dat"

FILE *fin1, *fout;

/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
int sorok_szama(int numout){

	char c;
	fin1 = fopen(filenev1,"r+");		
	numout = 0;

/*	A porreszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal				*/

		while((c = fgetc(fin1)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 							*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert, ezzel beolvastuk, hogy hany soros a file				*/
	fclose(fin1);	

	return numout;

}



/*	Adatok beolvasasa	*/
void data_be(double *r, int db) {

	double dat, dummy;
	int i;

	fin1 = fopen(filenev1,"r");

	for (i = 0; i < db; i++) {			
           	if(fscanf(fin1,"%lg %lg %lg",&dummy,&dummy,&dat) == 3) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 14, mert 14 oszlopot szeretnénk beolvasni.	*/
			r[i] = dat;

	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
//			getchar();
			exit(EXIT_FAILURE);
   	        }

	}
	
	fclose(fin1);
	
}


/*	random general szamokat locmin és locmax kozott --> ez egyenletes eloszlast ad	*/
double random_gen(double locmin,double locmax) {

	double f;
	f = (double)rand() / RAND_MAX;
	return (locmin + f * (locmax - locmin));

}




/* uniform distribution, (0..1] */

double drand() {
  	return (rand()+1.0)/(RAND_MAX+1.0);
}


 /* normal distribution, centered on 0, std dev 1 */
double random_normal() {

  	return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

int main() {

	srand(time(NULL));				/*	a randomszam generalast cpu idohoz koti, igy valoban random szamokat kapunk		*/
	int db = 0, out = 0;
	db = sorok_szama(out);

 	double locmin, locmax, phi[db], x[db], y[db], z[db], r[db];	/*	min max ertekek a random hely- es sebessegvektor generalashoz	*/
	locmax = 0.0;	/*	csillagaszati egysegben, r-hez sqrt(x*x + y*y) koordinatakhoz	*/
	locmin = 2.0 * M_PI;

  	int i;

	data_be(r,db);

  	for (i=0; i<db; i++) {
		phi[i] = random_gen(locmin,locmax);
		x[i] = r[i] * cos(phi[i]);
		y[i] = r[i] * sin(phi[i]);
		z[i] = 0.5*r[i]*random_normal();
	}

	fout = fopen("xyz_gen.dat","w");

	for(i=0; i<db; i++) fprintf(fout,"%lg %lg %lg\n",x[i],y[i],z[i]);	

	return 0;

}
