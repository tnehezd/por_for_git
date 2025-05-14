#include<stdio.h>
#include<stdlib.h>
#include<string.h>


#define filenev1 "temp2.dat"

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
void data_be(double data[][3], int db) {

	double dat1, dat2, dat3;
	int i;

	fin1 = fopen(filenev1,"r");

	for (i = 0; i < db; i++) {			
           	if(fscanf(fin1,"%lg  %lg %lg",&dat1,&dat2,&dat3) == 3) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 14, mert 14 oszlopot szeretnénk beolvasni.	*/
			data[i][0] = dat1;
			data[i][1] = dat2;
			data[i][2] = dat3;

	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
//			getchar();
			exit(EXIT_FAILURE);
   	        }

	}
	
	fclose(fin1);
	
}



int main(int argc, const char **argv) {

	int db = 0, out = 0, i;
	
	db = sorok_szama(out);

	double data[db][3];
	data_be(data,db);

	char *p = malloc(sizeof(argv[1]));
	strcpy(p, argv[1]);

	printf("a: %i  a: %s\n",argc,p);
	

	if(argc == 2) {
		fout=fopen(p,"w");
	
		for(i = 0; i < db; i++) {
			fprintf(fout,"%d   %lg  %lg\n",(int)data[i][0],data[i][1],data[i][2]);
				if(data[i][0] != data[i+1][0] && i+1 < db) fprintf(fout,"\n\n");
		}
	
		fclose(fout);

	} else {
		printf("******************ERROR*********************\n\nPlease type in output file name!\n");
	}
	
	return 0;

}
