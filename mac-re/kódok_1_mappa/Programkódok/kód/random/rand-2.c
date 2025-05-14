//
//generate random numbers with spectra like x^(-2)
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void main(){
    long N = 30;        //number of generated randoms
    float r=0;              //uniform distributed random in [0:1]
    float alpha=-2;         //power of spectra
    float x0=1;             //x min of generated randoms
    float x1=10;            //x max of generated randoms

    float p0=0, p1=0;       //temporary var for powers
    float E=0;              //generated random
    char* outname = "r.txt";  //name of file of generated randoms

    FILE* OUT;

    OUT = fopen(outname, "wt");

    p0 = pow(x0, alpha + 1);
    p1 = pow(x1, alpha + 1);


    for(long i=0; i<N; i++){
	r = (float)rand()/RAND_MAX;
	if(r==0){                     //awoid to take log(0)
	    continue;
	}
	E = pow(r * (p1 - p0) + p0, 1.0/(alpha + 1.0));
	fprintf(OUT, "%f\n", E);
    }
    fclose(OUT);
}
