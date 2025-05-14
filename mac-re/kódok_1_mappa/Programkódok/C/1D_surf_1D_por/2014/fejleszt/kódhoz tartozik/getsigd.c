void getSigDustInterpol(double sigin[][3], double sigmicin[][3], double *sigout, double *sigmicout) {

	int i, j;
	double intsig[PARTICLE_NUMBER], intsigmicr[PARTICLE_NUMBER];

	for(i=0;i<PARTICLE_NUMBER;i++) {

		double a, b;

		for(j=0; j<PARTICLE_NUMBER;j++) {

			double temp = 0;
			double simt = 0;

			if(sigin[i][1] <= sigmicin[j][1] && sigin[i+1][1] > sigmicin[j][1]){

/*	Megnezi, hogy a mikronos porszemcse hol helyezkedik el a cm-eshez kepest: ha a mikoronos reszecske tavolsaga 2 cm-es kozott van: 	*/
				a = (sigin[i+1][0] - sigin[i][0])/(sigin[i+1][1] - sigin[i][1]);	/*	akkor a mikronos tavolsagara egyenes illesztessel interpolalja a cm-es porszemcsekbol szarmazo feluletisuruseget	*/
				b = sigin[i][0] - a * sigin[i][1];
				simt = a * sigmicin[j][1] + b;
				temp = temp + simt;
			} 

			if(temp != 0) intsigmicr[j] = temp;	/*	ide rakja az interpolalt feluletisuruseget	*/		
			
			double temp2=0;
			double simt2 = 0;

/*	Az elozohoz hasonloan interpolalja a mikronos feluletisuruseget a cm-es porszemcsek tavolsagara		*/	
			if((sigmicin[i][1] <= sigin[j][1]  && sigmicin[i+1][1] > sigin[j][1])) {// && sigmicin[i][0] > RMIN) {
				a = (sigmicin[i+1][0] - sigmicin[i][0])/(sigmicin[i+1][1] - sigmicin[i][1]);
				b = sigmicin[i][1] - a * sigmicin[i][1];
				simt2 = a * sigin[j][0] + b;
//				printf("ri: %lg r: %lg  ro: %lg\n", sigdtemp[i][0], simt, sigdtemp[i+1][0]);
				temp2 = temp2 + simt2;
			}

			if(temp2 != 0) intsig[j] = temp2;	/*	elmenti az interpolalt feluletisuruseget	*/		
		}

	}

	for(i=0; i < PARTICLE_NUMBER; i++) {
/*	itt adja ossze az interpolalt feluletisuruseget az adott meretbol szarmazo feluletisuruseggel	*/	
		sigout[i] = sigin[i][0] + intsig[i];
		sigmicout[i] = sigmicin[i][0] + intsigmicr[i];
	}
}

