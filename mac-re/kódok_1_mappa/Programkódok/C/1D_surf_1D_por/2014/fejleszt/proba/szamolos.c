#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void gen(double x[][2], double t, int y) {

	int i;
	for(i = 0; i < 40; i++) {
		x[i][0] = y+i;
		if(t==0) x[i][1] = i;
	}	
}


void ize(double x[][2], int ii, int io, int in, int on){

	double index;
	int i;
	
	for(i = 0; i < 40; i++) {

		index = x[i][0];
		if(ii <= (int)index && (int)index <= io) {
			printf("ii: %i, index: %lg  io: %i\n",ii,index,io);	
		}

	}

	getchar();
}


int main(){

	double x[40][2], t,z;
	int y, ii, io, in, on;

	ii = 5;
	io = 25;
	
	t = 0;
	y = 0;
	z = 0;
	do {
		y = floor(t+z);
		gen(x,t,y);
		z = z+.2;

		ize(x,ii,io,in,on);
		
		t=t+.5;

	} while (t < 500);
	printf("mikell?\n");
	return 0;

}
