#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main() {

	double x[10] = {1,2,5,5,5,0,2,2,9,10}, y[10],p, ptemp;
	double lv[10] = {1,5,1,5,2,2,6,7,8,9}, sig, st[10];
	int i, j=0, k = 0, l = 0;
	
	for(i = 0; i < 10; i++) {
		y[i] = 0;
		st[i] = 0;

	}

	i = 0;
	do {
		if(x[i] != x[i+1]) {
			y[j] = x[i];
			st[j] = lv[i];
			sig = 0;
			k = 0;
			j++;
			i++;
		} else {
			do {
				y[j] = x[i];
				sig = sig + lv[i+k];
				st[j] = sig;
				k++;
			} while (x[i] == x[i+k]);
			printf("i: %d  j: %d\n",i, j);
			
			i = i+k;
			k =0;
			j++;


		}


	} while (i < 10);


	for(i = 0; i < 10; i++) printf("y: %lg   st: %lg\n",y[i],st[i]);
}
