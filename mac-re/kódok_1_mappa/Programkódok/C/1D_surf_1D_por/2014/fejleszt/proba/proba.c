#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main() {

	double x[10] = {1, 2, 5, 5, 5, 0, 2, 2, 9, 10}, y[10],p, ptemp, lv[10] = {1,5,1,5,2,2,6,7,8,9}, sig, st[10];
	int i, j=0, k = 0, inti, into, intiv[10], intov[10], l = 0;
	
	for(i = 0; i < 10; i++) {
		y[i] = 0;
		st[i] = 0;
		intiv[i] = 0;
		intov[i] = 0;
	}

	for(i = 0; i < 10; i++) {
		sig = 0;
		if(x[i] == x[i+1]) {
			p = x[i];
			if(k == 0) {
				ptemp = p;
			}

			if(p != ptemp) {
				k = 0;
			}
			k++;
			if(k == 1) {
				inti = i;
				l++;
				j++;
			}
			i++;
			into = i;
			ptemp = p;
			intiv[l-1] = inti;
			intov[l-1] = into;
		} else {
			p = x[i];
			j++;
		}
		y[j-1] = p;
	}

	int ki = 0, ko = 0;
	j = k = 0;
	for(i = 0; i < 10; i++) {
		if(intiv[i] != 0 && intov[i] != 0) {
			ki = intiv[i];
			ko = intov[i];
			sig = 0;
			printf("ki %d   ko %d   s %lg   s %lg\n",ki,ko, lv[ki],lv[ko]);
				for(j = ki; j <= ko; j++) {
					printf("\n\n j : %d,  lv: %lg\n",j,lv[j]);
					sig = sig + lv[j];
				}
			printf("s: %lg\n",sig);
			st[k] = sig;
			k++;
		} else {
			k++;
		}

	}

	for(i = 0; i < 10; i++) printf("y: %lg   st: %lg  ii: %d  io: %d\n",y[i],st[i], intiv[i], intov[i]);
}
