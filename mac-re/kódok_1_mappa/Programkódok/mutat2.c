#include<stdio.h>
#include<stdlib.h>

int main()
{
	int *p;	
	int x=7, y=3;
	int **q;

/*	Az int típus tárolására alkalmas terület lefoglalása a művelet eredményességének ellenőrzzésével	*/
//	p=(int *)malloc(sizeof(int));
//	if (p==NULL) exit(-1);

/*	Az összeg tárolása és kiírása	*/
//	*p=x+y;


//	printf("Az összeg: %d\n",*p);


/*	A lefoglalt terület felszabadítása	*/
//	free(p);


	p=&x;
	q=&p;
	x=x+ *p + **q;

	printf("x=%d\n",x);

}
