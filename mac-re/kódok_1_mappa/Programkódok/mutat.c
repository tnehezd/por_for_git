#include<stdio.h>
#include<stdlib.h>

int main()
{
	int *p;	
	int x=7, y=3;
	p=&x;
	*p=*p+x+y;
	printf("%d\n",x);

	return 0;
}
