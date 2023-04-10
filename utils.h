#include <stdlib.h>

double rnd() 
{
	return (double) rand()/RAND_MAX;
}

void zeroArray(double *x, int N) 
{
	int i;
	for (i = 0; i < N; i++)
		x[i] = 0;
}
