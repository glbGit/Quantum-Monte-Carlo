#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct 
{
	double x;
	double y;
	double z;
} vector, point;

vector vectorSum(vector a, vector b)
{
	vector v;
	v.x = a.x + b.x;
	v.y = a.y + b.y;
	v.z = a.z + b.z;
	return v;
}

double scalarProduct(vector a, vector b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
} 

vector vectorProduct(vector a, vector b)
{
	vector v;
	v.x  = a.y*b.z - a.z*b.y;
	v.y  = a.z*b.x - a.x*b.z;
	v.z  = a.x*b.y - a.y*b.x;
	return v;
} 

vector opposite(vector v)
{
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	return v;
}

double distance(vector a, vector b)
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}

vector randomPosition(double a, double b)
{
	vector v;
	v.x = (double) rand()/RAND_MAX*(b - a) + a;
	v.y = (double) rand()/RAND_MAX*(b - a) + a;
	v.z = (double) rand()/RAND_MAX*(b - a) + a;
	return v;
}

void printVector(vector v)
{
	printf("(%lg, %lg, %lg)\n", v.x, v.y, v.z);
}









