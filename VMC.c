#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "utils.h"

#define N 32			
#define rho 0.365		/* rho = N/V	*/
#define L 4.44         	/* V = L^3,  L = (N/rho)^1/3  */
#define Epsilon 10.22
#define D 0.92814        /* D (K sigma^2) = h_bar^2/(2m)  */
#define STEPS 32000
#define eps 1.e-5 
#define PI 3.1415926535
// delta 0.386
/**
 * UNITS:
 * Helium-4 mass = 4 a.u.
 * sigma = k_B = 1
 */
 
typedef struct  
{
	vector r[N];
} walker;


walker init(walker R)
{
	int i, j;
		for (j = 0; j < N; j++) 
			R.r[j] = randomPosition(0, L);
	return R;
}


void drawConfig(walker R)
{
	int i;
	FILE *out;
	out = fopen("config.dat", "w+");
	for (i = 0; i < N; i++)
		fprintf(out, "%lg\t%lg\t%lg\n", R.r[i].x, R.r[i].y, R.r[i].z);
}


double pbc(double c)
{
	double x = c;
	if (x < 0) {
		while (x < 0)
			x += L;
		return x;
	} else if (x > L) {
		while (x > L)
			x -= L;
		return x;
	} else return x;
}

walker newWalker(walker R, double delta)
{
	int i, d;
	walker w;
	d = rnd()*N;
	for (i = 0; i < N; i++) {
	if (i == d) {	
		w.r[i].x = pbc(R.r[i].x + (2.*rnd() - 1.)*delta);
		w.r[i].y = pbc(R.r[i].y + (2.*rnd() - 1.)*delta);
		w.r[i].z = pbc(R.r[i].z + (2.*rnd() - 1.)*delta);
	} else {
		w.r[i].x = R.r[i].x;
		w.r[i].y = R.r[i].y;
		w.r[i].z = R.r[i].z;
	}
	}
	return w;
}

double periodicDistance(vector a, vector b)
{
	vector d;
	d = vectorSum(a, opposite(b));
	if (d.x < -L/2.)
		d.x += L;
	else if (d.x > L/2.)
		d.x -= L;
	if (d.y < -L/2.)
		d.y += L;
	else if (d.y > L/2.)
		d.y -= L;	
	if (d.z < -L/2.)
		d.z += L;
	else if (d.z > L/2.)
		d.z -= L;
	return sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}

double V(double r)
{	
	if (r == 0)
		return 4.*Epsilon*(pow(eps, -12) - pow(eps, -6));
	else 
		return 4.*Epsilon*(pow(r, -12) - pow(r, -6));
}

double u(double b, walker R)
{
	int i, j;
	double r, sum;
	sum = 0;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			r = periodicDistance(R.r[i], R.r[j]);
			if (r == 0)
				r = eps;
			sum += pow(b/r, 5);
		}
	}
	return sum;
}

double quantities(double b, walker R)
{
	int i, j;
	double potential, kinetic, E, r, r_m;
	r_m = 0;
	potential = 0;
	kinetic = 0;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			r = periodicDistance(R.r[i], R.r[j]);
			if (r == 0)
				r = eps;
			potential += V(r);
			kinetic += 10.*D*pow(b/r, 5)/r/r;
			r_m += r;
		}
	}
	E = kinetic + potential;
	r_m /= 0.5*N*(N - 1);
	return E;
}



void VMC(walker R,double delta)
{
	int i, accepted, ns;
	double b, ratio, energy, d, q, potential_correction, kinetic_correction;
	walker temp;
	FILE *out;
	out = fopen("energy.dat", "w+");
	d = delta;
	for (b = 1.; b < 1.5; b += 0.01) {
	accepted = 0;
	energy = 0;
	ns = 9*STEPS/10;
	R = init(R);
	for (i = 0; i < STEPS; i++) {

		temp = newWalker(R, d);
		ratio = exp(-0.5*(u(b, temp) - u(b, R)));
		if (ratio*ratio > rnd()) {
			R = temp;
			accepted++;
		}
		q = quantities(b, R);
		if (i > STEPS/10.)
			energy += q;
		//fprintf(out, "%lg\n", q/N);
	}
	kinetic_correction = 10.*PI*D*rho*pow(b, 5)*pow(2./L, 4);
	potential_correction = 8./3.*PI*Epsilon*rho*(pow(2./L, 9)/3. - pow(2./L, 3));
	energy = energy/ns/N + potential_correction + kinetic_correction;
	fprintf(out, "%lg\t%lg\n", b, energy);
	printf("Acceptance ratio: %lg%%\n", (double) accepted/STEPS*100);
	printf("%lg\t%lg\n", b, energy);
	}
	fclose(out);
}

/******************************************************************************************/

int main()
{	
	double delta;
	walker R;
	printf("Delta = ");
	scanf("%lg", &delta);
	srand(29);
	VMC(R, delta);
	return 0;
}

























