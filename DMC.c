#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "utils.h"

#define N 108			
#define rho 0.365		
#define L 6.65       
#define Epsilon 10.22
#define b 1.2
#define D 0.92814 
#define PI 3.141592653
#define dTau 0.001
#define E_0 -7.1
#define RES 1000 
#define P 400
#define N_MAX 1000
#define eps 1.e-4


typedef struct  
{
	int isAlive;
	vector r[N];
} walker;

int OVERGROWTH;
int EXTINCTION;
double energy_correction = 10.*PI*D*rho*pow(b, 5)*pow(2./L, 4) + 8./3.*PI*Epsilon*rho*(pow(2./L, 9)/3. - pow(2./L, 3));

/**
 * Extracts a random number in 
 * range [-L,L] according to 
 * normal distribution.
 */
double rgauss(double m, double sigma2) 
{
	int FOUND = 0;
	double x;
	while (!FOUND) {
		x = (rnd() - 0.5)*L;
		if (exp(-0.5*(x - m)*(x - m)/sigma2)/sqrt(2*PI*sigma2) > rnd())
			FOUND = 1;
	}	
	return x;
}

/**
 * Histrogram of N_t values 
 * of v in range [-L,L]  
 * printed to "histo.dat".
 */
void histo(double *v)
{
	int i, j, k;
	FILE *out;
	double h[RES];
	out = fopen("histo.dat", "w+");
	zeroArray(h, RES);
	for (i = 0; i < RES; i++)
		for (j = 0; j < P; j++)
			if ((double) i*2*L/RES - L <= v[j] && v[j] < (double) (i + 1)*2*L/RES - L)
				h[i]++;
	k = 0;
	for (i = 0; i < RES; i++) {
		fprintf(out,"%g	%g\n", (double) i*2*L/RES - L, h[i]);	
		k += h[i];
	}	
	printf("hist_sum = %d\n", k);
	fclose(out);
}

/*************************************************************************************************/

void kill(walker *R)
{
	int i;
	for (i = 0; i < N_MAX; i++) 
		R[i].isAlive = 0;
}


void init(walker *R)
{
	int i, j;
	for (i = 0; i < P; i++) {
		R[i].isAlive = 1;
		for (j = 0; j < N; j++) {
			R[i].r[j] = randomPosition(0, L);
		}
	}
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
	int i, atom;
	double sigma2;
	walker w;
	sigma2 = 2.*D*dTau;
	w.isAlive = 1;
	atom = rnd()*N;
	for (i = 0; i < N; i++) {
	if (i == atom) {	
		w.r[i].x = pbc(R.r[i].x + sigma2*rgauss(0, sigma2));
		w.r[i].y = pbc(R.r[i].y + sigma2*rgauss(0, sigma2));
		w.r[i].z = pbc(R.r[i].z + sigma2*rgauss(0, sigma2));
	} else {
		w.r[i].x = R.r[i].x;
		w.r[i].y = R.r[i].y;
		w.r[i].z = R.r[i].z;
	}
	}
	return w;
}

double periodicDistance(vector r1, vector r2)
{
	vector d;
	d = vectorSum(r1, opposite(r2));
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

double E_L(walker R)
{
	int i, j;
	double r, v, t;
	t = 0;
	v = 0;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			r = periodicDistance(R.r[i], R.r[j]);
			if (r < 0.946) {
				v += 0;
				t += 0;
			} else {
				v += V(r);
				t += 10.*D*pow(b/r, 5)/r/r;
			}				
		}
	}
	return t + v;
}

double u(walker R)
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

//new set of configurations
void addToList(walker *a, walker w)
{
	int i = 0, searching = 1;
	while (i < N_MAX && searching) {
		if (a[i].isAlive == 0) {
			a[i] = w;
			searching = 0;
		} else i++;
		if (i == N_MAX)
			OVERGROWTH = 1;
	}
}

//Monte Carlo Step
int DMC_step(walker *R, double *E, int nw, int gen)
{
	int i, j, m, n = 0;
	walker newConfig;
	double W[N_MAX], E_l[N_MAX], sum1, sum2, E_T;
	E_T = (E_0 - energy_correction)*N + log((double)P/nw)/dTau;
	printf("factor = %lg\n", log((double)P/nw)/dTau);
	sum1 = 0;
	sum2 = 0;
	for (i = 0; i < nw; i++) {
		if (R[i].isAlive) {
			newConfig = newWalker(R[i], 5);					//diffusion
			E_l[i] = E_L(newConfig);								
			W[i] = exp(-dTau*(E_l[i] - E_T));					//weight
			sum1 += W[i]*E_l[i]; 
			sum2 += W[i];								
			m = (int) (W[i] + rnd());						//branching
			R[i].isAlive = 0;
			for (j = 0; j < m; j++)
				addToList(R, newConfig);
		}
		if (OVERGROWTH)
			break;
	}
	for (i = 0; i < N_MAX; i++)
		if (R[i].isAlive)
			n++;
	if (n == 0)
		EXTINCTION = 1;
	if (OVERGROWTH || EXTINCTION)
		return n;
	else {
		E[gen] = sum1/sum2;
		return n;
	}
	
}

//Diffusion Monte Carlo
double DMC(walker *R) 
{
	int i, nw = P, M = 3/dTau;
	double E[M], energy;
	FILE *energy_f;
	energy_f = fopen("energy.dat", "w+");
	energy = 0;
	for (i = 0; i < M; i++) {
		nw = DMC_step(R, E, nw, i);
		
		if (OVERGROWTH) {
			printf("Population has grown beyond limit. \nThe program will now exit.\n");
			fclose(energy_f);
			return energy/(i + 1);
		} 
		if (EXTINCTION) {
			printf("Walkers are all dead. \nThe program will now exit.\n");
			fclose(energy_f);
			return energy/(i + 1);
		} 
		
		printf("nw = %d\n\n", nw);
		fprintf(energy_f, "%lg\n", E[i]/N + energy_correction);
		
		if (i > M/10.)
			energy += E[i];
		
		if (nw < 0) {
			printf("ERROR: negative population.\n"); 
			fclose(energy_f);
			break;
		}
	}
	printf("Iteration successfully completed.\n");
	fclose(energy_f);
	printf("corr = %lg\n", energy_correction);
	return energy*10./9./M/N + energy_correction;
}


/**************************************************************************************************/

int main() 
{
	walker R[N_MAX];
	srand(321421);
	kill(R);
	init(R);
	printf("P = %d\n", P);
	printf("E = %g \n", DMC(R));
	return 0;
}











































