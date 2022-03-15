#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

#pragma warning(disable: 5045)			// Warning about Spectre mitigations
#pragma warning(disable: 4996)			// Warning about unsafe fopen function


double** createParticle(void);
void setAcc(double** data, double a, double b, double c);
void setVel(double** data, double a, double b, double c);
void setPos(double** data, double a, double b, double c);
double* linspace(double start, double end, unsigned int steps);
double*** initializeParticles(unsigned int num);
double*** initializeCube(unsigned const int n, const double L);
void writeToDataFile(FILE** infile, double*** p, unsigned const int N);
void writeToVelFile(FILE** infile, double*** p, unsigned const int N, const double endTime, const int timesteps);
void writeToPotFile(FILE** infile, double*** p, unsigned const int N, const double L);
double normalDist(double mean, double stdDeviation);
void initializeVelocities(double*** p, const int N, const double mean, const double deviation);



int main(void)
{
	unsigned int timesteps = 5000;
	const double endTime = 5;
	double* t = linspace(0, endTime, timesteps);
	unsigned int n_ = 2;
	const double L = n_ * 1.7f;
	double*** particles = initializeCube(n_, L);
	//particles[0][1][0] = 20.f;

	unsigned const int num = 4 * (int)pow(n_, 3);

	initializeVelocities(particles, num, 0, sqrt(300));

	
	FILE* dataFile = fopen("outputdata.txt", "w");				// Ovito-file
	FILE* potFile = fopen("outputdata_pot.txt", "w");			// Potential energy
	FILE* velFile = fopen("outputdata_vel.txt", "w");			// Velocities squared
	fprintf(potFile, "%lf %d\n", endTime, timesteps);


	for (unsigned int i = 0; i < timesteps - 1; i++)
	{
		double dt = t[i + 1] - t[i];

		// Write to output files:
		writeToDataFile(&dataFile, particles, num);
		writeToPotFile(&potFile, particles, num, L);
		writeToVelFile(&velFile, particles, num, endTime, timesteps);

		// Update new position:
		for (unsigned int j = 0; j < num; j++)
		{
			for (int n = 0; n < 3; n++)
			{
				double new_pos = particles[j][2][n] + particles[j][1][n] * dt + 0.5f * particles[j][0][n] * pow(dt, 2);
				//if (new_pos < 0) new_pos = L + fmod(new_pos, L);
				//else if (new_pos > L) new_pos = fmod(new_pos, L);
				particles[j][2][n] = fmod(new_pos, L);
			}
		}
		
		double** accumulated = (double**)calloc(num, sizeof(double*));
		if (accumulated == NULL) goto ErrorHandling;
		for (unsigned int j = 0; j < num; j++)
		{
			*(accumulated + j) = (double*)calloc(3, sizeof(double));
			if (*(accumulated + j) == NULL) goto ErrorHandling;
		}

		for (unsigned int j = 0; j < num; j++)
		{
			// Accumulate accelleration:
			for (unsigned int k = j + 1; k < num; k++)
			{
				double dr[3] = { 0. };
				double distsq = 0.;
				for (unsigned int d = 0; d < 3; d++)
				{
					dr[d] = particles[j][2][d] - particles[k][2][d];
					dr[d] -= round(dr[d] / L) * L;
					distsq += pow(dr[d], 2);
				}

				if (distsq < 9.)
				{
					for (int n = 0; n < 3; n++)
					{
							double acc_new = 24. * (2. * pow(distsq, -6) - pow(distsq, -3)) * dr[n] / distsq;
							accumulated[j][n] += acc_new;
							accumulated[k][n] -= acc_new;
					}
				}
			}
			
			// Update accelleration + velocity:
			for (int n = 0; n < 3; n++)
			{
				particles[j][0][n] = accumulated[j][n];
				particles[j][1][n] += 0.5f * (particles[j][0][n] + accumulated[j][n]) * dt;
			}
		}

		// De-allocate allocated memory:
		accumulated = NULL;
		free(accumulated);
	}

	writeToDataFile(&dataFile, particles, num);
	writeToPotFile(&potFile, particles, num, L);
	writeToVelFile(&velFile, particles, num, endTime, timesteps);

	fclose(dataFile);
	fclose(velFile);
	fclose(potFile);
	free(t);
	free(particles);
	return 0;

ErrorHandling:
	printf("Calloc returned NULL pointer in main. \n");
	return 1;
}

double** createParticle(void)
{
	double** p;
	p = calloc(3, sizeof(double*));
	if (p == NULL) goto ErrorHandling;

	for (int i = 0; i < 3; i++)
	{
		*(p + i) = calloc(3, sizeof(double));
		if (*(p + i) == NULL) goto ErrorHandling;
	}
	return p;

ErrorHandling:
	printf("CreateParticle returned NULL pointer.\n");
	exit(1);
}

void setAcc(double** data, double a, double b, double c)
{
	data[0][0] = a;
	data[0][1] = b;
	data[0][2] = c;
}

void setVel(double** data, double a, double b, double c)
{
	data[1][0] = a;
	data[1][1] = b;
	data[1][2] = c;
}

void setPos(double** data, double a, double b, double c)
{
	data[2][0] = a;
	data[2][1] = b;
	data[2][2] = c;
}

double* linspace(double start, double end, unsigned int steps)
{
	double* p;
	p = malloc(steps * sizeof(double));
	if (p == NULL) goto ErrorHandling;
	double dt = (end - start) / (steps-1);
	for (unsigned int i = 0; i < steps; i++)
	{
		*(p + i) = start + i * dt;
	}
	return p;

ErrorHandling:
	printf("Linspace returned NULL pointer. \n");
	exit(1);
}

double*** initializeParticles(unsigned int num)
{
	double*** p;
	p = calloc(num, sizeof(double**));
	if (p == NULL) goto ErrorHandling;
	for (unsigned int i = 0; i < num; i++)
	{
		*(p + i) = createParticle();
	}
	return p;

ErrorHandling:
	printf("InitializeParticles returned NULL pointer.\n");
	exit(1);
}

double*** initializeCube(unsigned const int n, const double L)
{
	double*** p = initializeParticles(4 * (int)pow(n, 3));
	double d = L / n;
	int index = 0;
	for (unsigned int i = 0; i < n; i++)
	{
		for (unsigned int j = 0; j < n; j++)
		{
			for (unsigned int k = 0; k < n; k++)
			{
				setPos(p[index], i*d, j*d, k*d);
				setPos(p[index + 1], i*d, (0.5f + j)*d, (0.5+k)*d);
				setPos(p[index + 2], (0.5f + i)*d, j*d, (0.5f + k)*d);
				setPos(p[index + 3], (0.5f + i)*d, (0.5f + j)*d, k*d);
				index += 4;
			}
		}
	}
	return p;
}


void writeToDataFile(FILE** infile, double*** p, unsigned const int N)
{
	fprintf(*infile, "%d\n", N);
	fprintf(*infile, "type x y z\n");

	for (unsigned int j = 0; j < N; j++)
	{
		fprintf(*infile, "Ar %lf %lf %lf\n", p[j][2][0], p[j][2][1], p[j][2][2]);
	}
}

void writeToVelFile(FILE** infile, double*** p, unsigned const int N, const double endTime, const int timesteps)
{
	fprintf(*infile, "*%d %lf %d\n", N, endTime, timesteps);

	for (unsigned int j = 0; j < N; j++)
	{
		fprintf(*infile, "%lf %lf %lf\n", p[j][1][0], p[j][1][1], p[j][1][2]);
	}
}

void writeToPotFile(FILE** infile, double*** p, unsigned const int N, const double L)
{
	double U_tot = 0;
	for (unsigned int j = 0; j < N; j++)
	{
		for (unsigned int k = j + 1; k < N; k++)
		{
			double drSq = 0;
			for (int n = 0; n < 3; n++)
			{
				double dn = p[j][2][n] - p[k][2][n];
				dn = pow(dn - round(dn / L) * L, 2);
				drSq += dn;
			}
			U_tot += 4 * (pow(drSq, -6) - pow(drSq, -3));
		}

	}
	fprintf(*infile, "%lf\n", U_tot);
}

double normalDist(const double mean, const double deviation)
{
	// Following the Box-Miller transformation
	double U1 = ((double) rand() / (RAND_MAX));
	double U2 = ((double) rand() / (RAND_MAX));
	return mean + deviation * sqrt(-2 * log(U1)) * cos(2 * M_PI * U2);
}

void initializeVelocities(double*** p, const int N, const double mean, const double deviation)
{
	srand((unsigned int)clock());
	for (int i = 0; i < N; i++)
	{
		setVel(p[i], normalDist(mean, deviation), normalDist(mean, deviation), normalDist(mean, deviation));
	}
}