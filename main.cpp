#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <random>
#include <iostream>

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
	
	
	unsigned int n_ = 4;
	const double L = n_ * 1.7;
	double*** particles = initializeCube(n_, L);
	const unsigned int num = 4 * (int)pow(n_, 3);
	initializeVelocities(particles, num, 0, sqrt(300 / 119.7));
	
	
	/*
	const unsigned int num = 2;
	const double L = 3 * 1.7;
	double*** particles = (double***)calloc(2, sizeof(double*));
	particles[0] = createParticle();
	particles[0][2][0] = 1.5;
	particles[1] = createParticle();
	initializeVelocities(particles, num, 0, sqrt(300 / 119.7));
	*/
	
	FILE* dataFile = fopen("outputdata.txt", "w");				// Ovito-file
	FILE* potFile = fopen("outputdata_pot.txt", "w");			// Potential energy
	FILE* velFile = fopen("outputdata_vel.txt", "w");			// Velocities 
	fprintf(potFile, "%lf %d\n", endTime, timesteps);


	for (unsigned int i = 0; i < timesteps - 1; i++)
	{
		double dt = t[i + 1] - t[i];
		
		/*
		// Write to output files:
		writeToDataFile(&dataFile, particles, num);
		writeToPotFile(&potFile, particles, num, L);
		writeToVelFile(&velFile, particles, num, endTime, timesteps);
		*/

		// Update new position:
		for (unsigned int j = 0; j < num; j++)
		{
			for (int n = 0; n < 3; n++)
			{
				double new_pos = particles[j][2][n] + particles[j][1][n] * dt + 0.5 * particles[j][0][n] * pow(dt, 2);
				if (new_pos < 0) particles[j][2][n] = L + fmod(new_pos, L);
				else particles[j][2][n] = fmod(new_pos, L);

			}
		}
		
		double** accumulated = (double**)calloc(num, sizeof(double*));
		//if (accumulated == NULL) goto ErrorHandling;
		for (unsigned int j = 0; j < num; j++)
		{
			*(accumulated + j) = (double*)calloc(3, sizeof(double));
			//if (*(accumulated + j) == NULL) goto ErrorHandling;
		}

		for (unsigned int j = 0; j < num; j++)
		{
			// Accumulate accelleration:
			for (unsigned int k = j + 1; k < num; k++)
			{
				double dr[3] = { 0. };
				double distsq = 0.;
				for (unsigned int n = 0; n < 3; n++)
				{
					double distn = particles[j][2][n] - particles[k][2][n];
					distn -= round(distn / L) * L;
					distsq += pow(distn, 2);
					dr[n] = distn;
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
				particles[j][1][n] += 0.5 * (particles[j][0][n] + accumulated[j][n]) * dt;
				particles[j][0][n] = accumulated[j][n];
			}
		}

		// Write to output files:
		writeToDataFile(&dataFile, particles, num);
		writeToPotFile(&potFile, particles, num, L);
		writeToVelFile(&velFile, particles, num, endTime, timesteps);

		// Free allocated memory:
		for (unsigned int j = 0; j < num; j++)
		{
			free(accumulated[j]);
			accumulated[j] = NULL;
		}
		free(accumulated);
		accumulated = NULL;
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
/*
ErrorHandling:
	printf("Calloc returned NULL pointer in main. \n");
	return 1;
*/
	}

double** createParticle(void)
{
	double** p;
	p = (double**)calloc(3, sizeof(double*));
	//if (p == NULL) goto ErrorHandling;

	for (int i = 0; i < 3; i++)
	{
		*(p + i) = (double*)calloc(3, sizeof(double));
		//if (*(p + i) == NULL) goto ErrorHandling;
	}
	return p;
/*
ErrorHandling:
	printf("CreateParticle returned NULL pointer.\n");
	exit(1);'
*/
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
	double* p = NULL;
	p = (double*)malloc(steps * sizeof(double));
	//if (p == NULL) goto ErrorHandling;
	double dt = (end - start) / (steps-1);
	for (unsigned int i = 0; i < steps; i++)
	{
		*(p + i) = start + i * dt;
	}
	return p;
/*
ErrorHandling:
	printf("Linspace returned NULL pointer. \n");
	exit(1);
*/
	}

double*** initializeParticles(unsigned int num)
{
	double*** p;
	p = (double***)calloc(num, sizeof(double**));
	//if (p == NULL) goto ErrorHandling;
	for (unsigned int i = 0; i < num; i++)
	{
		*(p + i) = createParticle();
	}
	return p;
/*
ErrorHandling:
	printf("InitializeParticles returned NULL pointer.\n");
	exit(1);
*/
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
				setPos(p[index + 1], i*d, (0.5 + j)*d, (0.5+k)*d);
				setPos(p[index + 2], (0.5 + i)*d, j*d, (0.5 + k)*d);
				setPos(p[index + 3], (0.5 + i)*d, (0.5 + j)*d, k*d);
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
	double Z = sqrt(- 2 * log(U1)) * cos(2 * M_PI * U2);
	return mean + deviation * Z;
}

void initializeVelocities(double*** p, const int N, const double mean, const double deviation)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(mean, deviation);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			p[i][1][j] = distribution(generator);
		}
	}
}