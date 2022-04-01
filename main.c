#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define _USE_MATH_DEFINES
#include <math.h>

#pragma warning(disable: 5045)			// Warning about Spectre mitigations
#pragma warning(disable: 4996)			// Warning about unsafe fopen function
#pragma warning(disable: 4100)			// Warning about unreferenced parameter (argv)


double** createParticle(void);
double* linspace(double start, double end, unsigned int steps);
double*** initializeParticles(unsigned int num);
double*** initializeCube(unsigned const int n, const double L);
void writeToDataFile(FILE** infile, double*** p, unsigned const int N);
void writeToVelFile(FILE** infile, double*** p, unsigned const int N, const double endTime, const int timesteps);
void writeToPotFile(FILE** infile, double*** p, unsigned const int N, const double L);
double normalDist(double mean, double stdDeviation);
void initializeVelocities(double*** p, const int N, const double mean, const double deviation);
void initializePreviousInstance(double*** p, const unsigned int N);
void rdf(double distsq, double binsize, int* bin, int maxBinIndex);



int main(int argc, char* argv)
{
	unsigned int timesteps = 3000;
	const double endTime = 3;
	double* t = linspace(0, endTime, timesteps);


	unsigned int n_ = 6;								// Unit-cells along each axis
	const double L = n_ * 1.7;							// Size of cube
	
	double*** particles;								// Declaration of particle-array
	const unsigned int num = 4 * (int)pow(n_, 3);		// Total number of particles

	double binSize = 0.005;
	static int maxBinIndex = 1000;
	int bin[1000] = { 0 };

	// Pass some random Command-line argument to initialize particles with last data from previous simulation
	if (argc > 1)
	{
		particles = initializeParticles(num);
		initializePreviousInstance(particles, num);
		printf("Initialized previous instance.\n");
	}
	else
	{
		particles = initializeCube(n_, L);
		initializeVelocities(particles, num, 0, sqrt(112));
	}
	
	FILE* dataFile = fopen("outputdata.txt", "w");				// Ovito-file
	FILE* potFile = fopen("outputdata_pot.txt", "w");			// Potential energy
	FILE* velFile = fopen("outputdata_vel.txt", "w");			// Velocities 
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
				if (new_pos < 0) new_pos = L + fmod(new_pos, L);
				else if (new_pos > L)  new_pos = fmod(new_pos, L);
				particles[j][2][n] = new_pos;

			}
		}
		
		// Allocate arrays to accumulate forces for each axis:
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
				
				for (int n = 0; n < 3; n++)
				{
					double distn = particles[j][2][n] - particles[k][2][n];
					distn -= round(distn / L) * L;
					distsq += pow(distn, 2);
					dr[n] = distn;
				}

				if (i == timesteps - 2) rdf(distsq, binSize, bin, maxBinIndex);


				if (distsq < 9.)
				{
					for (int n = 0; n < 3; n++)
					{
						double acc_new = 24. * (2. * pow(distsq, -6) - pow(distsq, -3)) * dr[n] /distsq;
						accumulated[j][n] += acc_new;
						accumulated[k][n] -= acc_new;
					}
				}
			}
		}

		// Update accelleration + velocity:
		for (unsigned int j = 0; j < num; j++)
		{
			for (int n = 0; n < 3; n++)
			{
				particles[j][1][n] += 0.5f * (particles[j][0][n] + accumulated[j][n]) * dt;
				particles[j][0][n] = accumulated[j][n];
			}
		}

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
 


	FILE* lastPos = fopen("outputdata_lastPos.txt", "w");
	FILE* lastVel = fopen("outputdata_lastVel.txt", "w");

	fprintf(lastPos, "%u\n", num);						// Provide number of atoms to data file

	// Write last position and velocity to files:
	for (unsigned int i = 0; i < num; i++)
	{
		fprintf(lastPos, "%lf %lf %lf\n", particles[i][2][0], particles[i][2][1], particles[i][2][2]);
		fprintf(lastVel, "%lf %lf %lf\n", particles[i][1][0], particles[i][1][1], particles[i][1][2]);
	}

	// Write normalized bin values from last iteration in main to file:
	FILE* rdfFile = fopen("outputdata_rdf.txt", "w");
	for (int i = 0; i < maxBinIndex; i++)
	{
		double normalized = pow(L, 3) / pow(num, 2) * bin[i] / (4 * M_PI * pow(0.5 * binSize + i * binSize, 2) * binSize);
		fprintf(rdfFile, "%lf\n", normalized);
	}	


	fclose(lastPos);
	fclose(lastVel);
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
	p = (double**)calloc(3, sizeof(double*));
	if (p == NULL) goto ErrorHandling;

	for (int i = 0; i < 3; i++)
	{
		*(p + i) = (double*)calloc(3, sizeof(double));
		if (*(p + i) == NULL) goto ErrorHandling;
	}
	return p;

ErrorHandling:
	printf("CreateParticle returned NULL pointer.\n");
	exit(1);
}

double* linspace(double start, double end, unsigned int steps)
{
	double* p = NULL;
	p = (double*)malloc(steps * sizeof(double));
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
	p = (double***)calloc(num, sizeof(double**));
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
				p[index][2][0]   = i*d;			p[index][2][1]   = j*d;			p[index][2][2]	 = k*d;
				p[index+1][2][0] = i*d;			p[index+1][2][1] = (0.5+j)*d;	p[index+1][2][2] = (0.5+k)*d;
				p[index+2][2][0] = (0.5+i)*d;	p[index+2][2][1] = j*d;			p[index+2][2][2] = (0.5+k)*d;
				p[index+3][2][0] = (0.5+i)*d;	p[index+3][2][1] = (0.5+j)*d;	p[index+3][2][2] = k*d;
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
				dn -= round(dn / L) * L;
				drSq += pow(dn, 2);
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
	srand(clock());
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			p[i][1][j] = normalDist(mean, deviation);
		}
	}
}

void initializePreviousInstance(double*** p, const unsigned int N)
{
# pragma warning(disable: 6031)			// Warning about sscanf return value ignored

	FILE* dataFile = fopen("outputdata_lastPos.txt", "r");
	FILE* velFile = fopen("outputdata_lastVel.txt", "r");
	if (dataFile == NULL || velFile == NULL)
	{
		printf("Opening datafiles returned NULL pointer, run a simulation first.\n");
		exit(1);
	}
	char buffer[128];
	if (fgets(buffer, 128, dataFile) == NULL)
	{
		printf("Failed to read header from dataFile.\n");
		exit(1);
	}
	unsigned int prevNum;
	sscanf(buffer, "%u", &prevNum);
	if (prevNum != N)
	{
		printf("Atoms from previous simulation does not match current simulation.\nPrev: %u\nCurr: %u\n", prevNum, N);
		exit(1);
	}
	for (unsigned int i = 0; i < prevNum; i++)
	{
		if (fgets(buffer, 128, dataFile) == NULL)
		{
			printf("Failed to read header from dataFile.\n");
			exit(1);
		}
		sscanf(buffer, "%lf %lf %lf", &p[i][2][0], &p[i][2][1], &p[i][2][2]);

		if (fgets(buffer, 128, velFile) == NULL)
		{
			printf("Failed to read header from dataFile.\n");
			exit(1);
		}
		sscanf(buffer, "%lf %lf %lf", &p[i][1][0], &p[i][1][1], &p[i][1][2]);
	}

	fclose(dataFile);
	fclose(velFile);
}

void rdf(double distsq, double binSize, int* bin, int maxBinIndex)
{
	double dist = sqrt(distsq);
	int binIndex = floor(dist / binSize);
	if (binIndex < maxBinIndex) bin[binIndex] ++;
}