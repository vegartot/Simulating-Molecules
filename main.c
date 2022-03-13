#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#pragma warning(disable: 5045)			// Warning about Spectre mitigations
#pragma warning(disable: 4996)			// Warning about unsafe fopen function


double** createParticle(void);
void setAcc(double** data, double a, double b, double c);
void setVel(double** data, double a, double b, double c);
void setPos(double** data, double a, double b, double c);
void pprint(double** data);
double* linspace(double start, double end, unsigned int steps);
double distance(double** p1, double** p2);
double*** initializeParticles(unsigned int num);
double*** initializeCube(unsigned const int n, const double L);


int main(void)
{
	unsigned int timesteps = 5000;
	double* t = linspace(0, 5, timesteps);
	unsigned int n_ = 4;
	const double L = n_ * 1.7f;

	double*** particles = initializeCube(n_, L);
	
	unsigned const int num = 4 * (int)pow(n_, 3);

	
	FILE* infile = fopen("outputdata.txt", "w");


	for (unsigned int i = 0; i < timesteps - 1; i++)
	{
		double dt = t[i + 1] - t[i];

		// Update new position:
		for (unsigned int j = 0; j < num; j++)
		{
			for (int n = 0; n < 3; n++)
			{
				double new_pos = particles[j][2][n] + particles[j][1][n] * dt + 0.5f * particles[j][0][n] * pow(dt, 2);
				
				if (new_pos < 0) new_pos += L;
				else if (new_pos > L) new_pos -= L;
				
				particles[j][2][n] = new_pos;
			}
		}
		
		double** accumulated = calloc(num, sizeof(double*));
		if (accumulated == NULL) goto ErrorHandling;
		for (unsigned int j = 0; j < num; j++)
		{
			*(accumulated + j) = calloc(3, sizeof(double));
			if (*(accumulated + j) == NULL) goto ErrorHandling;
		}


		for (unsigned int j = 0; j < num; j++)
		{

			// Accumulate forces:
			for (unsigned int k = j + 1; k < num; k++)
			{
				double dr = distance(particles[j], particles[k]);
				//dr -= ceil(dr / L) * L;

				{
				if (dr < 3)
					for (int n = 0; n < 3; n++)
					{
						double direction = particles[j][2][n] - particles[k][2][n];
						double acc_new = 24 * (2 * pow(dr, -12) - pow(dr, -6)) * direction / pow(dr, 2);
						accumulated[j][n] += acc_new;
						accumulated[k][n] -= acc_new;
					}
				}
			}

			// Update velocity:
			for (int n = 0; n < 3; n++)
			{
				particles[j][1][n] += 0.5f * (particles[j][0][n] + accumulated[j][n]) * dt;
			}

			// Update acceleration:
			particles[j][0] = accumulated[j];
		}
		// De-allocate allocated memory:
		accumulated = NULL;
		free(accumulated);

		fprintf(infile, "%d\n", num);
		fprintf(infile, "type x y z\n");
		for (unsigned int j = 0; j < num; j++)
		{
			fprintf(infile, "Ar %lf %lf %lf\n", particles[j][2][0], particles[j][2][1], particles[j][2][2]);
		}
	}
	fclose(infile);
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

void pprint(double** data)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf("%lf ", data[i][j]);
		}
		printf("\n");
	}
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

double distance(double** p1, double** p2)
{
	double d[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		d[i] = pow(fabs(p1[2][i] - p2[2][i]), 2);
	}
	return sqrt(d[0] + d[1] + d[2]);
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

void initializeVelocities(double*** p)
{

}