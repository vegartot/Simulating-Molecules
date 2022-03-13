#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#pragma warning(disable: 5045)

double** createParticle(void);
void setAcc(double** data, double a, double b, double c);
void setVel(double** data, double a, double b, double c);
void setPos(double** data, double a, double b, double c);
void pprint(double** data);
typedef double** particle;

double* linspace(double start, double end, unsigned int steps);
double distance(particle p1, particle p2);


void main(void)
{
	unsigned int timesteps = 5000;
	double* t = linspace(0, 5, timesteps);
	
	particle p1 = createParticle();
	particle p2 = createParticle();
	setPos(p2, 1.5, 0, 0);
	setVel(p2, 0, 0, 0);

#pragma warning(suppress: 4996)
	FILE* infile = fopen("outputdata.txt", "w");

	for (unsigned int i = 0; i < timesteps-1; i++)
	{
		double dr = distance(p1, p2);
		double dt = t[i + 1] - t[i];

		double p1_new_pos[3];
		double p1_new_acc[3];
		double p1_new_vel[3];

		double p2_new_pos[3];
		double p2_new_acc[3];
		double p2_new_vel[3];

		for (int j = 0; j < 3; j++)
		{
			double dir = p1[2][j] - p2[2][j];

			p1_new_pos[j] = p1[2][j] + p1[1][j] * dt + 0.5f * p1[0][j] * pow(dt, 2);
			p1_new_acc[j] = 24.f * (2 * pow(dr, -12) - pow(dr, -6)) * dir / pow(dr, 2);
			p1_new_vel[j] = p1[1][j] + 0.5f * (p1[0][j] + p1_new_acc[j]) * dt;

			p2_new_pos[j] = p2_new_pos[j] = p2[2][j] + p2[1][j] * dt + 0.5f * p2[0][j] * pow(dt, 2);
			p2_new_acc[j] = - p1_new_acc[j];
			p2_new_vel[j] = p2[1][j] + 0.5f * (p2[0][j] + p2_new_acc[j]) * dt;
		}
		p1[0] = p1_new_acc;
		p1[1] = p1_new_vel;
		p1[2] = p1_new_pos;

		p2[0] = p2_new_acc;
		p2[1] = p2_new_vel;
		p2[2] = p2_new_pos;

		fprintf(infile, "2\n");
		fprintf(infile, "type x y z\n");
		fprintf(infile, "Ar %lf %lf %lf\n", p1[2][0], p1[2][1], p1[2][1]);
		fprintf(infile, "Ar %lf %lf %lf\n", p2[2][0], p2[2][1], p2[2][2]);

	}
	fclose(infile);
}

double** createParticle(void)
{
	double** p;
	p = calloc(3, sizeof(double*));
	if (p == NULL) goto HandleError;

	for (int i = 0; i < 3; i++)
	{
		*(p + i) = calloc(3, sizeof(double));
		if (*(p + i) == NULL) goto HandleError;
	}
	return p;

HandleError:
	printf("Calloc returned NULL pointer.\n");
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
	printf("Malloc returned NULL pointer. \n");
	exit(1);
}

double distance(particle p1, particle p2)
{
	double d[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		d[i] = pow(fabs(p1[2][i] - p2[2][i]), 2);
	}
	return sqrt(d[0] + d[1] + d[2]);
}
