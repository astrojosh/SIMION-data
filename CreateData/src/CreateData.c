/*
 ============================================================================
 Name        : CreateData.c
 Author      : Josh Morrow
 Version     : 1.0
 Copyright   :
 Description : Create ion data for SIMION
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

double randomNum(int min, int max)
{
  return  min + (double) (rand() / (double) (RAND_MAX + 1) * (max - min));
}

int main(void) {

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	int num;

	printf("Number of Ions\n");
	scanf("%d", &num);

	int tob = 0,
		cwf = 1,
		colour = 1;

	double mass, charge, x, y, z;

	printf("Mass / amu\n");
	scanf("%lf", &mass);

	printf("Charge / e\n");
	scanf("%lf", &charge);

	printf("x / gu\n");
	scanf("%lf", &x);

	printf("y / gu\n");
	scanf("%lf", &y);

	printf("z / gu\n");
	scanf("%lf", &z);

	int i;
	double azimuth[num], elevation[num], kinetic_energy[num],
			kinetic_energy_x[num], kinetic_energy_y[num], kinetic_energy_z[num];

	FILE *finout;
	finout = fopen("data.ion","w");

	for(i = 0; i < num; i++) {

		azimuth[i] = randomNum(0, 360);
		elevation[i] = 90 - acos(randomNum(-1, 1)) * 180.0 / PI;

		kinetic_energy_x[i] = -0.02585*log(randomNum(0, 1));
		kinetic_energy_y[i] = -0.02585*log(randomNum(0, 1));
		kinetic_energy_z[i] = -0.02585*log(randomNum(0, 1));

		kinetic_energy[i] = kinetic_energy_x[i] + kinetic_energy_y[i] + kinetic_energy_z[i];

		fprintf(finout,"%d,%g,%g,%g,%g,%g,%f,%f,%f,%d,%d\n", tob, mass, charge, x, y, z, azimuth[i], elevation[i], kinetic_energy[i], cwf, colour);

	}

	fclose(finout);

	return 0;
}
