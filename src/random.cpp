/*
 * random.cpp
 *
 *  Created on: Sep 8, 2010
 *      Author: Marzuk M Kamal
 *      Adapted from Numerical Recipes in C
 */

#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "random.h"
#include <unistd.h>

/************
 *
 * Random number generator
 *
 * call init_random() before using random_num() function.
 *
 * call delete_random(), at the end of the program
 *
 */


#define IM1 	2147483563
#define IM2 	2147483399
#define AM 		(1.0/IM1)
#define IMM1 	(IM1-1)
#define IA1 	40014
#define IA2 	40692
#define IQ1 	53668
#define IQ2 	52774
#define IR1 	12211
#define IR2 	3791
#define NTAB 	32
#define NDIV 	(1+IMM1/NTAB)
#define EPS 	1.2e-7
#define RNMX 	(1.0-EPS)


#define MAX_MEMORY_FOR_RANDOM_DATA		(1*1024) // KB
#define	MAX_RANDOM_DATA					(MAX_MEMORY_FOR_RANDOM_DATA*1024/8)


double *random_data = NULL;
int		data_index = 0;
long 	idummy = -1;


double ran2(long *idum){
	int j;
	long k;
	double temp;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];

	if (*idum <= 0){
		if (-(*idum) < 1)
			*idum=1;
		else
			*idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--){
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;

			if (*idum < 0)
				*idum += IM1;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;

	if (*idum < 0)
		*idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;

	if (idum2 < 0)
		idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;

	if (iy < 1)
		iy += IMM1;
	if ((temp=AM*iy) > RNMX)
		return (double) RNMX;
	else
		return temp;
}


void init_random(){

	if(random_data == NULL){
		idummy = -1;
		random_data = (double *) calloc(MAX_RANDOM_DATA,sizeof(double));

		if(random_data == NULL){
			printf("ERROR: Cannot allocate memory for random_data.\n");
			printf("Try reducing the size of MAX_RANDOM_DATA.\n");
			printf("Exiting program...\n");
			exit(-1);
		}

		ran2(&idummy);

		srand(time(NULL)*getpid());
	}

	data_index = rand() % MAX_RANDOM_DATA; // select the starting point

	for(int i = data_index+1; i < MAX_RANDOM_DATA; i++)
		random_data[i] = ran2(&idummy);

	for(int i = 0; i <= data_index; i++)
		random_data[i] = ran2(&idummy);

	data_index = 0;
}


void delete_random(){
	free(random_data);
	random_data = NULL;
	data_index = 0;
}

double random_num(){
	if(data_index >= MAX_RANDOM_DATA)
		init_random();

	return random_data[data_index++];
}

