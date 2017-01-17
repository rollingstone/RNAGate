/*
 * random.h
 *
 *  Created on: Sep 8, 2010
 *      Author: kamal
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include<cstdlib>
#include<cstdio>

#ifdef __cplusplus
extern "C"{
#endif

double 	ran2(long *idum);
void 	init_random();
void 	delete_random();
double 	random_num();


#ifdef __cplusplus
}
#endif

#endif /* RANDOM_H_ */
