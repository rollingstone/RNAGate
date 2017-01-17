/*
 * permutation_combination.cpp
 *
 *  Created on: Oct 8, 2010
 *      Author: kamal
 */


#include "permutation_combination.h"


int fact(int n){

	if(n == 0)
		return 1;

	if(n < 0)
		return -1;

	int	fact_n = 1;

	while(n)
		fact_n *= n--;

	return fact_n;
}


int combination(int n, int r){
	if(n < r)
		return 1;

	return fact(n)/(fact(n-r) * fact(r));
}

int permutation(int n, int r){
	if(n < r)
		return 1;

	return fact(n)/(fact(n-r));
}
