/*
 * global_variables.cpp
 *
 *  Created on: Aug 23, 2010
 *      Author: kamal
 */

#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H


int 	ITERATIONS; //number of iterations to optimize RNAs
int 	SUBITER; //frequency to plot in number of iterations
double 	TEMPERATURE; //initial temperature for Monte Carlo
double 	WEIGHT; //chances for kinetics instead of structure (related to weighting factor for multiobjective)
int 	WORDMUT; //number of numcleotides to mutate
double 	MUT_DIR; //probability to do a directed mutation instead of random .. not used!
int	  	SCORING_FUNCTION_MODE = 1; // default is RNAgate scoring function
long	TOTAL_CPU_COUNT; // total number of processors in the computer.
int		RNA_PROB_CALC_ENFORCE_FLAG = 1;

double	RBS_FREE_THRESHOLD_VALUE = 0.001;	// upper limit from which probabilistic calculation starts
double	RBS_BOUND_THRESHOLD_VALUE = 0.8;

int		ARGC_VALUE;
char	**ARGV_VALUE;


double	PROBABILITY_SUM_WEIGHT = 100.0;
double	CONCENTRATION_SUM_WEIGHT = 100.0;

int		WAIT_FOR_PROBABILITY_MODE = 1;
int		NUPACK_ENERGY_LOADED = 0;

double DG_PARAM = 15;

#endif
