/*
 * global_variables.h
 *
 *  Created on: Aug 23, 2010
 *      Author: kamal
 */

#ifndef GLOBAL_VARIABLES_H_
#define GLOBAL_VARIABLES_H_

extern int 		ITERATIONS; //number of iterations to optimize RNAs
extern int 		SUBITER; //frequency to plot in number of iterations
extern double 	TEMPERATURE; //initial temperature for Monte Carlo
extern double 	WEIGHT; //chances for kinetics instead of structure (related to weighting factor for multiobjective)
extern int 		WORDMUT; //number of numcleotides to mutate
extern double 	MUT_DIR; //probability to do a directed mutation instead of random
extern double 	DG_PARAM;
extern int	    SCORING_FUNCTION_MODE;
extern int		TOTAL_CPU_COUNT;
extern int		RNA_PROB_CALC_ENFORCE_FLAG;

extern double	RBS_FREE_THRESHOLD_VALUE;	// upper limit from which probabilistic calculation starts
extern double	RBS_BOUND_THRESHOLD_VALUE;


extern int		ARGC_VALUE;
extern char		**ARGV_VALUE;



extern double	PROBABILITY_SUM_WEIGHT;
extern double	CONCENTRATION_SUM_WEIGHT;

extern int		WAIT_FOR_PROBABILITY_MODE;
extern int		NUPACK_ENERGY_LOADED;


#endif /* GLOBAL_VARIABLES_H_ */
