/*
 * nupack_extension.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: kamal
 */

/*
  CalcConc.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  CALCCONC.C

  For use with Concentrations.c.

  Computes the equilibrium mole fractions of the products of
  aggregation reactions in dilute solution, given the identity of the
  aggregates (called "complexes") and their repsective free energies.

  This program solves the problem presented in Dirks, et al.,
  "Thermodynamic analysis of interacting nucleic acid strands", SIAM
  Review, in press.  All variable names are chosen to match those in
  that paper.

  This file contains CalcConc, the main function that computes the
  equilibrium mole fractions of the complexes, and the auxillary
  functions it calls.  Some of the functions it calls are standard
  utility functions, such as functions to sum entries in an array,
  etc.  These are included in utils.c.

  The trust region algorithm for solving the dual problem is that in
  Nocedal and Wright, Numerical Optimization, 1999, page 68, with the
  dogleg method on page 71.  There are some inherent precision issues.
  For most systems involving nucleic acids, this isn't a problem, but
  for some problematic cases, adjustments are made.  For some initial
  conditions, these precision issues cannot be overcome and a new
  initial condition must be generated.  This is done by randomly
  perturbing the standard initial condition (see comments in the
  function getInitialCondition, below), and re-running the trust
  region optimization.

  The inputs are as follows:
  A: A 2-D array; A[i][j] is the number of monomers of type i in complex j
  G: An array containing the corresponding complex free energies in units of kT.
     G[j] corresponds to the entries A[..][j].
  x0: Initial mole fractions of the unit-size complexes as mole fractions.
  numSS: The number of single-species (monomers) in the system.
  numTotal: The total number of complexes.
  maxIters: Maximum number of interations allowed in trust region method.
  tol: The tolerance for convergence.  The absolute tolerance is tol*(mininium
       single-species initial mole fraction)
  deltaBar: The maximum step size allowed in the trust region method
  eta: The value for eta in the trust region method, 0 < eta < 1/4
  outputFile: The file to which the output is written.
  kT: kT in kcal/mol.
  MaxNoStep: The maximum number of iterations without a step being taken before
             the initial conditions are regenerated.
  MaxTrial: The maximum number of initial conditions to be tried.
  PerturbScale: The multiplier on the random perturbations to the initial conditions
                as new ones are generated.
  quiet: = 1 for no printing of messages (except error messages) to screen.
  WriteLogFile: = 1 if log file is to be written
  logFile: file for printing information about the run.
  MolesWaterPerLiter: Number of moles of water per liter
*/


#include "nupack_conc_object.h"
#include<fstream>
#include<iostream>
#include <iomanip>

using namespace std;
extern "C"{
#include "ConcentrationsHeaderFile.h" // Header file for concentrations
}


Concentration::Concentration(){
	TolarenceValue = 0.0000001; // in percentage
	DeltaBarValue = TRUST_REGION_DELTABAR; // kT
	TemperatureKTValue = kB*(37.0 + ZERO_C_IN_KELVIN);

	MaxIterationValue = 1000*(1+10*1);
	MaxTrialValue = 100*(1+10*1);

	DoNotShowErrorMessage = 0;
	MaxNoStepValue = 50;
	MolesWaterPerLiterValue = 1;
	PerturbScaleValue = 100;
	ETA_Value = TRUST_REGION_ETA;

	calculation_flag = false;
	is_initialized = false;
}

Concentration::Concentration(int numSS, int numTotal){
	init_object(numSS,numTotal);
}

Concentration::~Concentration(){
	cout << "Destructor ~Concentration called..." << endl;
	free_object();
}


void Concentration::init_parameters(int total_rna, int total_complex,double *rna_conc, int **A1,  int *complex_index,int main_cpos){

	init_object(total_rna,total_complex);

	A_matrix_main = A1;
	rna_concentration_vector = rna_conc;
	main_complex_position = main_cpos;
	complex_index_vector = complex_index;

	for(int i = 0; i < TotalNumberOfComplex; i++){
		complex_ptr[complex_index_vector[i]] = i;
	}

}



void Concentration::print_A_matrix(){

	cout << "AMatrix:\n";
	for(int y = 0; y < NumberOfSingleStrand; y++){
		for(int x = 0; x < TotalNumberOfComplex; x++){
			cout << A_matrix_main[y][x] << " ";
		}
		cout << endl;
	}
}

void Concentration::print_text(char *t1){
	printf("\nThis is a test:\t%s",t1);
}

void Concentration::init_object(int numSS, int numTotal){
/*
 * Allocates memory for the calculation;
 *
 *
 */
	TolarenceValue = 1e-10;//0.000000001; // in percentage
	DeltaBarValue = TRUST_REGION_DELTABAR;///3; // kT
//	TemperatureKTValue = 0.6; // 37C in kT

	TemperatureKTValue = kB*(37.0 + ZERO_C_IN_KELVIN);

	MaxIterationValue = 1000*10;
	MaxTrialValue = 100*10;

	DoNotShowErrorMessage = 1;
	MaxNoStepValue = 50;
	MolesWaterPerLiterValue = 1;
	PerturbScaleValue = 100;
	ETA_Value = TRUST_REGION_ETA;

	calculation_flag = false;
	is_initialized = false;


	AT = (int **) malloc(numTotal * sizeof(int *));
	for(int j = 0; j < numTotal; j++){
		AT[j] = (int *) malloc(numSS * sizeof(int));
	}
	Hes = (double **) malloc(numSS * sizeof(double *));
	for(int i = 0; i < numSS; i++){
		Hes[i] = (double *) malloc(numSS * sizeof(double));
	}
	AbsTol = (double *) malloc(numSS * sizeof(double));
	Grad = (double *) malloc(numSS * sizeof(double));
	lambda = (double *) malloc(numSS * sizeof(double));
	p = (double *) malloc(numSS * sizeof(double));
	idum = (long *) malloc (sizeof(long));

	NumberOfSingleStrand = numSS;
	TotalNumberOfComplex = numTotal;


	rna_complex_concentration = (double *) malloc(numTotal * sizeof(double));
	rna_complex_free_energy   = (double *) malloc(numTotal * sizeof(double));

	is_initialized = true;
}

void Concentration::free_object(){
/*
 * free up memory;
 */

	is_initialized = false;
}


bool Concentration::calculate_concentration(){

	calculation_flag = CalcConc(rna_complex_concentration,A_matrix_main,rna_complex_free_energy,rna_concentration_vector, NumberOfSingleStrand, TotalNumberOfComplex,
			MaxIterationValue, TolarenceValue, DeltaBarValue, ETA_Value, TemperatureKTValue,
		     MaxNoStepValue, MaxTrialValue, PerturbScaleValue, DoNotShowErrorMessage,
		     0, NULL, 55.12, clock());

	return calculation_flag;
}


bool Concentration::calculate_concentration(double *rna_conc,double *complex_conc){
	return (calculation_flag = CalcConc(complex_conc,A_matrix_main,rna_complex_free_energy,rna_conc, NumberOfSingleStrand, TotalNumberOfComplex,
			MaxIterationValue, TolarenceValue, DeltaBarValue, ETA_Value, TemperatureKTValue,
		     MaxNoStepValue, MaxTrialValue, PerturbScaleValue, DoNotShowErrorMessage,
		     0, NULL, 55.12, clock()) );
}

bool Concentration::calculate_concentration(double *complex_free_energy){
	return (calculation_flag = calculate_concentration(rna_complex_concentration,A_matrix_main,complex_free_energy,rna_concentration_vector));
}


//bool Concentration::calculate_concentration(double *complex_free_energy) {
//	return calculate_concentration(rna_complex_concentration,A_matrix_main,complex_free_energy,rna_concentration_vector);
//}


void Concentration::calculate_combinatorial_concentration(double rna_conc_min,double rna_conc_max, double min_rna_conc, ofstream &fconc){
	double	conc_vector[3000];
	double	complex_cvector[3000];
	double	rmin = 999999;
	int		chunk = 2;

	if(NumberOfSingleStrand == 2){
		fconc << "RNA0\t" <<"RNA1\t" << "R1/R0\t" <<"CplxConc\t" << "CplxCfrac\t"<<"flag" << endl;
//#pragma omp parallel for private(x,y)
		for(double x = rna_conc_min; x <= rna_conc_max; x++){
			conc_vector[0] = x;
			for(double y = rna_conc_min; y <= rna_conc_max; y++){
				conc_vector[1] = y;

				rmin = min(x,y);

				bool flag = calculate_concentration(conc_vector,complex_cvector);
				fconc << setprecision(3) << x << '\t' << y << '\t' << y/x << '\t' << complex_cvector[main_complex_position] << '\t'
					  << complex_cvector[main_complex_position]/rmin << '\t' << flag << endl;
			}
		}
	}
	else if(NumberOfSingleStrand == 3){
		fconc << "RNA0\t" <<"RNA1\t" <<  "RNA2\t" << "R1/R0\t" << "R2/R0\t" << "CplxConc\t" << "CplxCfrac\t" << "flag" << endl;

		for(double x = rna_conc_min; x <= rna_conc_max; x++){
			conc_vector[0] = x;
			for(double y = rna_conc_min; y <= rna_conc_max; y++){
				conc_vector[1] = y;
				for(double z = rna_conc_min; z <= rna_conc_max; z++){
					conc_vector[2] = z;

					rmin = min(min(x,y),z);

					bool flag = calculate_concentration(conc_vector,complex_cvector);
					fconc << setprecision(3) << x << '\t' << y << '\t' << z << '\t'
						  << y/x << '\t' << z/x << '\t'
						  << complex_cvector[main_complex_position] << '\t'
						  << complex_cvector[main_complex_position]/rmin << '\t'
						  << flag << endl;
				}
			}
		}
	}
}

bool Concentration::calculate_concentration(double *x, int **A, double *G, double *x0) {
	/*
    Computes the equilibrium mole fractions of species in dilute
    solution using a trust region algorithm on the dual problem.
    Discussion of the method is in Dirks, et al., Thermodynamic
    analysis of interacting nucleic acid strands, SIAM Review, (2006),
    in press.  The trust region algorithm for solving the dual problem
    is that in Nocedal and Wright, Numerical Optimization, 1999, page
    68, with the dogleg method on page 71.

    Returns true if converged and false otherwise.
	 */
/*
	  A: A 2-D array; A[i][j] is the number of monomers of type i in complex j
	  G: An array containing the corresponding complex free energies in units of kT.
	     G[j] corresponds to the entries A[..][j].
	  x0: Initial mole fractions of the unit-size complexes as mole fractions.
*/
	int i,j; // Counters, i is over single-species and j is over all complexes
	int iters; // Number of iterations

	int 	nNoStep; // Number of iterations without taking a step
	int 	nTrial; // Number of times we've perturbed lambda
	double 	FreeEnergy; // The free energy of the solution
	double 	tol;
	double 	deltaBar;
	double 	kT;
	int 	MaxTrial;
	int 	quiet;
	int 	MaxIters;
	int		MaxNoStep;
	double 	MolesWaterPerLiter;
	double 	eta;
	long 	seed;
	double 	PerturbScale;

	unsigned long 	rand_seed = 0;

	int 	numSS;
	int		numTotal;

	if(is_initialized == false){
		printf("\nERROR: Function calculate_concentration() is not initialized.\nCall init_conc_data() first...");
		exit(-1);
	}

	tol = 		TolarenceValue;
	deltaBar = 	DeltaBarValue;
	kT = 		TemperatureKTValue;
	MaxTrial = 	MaxTrialValue;
	MaxNoStep =	MaxNoStepValue;
	MaxIters =	MaxIterationValue;
	MolesWaterPerLiter = MolesWaterPerLiterValue;
	numSS = 	NumberOfSingleStrand;
	numTotal = 	TotalNumberOfComplex;
	eta	=		ETA_Value;
	PerturbScale  = PerturbScaleValue;

	quiet = 	DoNotShowErrorMessage;
	seed = 		(long) clock();


	return CalcConc(x, A, G, x0, numSS, numTotal,
		      MaxIters, tol, deltaBar, eta, kT,
		      MaxNoStep, MaxTrial, PerturbScale, quiet,
		     0, NULL,  MolesWaterPerLiter, seed);


//
//	// Initialize iters just so compiler doesn't give a warning when optimization is on
//	iters = 0;
//
//	// The absolute tolerance is a percentage of the entries in x0
//	for (i = 0; i < numSS; i++) {
//		AbsTol[i] = tol * x0[i];
//	}
//
//	// Compute AT (transpose of A), useful to have around.
//	IntTranspose(AT,A,numSS,numTotal);
//
//	nTrial = 0;
//	for (i = 0; i < numSS; i++) {
//		Grad[i] = AbsTol[i] + 1.0; // Initialize just to get started.
//	}
//	while (CheckTol(Grad,AbsTol,numSS) == 0 && nTrial < MaxTrial) {
//
//		if (nTrial == 1) {
//			// Seed the random number generator if necessary
//			rand_seed = GetRandSeed(seed);
//			init_genrand(rand_seed);
//		}
//
//		// Set initial guess
//		getInitialGuess(x0,lambda,G,AT,A,numSS,numTotal,PerturbScale,rand_seed);
//
//		// Calculate the counts of the species based on lambda
//		if (getx(x,lambda,G,AT,numSS,numTotal) == 0) { // Should be fine; checked prev.
//			if (quiet == 0) {
//				printf("Overflow error in calcution of mole fractions.\n\n");
//				printf("Exiting....\n");
//			}
//			exit(ERR_OVERFLOW);
//		}
//
//		// Calculate the gradient
//		getGrad(Grad,x0,x,A,numSS,numTotal);
//
//		// Initialize delta to be just less than deltaBar
//		delta = 0.99 * deltaBar;
//
//		// Initializations
//		iters = 0;
//		nNoStep = 0;
//		RunStats[0] = 0; // Number of pure Newton steps (didn't hit trust region boundary)
//		RunStats[1] = 0; // Number of pure Cauchy steps (hit trust region boundary)
//		RunStats[2] = 0; // Number of dogleg steps (part Newton and part Cauchy)
//		RunStats[3] = 0; // Number of steps with Cholesky failure forcing Cauchy step
//		RunStats[4] = 0; // Number of steps with irrelovent Cholesky failures
//		RunStats[5] = 0; // Number of failed dogleg calculations
//
//		// Run trust region with these initial conditions
//		while (iters < MaxIters && CheckTol(Grad,AbsTol,numSS) == 0
//				&& nNoStep < MaxNoStep) {
//
//			// Compute the Hessian (symmetric, positive, positive definite)
//			getHes(Hes,x,A,numSS,numTotal);
//
//			// Solve for the search direction
//			(RunStats[getSearchDir(p,Grad,Hes,delta,numSS) - 1])++;
//
//			// Calculate rho, ratio of actual to predicted reduction
//			rho = getRho(lambda,p,Grad,x,Hes,x0,G,AT,numSS,numTotal);
//
//			// Adjust delta and make step based on rho
//			if (rho < 0.25) {
//				delta = norm(p,numSS)/4.0;
//			}
//			else if (rho > 0.75 && fabs(norm(p,numSS) - delta) < NUM_PRECISION) {
//				delta = min2(2.0*delta,deltaBar);
//			}
//			if (rho > eta) {
//				for (i = 0; i < numSS; i++) {
//					lambda[i] += p[i];
//				}
//				nNoStep = 0;
//			}
//			else {
//				nNoStep++;
//			}
//
//			// Calculate the mole fractions of the complexes based on lambda
//			if (getx(x,lambda,G,AT,numSS,numTotal) == 0) {// Should be fine;checked prev.
//				if (quiet == 0) {
//					printf("Overflow error in calcution of mole fractions.\n\n");
//					printf("Exiting....\n");
//				}
//				exit(ERR_OVERFLOW);
//			}
//
//			// Calculate the gradient
//			getGrad(Grad,x0,x,A,numSS,numTotal);
//
//			// Advance the iterations count
//			iters++;
//		}
//
//		// Advance the number of perturbations we've tried
//		nTrial++;
//
//	}
//
//	// Compute the free energy
//	FreeEnergy = 0;
//	// First the reference free energy
//	for (i = 0; i < numSS; i++) {
//		FreeEnergy += x0[i]*(1.0 - log(x0[i]));
//	}
//	// Now the free energy
//	for (j = 0; j < numTotal; j++) {
//		if (x[j] > 0) {
//			FreeEnergy += x[j]*(log(x[j]) + G[j] - 1.0);
//		}
//	}
//	// Convert to kcal/liter of solution
//	FreeEnergy *= kT*MolesWaterPerLiter;
//
//
//
//	// Return convergence
//	if (nTrial == MaxTrial) {
//		return 0;
//	}
//	else {
//		return 1;
//	}

}

void Concentration::update_complex_free_energy_vector(double *free_energy){
/*
 *
 *
 */

	for(int k = 0; k < TotalNumberOfComplex; k++){
		rna_complex_free_energy[k] = free_energy[complex_index_vector[k]];
	}
}

/***************************Getters******************************************/


void Concentration::get_complex_free_energy_vector(double *free_energy){
	for(int k = 0; k < TotalNumberOfComplex; k++){
		free_energy[complex_index_vector[k]] = rna_complex_free_energy[k];
	}
}

double Concentration::get_total_free_energy(){
	return TotalFreeEnergy;
}

double Concentration::get_concentration(int pos){
	if(pos < 0){
		return rna_complex_concentration[main_complex_position];
	}

	return rna_complex_concentration[pos];
}

bool Concentration::get_claculation_status(){
	return calculation_flag;
}



/***************************Setters******************************************/

void Concentration::reset_init_flag(){
	is_initialized = false;
}

void Concentration::set_iteration_values(int max_iteration,int max_trial){
	MaxIterationValue = max_iteration;
	MaxTrialValue = max_trial;
}

int Concentration::get_max_iteration(){
	return MaxIterationValue;
}

int Concentration::get_trial_iteration(){
	return MaxTrialValue;
}

int Concentration::get_total_complex(){
	return TotalNumberOfComplex;
}

int	Concentration::get_complex_index(int idx){
	return complex_index_vector[idx];
}

int	Concentration::get_total_rna(){
	return NumberOfSingleStrand;
}

