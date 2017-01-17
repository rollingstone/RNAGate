/*
 * nupack_extension.h
 *
 *  Created on: Nov 19, 2010
 *      Author: kamal
 */

#ifndef NUPACK_EXTENSION_H_
#define NUPACK_EXTENSION_H_

//#ifdef __cplusplus
//extern "C"{
//#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

//#include "ConcentrationsHeaderFile.h" // Header file for concentrations


//#include "pfuncUtilsHeader.h"
//#include "DNAExternals.h"


//void 	init_conc_data(int numSS, int numTotal);
//bool 	calculate_concentration(double *x, int **A, double *G, double *x0);
//void 	reset_init_flag();
//void 	set_iteration_values(int max_iteration,int max_trial);
//int  	get_max_iteration();
//int  	get_trial_iteration();
//double 	get_total_free_energy();

//#ifdef __cplusplus
//}
//#endif

using namespace std;


#define RNA_MAX			30
#define COMPLEX_MAX		30

class Concentration{

private:
//*************** variables required for calculate_concentration ******************

	double 		TolarenceValue;
	double 		DeltaBarValue;
	double 		TemperatureKTValue;
	int			MaxTrialValue;
	double 		PerturbScaleValue;
	int			DoNotShowErrorMessage;
	int 		MaxNoStepValue;
	int			MaxIterationValue;
	double 		MolesWaterPerLiterValue;
	double 		ETA_Value;
	double		TotalFreeEnergy;

	double		NumberOfSingleStrand;
	int 		TotalNumberOfComplex;


	bool		is_initialized;

	double 		*AbsTol; // The absolute tolerance on all values of gradient
	double 		rho; // Ratio of actual to predicted reduction in trust region method
	double 		delta; // Radius of trust region
	double 		*Grad; // The gradient of -g(lambda)
	double 		*lambda; // Lagrange multipliers (dual variables),x[j] = Q[j]*exp(lambda[j])
	// for j \in \Psi^0
	double 		*p; // The step we take toward minimization
	double 		**Hes; // The Hessian
	long 		*idum; // Input for random number generator
	int 		**AT; // Transpose of A
	int 		RunStats[6]; // Statistics on results from getSearchDir (see comments below)
// ***********************

	double		*rna_concentration_vector;
	double		*rna_complex_concentration;
	double		*rna_complex_free_energy;
	int			main_complex_position;
	int			*complex_index_vector;
	int			**A_matrix_main;
	int			complex_ptr[COMPLEX_MAX];
	bool		calculation_flag;

public:
	Concentration();
	Concentration(const Concentration &conc);
	Concentration(int numSS, int numTotal);
	~Concentration();

	bool 	calculate_concentration();
	bool 	calculate_concentration(double *complex_free_energy);
	bool 	calculate_concentration(double *rna_conc,double *complex_conc);

//	bool 	calculate_concentration(double *complex_conc, double *complex_free_energy);
	void 	calculate_combinatorial_concentration(double rna_conc_min,double rna_conc_max,  double min_rna_conc, ofstream &fconc);
	bool 	calculate_concentration(double *x1, int **A1, double *G1, double *x0);
//	void	calculate_combinatorial_concentration(double rna_cocnc_range[2],double *complex_concentration);


	void 	init_object(int numSS, int numTotal);
	void 	init_parameters(int total_rna, int total_complex,double *rna_conc, int **A1, int *complex_index, int main_cpos);
	void	free_object();
	void	print_A_matrix();

	double 	get_total_free_energy();
	double	get_concentration(int pos = -1);
	int		get_total_complex();
	int		get_total_rna();
	int		get_complex_index(int idx);
	bool 	get_claculation_status();
	void	update_complex_free_energy_vector(double *free_energy);
	void 	get_complex_free_energy_vector(double *free_energy);


	void 	reset_init_flag();
	void 	print_text(char *t);
	void 	set_iteration_values(int max_iteration,int max_trial);
	int 	get_max_iteration();
	int 	get_trial_iteration();
};


#endif /* NUPACK_EXTENSION_H_ */
