/*
 * RNAObject.h
 *
 *  Created on: Aug 25, 2010
 *      Author: kamal
 */

#ifndef RNAOBJECT_H_
#define RNAOBJECT_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <ctime>
#include <cmath>
#include <streambuf>
#include <ctime>
#include <omp.h>


// *************************** LocARNA **************************

#include <iostream>
#include <fstream>
#include <vector>

#include <memory> // for auto_ptr


// *************************** End of LocARNA **************************

//#include  <stdlib.h>

#include "global_variables.h"



extern "C"{
#include "fold.h"
#include "inverse.h"
#include "utils.h"
}

#include "global_variables.h"
#include "random.h"

//#include "nupack_object.h"
//#include "nupack_mfe_pfunc.h"

extern "C"{
#include "pfuncUtilsConstants.h"
#include "pfuncUtilsHeader.h"
//#include "pfuncUtilsHeader.h"
#include "DNAExternals.h"
#include "ConcentrationsHeaderFile.h"
}

//#include <gc.h>
//#include <leak_detector.h>


using namespace std;


#define	ABSOLUTE_ZERO_IN_KELVIN			273.15	// K
#define	BOLTZMANN_Kb					3.297626e-27 // kcal/K
#define	GAS_CONSTANT_R					1.985877e-3 // Kcal/K/mol

#define	KT37							((37 + ABSOLUTE_ZERO_IN_KELVIN)*BOLTZMANN_Kb)

#define RANDOM 			( (double) (rand()/(RAND_MAX+1.0)) ) //random number 0-1
#define MAXCHAR 		5000 //maximum number of chars in a string (do not change)
#define	MAX_CHAR		5000
#define KT 				0.616 //energy constant in Kcal/mol at 37 C
#define TOEHOLD 		6 //toehold of saturation for scoring
#define DELTAG 			15 //DG of saturation for scoring
#define DG_BP 			1.28 //effective DG per bp
#define DG_EFF			1.28 // kcal/M
#define KON 			0.27 //minimal k_on
#define DIST 			3 //maximal distance in structure
#define REPEATS 		3 //maximal number of repeats


#define Kb				0.00198721	// Kcal/mol/K


#define MAX_RNA_NUMBER				10
#define MAX_RNA_DIMER_NUMBER		(MAX_RNA_NUMBER*1)	// needs adjusting Permutation(MAX_RNA_NUMBER,2)
#define MAX_RNA_TRIMER_NUMBER		(MAX_RNA_NUMBER*1)  // needs adjusting Permutation(MAX_RNA_NUMBER,2)
#define REPEATS						3 // Consecutive nucleotides of sae kind

#define	OPTIMAL_TOEHOLD				6


#define	MAX_RNA_COMPLEX				MAX_RNA_NUMBER*3	// roughly
#define MAX_TRUTH_TABLE_COLUMNS		5

//********** The wildcards and symbols for target RNA structures ***************
#define	RBS_SITE			'R'
#define RBS_T1_FREE			'Q'
#define RBS_T0_BOUND		'P'
#define	RBS_BOUND			'B'
#define RBS_FREE			'F'
#define	SITE_BOUND			'|'
#define	SITE_FREE			'.'
#define SITE_LEFT_BOUND		'('
#define SITE_RIGHT_BOUND	')'
#define	ANYTHING			'@'
#define ANY_NUCLEOTIDE		'N'



#define	LINE_BREAK			'\n'
#define	LB					'\n'



#define EXP_FUNC(x)			exp((double)(x))

////********************** Input file tags ***********************
//
//#define		CONFIGURATION_HEADER		"RNA_Configuration"
//#define		TRUTH_TABLE_HEADER			"Truth_Table"
//#define		SEQUENCE_LABEL				"Sequence"
//#define		SEQUENCE_CONSTARINT_LABEL	"Sequence_Constraint"
//#define		STRUCTURE_LABEL				"Structure"
//#define		TARGET_STRUCTURE_LABEL		"Target_Structure"
//#define		ROW_LABEL					"Row"
//#define		COLUMN_LABEL				"Column"
//#define		MAX_ITERATION_LABEL			"Maximum_Iteration"
//
////*************************************************************


#define		FOUR_REPEAT_MODE		true
#define		SIX_EXTRA_CONFIG_MODE	false

#define		BASIC_REPEAT_CHECK			1
#define		BASIC_SEQUENCE_CHECK		1
#define		EXTENDED_SEQUENCE_CHECK		0


#define		STEM_PRESERVE_ON			1
#define		STEM_PRESERVE_OFF			0


/*
 *
 * RNAObject class contains the sequence, sequence constraint, structure and target structure information of the
 * RNAs.
 *
 * all the variables are names are self explanatory.
 *
 *
 */

#define		MAX_NUPACK_SEQLEN		1000
#define		EVOLVE_FOR_PROBABILITY  0xCAFE
#define		EVOLVE_FOR_STRUCTURE	0xBEE

#define		WOBBLE_PROBABILITY		0.4


class RNAObject{
private:
	int			id;
	string		rna_name;

	string		sequence;
	string		structure;

//	NUPACKObject	*nu_obj; // not used
	dnaStructures	mfe_structure;
	bool			nu_obj_initialized;

	int			num_sequence[MAX_NUPACK_SEQLEN];
	int			num_structure[MAX_NUPACK_SEQLEN];

	string		sequence_constraint;
	string		structure_constraint;
	string		target_structure;

	double		str_weight[MAX_CHAR];
	int			str_coords[MAX_CHAR];


	int			length;
	int			base_pair_length; // number of base pairs present in the structures
	double		energy;
	double		energy_per_base; // average binding energy per base
	double		temperature;
	double		concentration;
	int			total_rna_starnds;

	double		probability_value;
	DBL_TYPE	pfunc_z_value;

	int			total_distance;
	int			rbs_distance;		// distance between the target and current structure in the RBS region

	// save_ and old_ variables are used for storing purpose only
	string		save_sequence;
	string		save_structure;
	double		save_energy;
	double		save_energy_per_base;
	double		save_concentration;
	int			stack_height;
	int			save_total_distance;
	int			save_rbs_distance;
	double		save_probability_value;
	DBL_TYPE	save_pfunc_z_value;

	double		old_energy;
	double		old_score;
	double		old_concentration;

	int			G_count;
	int			C_count;
	int			A_count;
	int			U_count;
	int			total_bases;

	double		GC_fraction;


	bool		is_mutated; // true is the mutation occurred in the RNA
	bool		save_is_mutated;

	bool		is_rbs_site_unrestricted;

	int			max_mutation_locations;

	int			base_repeat_limit;
	int			distance_limit;

	int			sequence_check_mode;

	int 		do_symetry_correction;
	int 		do_calculate_pairs;
	int			subopt_energy_gap;
	int			complexity_num;
	int			perm_symmetry;
	int			temp_val;

	int 		evolve_mode; // EVOLVE_FOR_STRUCTURE or EVOLVE_FOR_PROBABILITY
	int			rna_structure_mutation_mode; // with or without preserving stems

//	static int	rna_object_count;


public:
	RNAObject();
	RNAObject(string &seqin, string &seq_constraint_in, string &strin, string &target, int id_number);
	~RNAObject();

	RNAObject& 	operator=(RNAObject &rna);
	RNAObject&	operator<<(RNAObject &rna);
	RNAObject&	operator>>(RNAObject &rna);


	void 		init_nupack_object(string &seq,bool do_reinitialization = false);
	void 		init_nupack_object(int *num_seq,int seq_len);

	bool		check_sequence();
	bool		check_structure();

	bool 		check_structure(const string &str);
	string 		reverse_complementary (string seq);
	string 		to_upper(string sequence);
	void 		string_upper(string &seq);

	void		generate_str_weight();

	int			calculate_distance();
	void		calculate_structure_weight();

	bool		directed_mutation(int word_len);
	bool		mutation(const char *word, int wlen);
	bool 		mutation(const string &word);
	bool		mutation(const string &word,int *position);

	double 		fold_mfe(string &seq, string &str,int calc_evolve_mode,DBL_TYPE &pf_z_value,int calc_one_mfe = 1);

	DBL_TYPE	calculate_pfunc();
	DBL_TYPE 	calculate_pfunc(string &seq);
	double		calculate_probability(DBL_TYPE z_value = -1);

	bool 		single_mutation(char w, int p, int check_more = BASIC_REPEAT_CHECK,int mutation_mode = STEM_PRESERVE_ON);

	bool		single_mutation_unconstrained(char w, int pos);
	void		generate_structure(const string &seq);
	void		generate_sequence(const string &str);
	int 		calculate_target_distance(const string &source,const string &target, int &rbs_distance);

	void 		update_sequence(const string &seq);
	void 		update_structure(const string &str);

	void		generate_random_word(char *word,int word_len);

	int 		compute_toehold(const string &str1, const string &str2, const string &str12);
	int 		distance_structures(const string &str, const string &target);

	int 		calculate_structure_mismatch(const char *str1, const char *str2,int str_len);

	double 		calculate_energy_per_base();
	double 		calculate_energy_per_base(char *str,int len, int &bp_length,double en);

	void		calculate_base_fractions();

	RNAObject*	duplicate();

	void 		update(int calculation_mode = EVOLVE_FOR_STRUCTURE);
	void		update_rna_sequence_to_complex();

	void		save_parameters();
	void		restore_parameters();

	bool		check_repeats(const string &seq);
	bool		check_repeats(const char *seq, int len = -1);

	/***** Define some static member functions **************/
	static bool	check_local_repeats(const string &seq,int pos = -1);
	static bool	check_no_go_sequences(const string &seq, int pos = -1);
	static int	find_complementary_site(int pos, const string &strin);


	static bool	check_valid_sequences(const string &seq, int pos = -1, int check_mode = BASIC_REPEAT_CHECK);
	static void	generate_random_sequence(char *seq, int len);
	static void	generate_unique_random_id(int max_val, int *rnac_indeces);
    static double getcputime();
	static char calc_non_comp_base(char nt);
	static char calc_comp_base(char nt);
	static void	make_rna_structure(const dnaStructures &ds,const char *seq,char *str,int index = 0);
	static char **make_multiple_rna_structures(const dnaStructures &ds,const char *seq,int &number_of_structures);
	static double subsequence_free_probability(string &seq, DBL_TYPE &pf,DBL_TYPE &pf_subopt_value, int &total_str,int &total_str_free, DBL_TYPE subopt_gap, int start_pos,int end_pos,int natype = 2, int dangle_type = 2);
	static void	calculate_subseq_free_probability(ifstream &infile,string &output_file_name);


    /********************** getters and setters ****************/
	int			get_id();
	string&		get_rna_name();
	string&		get_sequence();
	string&		get_structure();
	string&		get_sequence_constraint();
	string&		get_structure_constraint();
	string&		get_target_structure();
	int			get_sequence_length();
	int			get_base_pair_length();
	double		get_energy();
	double		get_energy_per_base();
	double		get_save_energy();
	double		get_concentration();
	double		get_temperature();
	RNAObject*	get_rna_pointer();
	int			get_total_rna_number();
	int			get_rbs_distance();
	int			get_total_distance();
	int	 		get_distance_limit();
	double		get_energy_difference();
	bool 		get_mutation_flag();
	int			get_max_mutation_locations();
	int			get_repeat_limit();
	double		get_GC_fraction();
	double		get_AU_fraction();
	int			get_sequence_check_mode();
	void		get_base_content_fraction(double &A_frac,double &U_frac,double &G_frac,double &C_frac);
//	NUPACKObject& get_nupack_object();
	double 		get_probability();
	DBL_TYPE 	get_pfunc_z_value();
	int			get_evolve_mode();
	void		get_total_free_distance(int &distance, int &free_site_count);
	void		get_total_bound_distance(int &distance, int &bound_site_count);
	int			get_subopt_energy_gap();
	int			get_rna_structure_mutation_mode();

	void		set_id(int id_value);
	void		set_rna_name(string &name);
	void		set_energy(double energy_value);
	void		set_energy_per_base(double energy_pb_value);
	void 		set_sequence(const string &seq);
	void 		set_structure(const string &str);
	void		set_sequence_constraint(const string &seq);
	void		set_structure_constraint(const string &str);
	void		set_target_sequence(const string &seq);
	void		set_target_structure(const string &str);
	void		set_base_pair_length(int bp_length);
	void		set_rbs_distance(int rbs_dist);
	void		set_total_distance(int total_dist);
	void 		set_distance_limit(int dist_limit);
	void 		set_concentration(double conc_value);
	void		set_old_energy(double old_energy);
	void 		set_mutation_flag(bool flag);
	void		set_repeat_limit(int rlimit);
	void		set_max_mutation_locations(int mutation_locaitons);
	void		set_sequence_check_mode(int cmode);
	void 		set_probability(double prob_val);
	void 		set_pfunc_z_value(DBL_TYPE zvalue);
	void		set_evolve_mode(int e_mode);
	void		set_subopt_energy_gap(int egap);
	void		set_rna_structure_mutation_mode(int mode);

};


#endif /* RNAOBJECT_H_ */
