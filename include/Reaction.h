/*
 * Reaction.h
 *
 *  Created on: Aug 30, 2010
 *      Author: kamal
 */

#ifndef REACTION_H_
#define REACTION_H_

#include "RNAComplex.h"
#include "permutation_combination.h"
#include <streambuf>
#include "ReadWrite.h"
#include "nupack_conc_object.h"
#include "nupack_object.h"


#define		CONFIGURATION_HEADER		"RNA_Configuration"
#define		TRUTH_TABLE_HEADER			"Truth_Table"
#define		SEQUENCE_LABEL				"Sequence"
#define		SEQUENCE_CONSTARINT_LABEL	"Sequence_Constraint"
#define		STRUCTURE_LABEL				"Structure"
#define		TARGET_STRUCTURE_LABEL		"Target_Structure"
#define		ROW_LABEL					"Row"
#define		COLUMN_LABEL				"Column"
#define		MAX_ITERATION_LABEL			"Maximum_Iteration"


#define		MONTE_CARLO_MODE			false
#define		SIMULATED_ANNEALING_MODE	true

#define		BASIC_REPEAT_CHECK			1
#define		EXTENDED_SEQUENCE_CHECK		0

#define		WRITE_TO_FILE				1
#define		WRITE_TO_STDOUT				2
#define		WRITE_TO_FILE_AND_STD		(WRITE_TO_FILE | WRITE_TO_STDOUT) // 3 :-)
#define		CALCULATE_CONC_RANGE		4
#define		CALCULATE_SUBSEQUENCE_FREE_PROBABILITY	8

#define 	MAX_LINES_PER_RNA_DATA		5

#define		RBS_FREE_THRESHOLD			0.99//RBS_FREE_THRESHOLD_VALUE	// upper limit from which probabilistic calculation starts
#define		RBS_BOUND_THRESHOLD			0.8//RBS_BOUND_THRESHOLD_VALUE

#define		DENSITY_OF_WATER			55.14 // mol/L


#define		ADAPTIVE_SCORE_WEIGHTING	1
#define		FIXED_SCORE_WEIGHTING		0


extern int truth_table_OR_gate[14];
extern int truth_table_AND_gate[14];
extern int truth_table_XOR_gate[14];

extern int truth_table_NOR_gate[14];
extern int truth_table_NAND_gate[14];
extern int truth_table_XNOR_gate[14];

extern int truth_table_YES_gate[6];
extern int truth_table_NOT_gate[6];


typedef struct{
	int		row;
	int		col;
	bool	table[1 << (MAX_TRUTH_TABLE_COLUMNS - 1)][MAX_TRUTH_TABLE_COLUMNS];
	int		interaction_table[1 << (MAX_TRUTH_TABLE_COLUMNS - 1)][MAX_TRUTH_TABLE_COLUMNS];
}TruthTable;





class Reaction{
private:
	int				id;
	string			gate_name;

	Concentration	*concentration_obj[MAX_RNA_COMPLEX];

	RNAComplex		dimer;
	RNAComplex		trimer;
	RNAComplex		*main_rna_complex;
	RNAComplex		*rna_complex[MAX_RNA_COMPLEX];

// *** variables used for RNA complex concentration calculation ***

	int				total_rna;			// Total number of RNA single strands are shared by all the complexes
	int				total_truth_table_rna_complex;	// Total number of complexes.
	int				extra_rna_complex;		// extra RNA complexes generated for concentration calculation
//	int				total_complex;

	int				max_iteration;		// Maximum number of Monte-Carlo steps
	double			weight;

	int				mutation_word[MAX_CHAR];
	int				mutation_word_len;
	int				mutation_position;

	double			total_free_energy;
	double			energy;
	double 			concentration;

	double			sscore;
	double			score_kin;
	double			score_dist;
	double			score_weight;

	TruthTable		truth_table;

	int				truth_table_row;
	int				truth_table_column;

	ifstream		input;					// input and output files
	ofstream		output;
	ofstream		energy_output;
	char			input_file_name[300];
	char			output_file_name[300];

	int				output_mode;	// 1 = write to file only
									// 2 = write to std only
									// 3 = 1 | 2 = write to both file and display

	time_t			start_time;
	time_t			end_time;
	double			total_time_taken;
	int 			mc_steps_count;
	int				total_mc_count;
	bool			calculation_mode;
	bool			score_weighting_mode;
	int				calculation_evolution_mode; // 1 for structure and 2 for probability

	long			total_cpu_count;	// number of processors in the host computer


public:
	Reaction();
	Reaction(RNAComplex rnac, RNAComplex target);
	Reaction(RNAComplex &rnac);
	~Reaction();

	void		add_to_dimer(string seq,string str);
	void		add_to_trimer(string seq,string str);
	void		add_to_target(string str);
	void		add_rna_complex(RNAComplex &rnac);
	void 		add_extra_rna_complex(RNAComplex &rnac);


	double		calc_score_kin();	// not used
	double		calc_score_struct(); // not used

	int			get_id();
	double		get_energy();
	double		get_concentration();
	double 		get_concentration(int complex_idx);
	double 		get_concentration_fraction(int complex_idx);

	double		get_temperature();
	double		get_total_free_energy();
	double		get_kon_rate();
	double		get_keq_rate();
	void		get_rna_and_complex_vectors(double *rna_energy_vector, double *complex_energy_vector);
	void 		get_combinatorial_complexes(int complex1, int *idx_data, int &idx_len,int &main_complex_pos);
	int			get_mc_steps_count();


	void 		set_max_iteration(int max_itr = 1000);
	void		set_output_file_name(string &file_name);
	void		read_input_file(string &file_name);

	bool		evolve(int word_len);

	void		create_reaction_object_from_file(char *infile,char *outfile);
	void 		create_reaction_object_from_file_v2(char *infile,char *outfile);

	void 		create_rna_complex_from_truth_table(RNAComplex *full_rna_complex,int *truth_table, int mrna_position = 0,bool place_mrna_first = true);
	void		create_mass_matrix();

	void 		calculate_shared_rna_count();
	void		rna_complex_random_select(int rnac_num, int *rnac_indeces);
	bool 		check_complexes_with_shared_rna(int complex1, int complex2);
	bool		check_combinatorial_complexes(int complex1, int complex2);
	void		create_vector_from_complex(int complex1,int *complex_vec);

	double		calculate_total_rnagate_score(double &prob_score, double &conc_score); // calculates the total free energy of all the RNAComplexes involved
//	double		calculate_free_energy();
	double		calculate_energy();
	void		calculate_all_complex_concentration(bool force_all_calculation = false);
	double 		calculate_concentration(int complex_id);

	void 		calculate_total_rbs_distance(double &rbs_free_dist,double &rbs_bound_distance);
	double 		calculate_complex_mean_pvalue();

	void 		calculate_combinatorial_concentration(double cmin,double cmax,string &filename);

	static int	sequence_distance(string &seq1, string &seq2);

	double 		sequence_alignment_score(string &seq1,string &seq2);
	int 		calculate_2D_sequence_map(string *seq, int total_seq);
	static void		calculate_seqence_entropy(vector<string> &seq, double *entropy_out,double frequency_aa[][4]);





	void		copy_rna_complex(RNAComplex &rnac);
	void		update_rna_complex(const RNAComplex &rnac);
	int**		generate_A_matrix(int max_rna, int max_complex, int max_rna_in_complex);

	void 		update();

	void		read_input(char *infile);

	void		set_input_file(const char *fname);
	void		set_output_file(const char *fname);
	void		set_output_mode(int mode);
	void		set_weight(double new_weight);
	void		set_calculation_mode(bool cmode);
	void		set_sequence_check_mode(int smode);
	void		set_calculation_evolution_mode(int cmode);
	void		set_score_weighting_mode(bool mode);


	double		get_sscore();

	double		get_weight();
	bool		get_calculation_mode();
	int			get_calculation_evolution_mode();
	bool		get_score_weighting_mode();



	void		write_output();

	void		write_results();
	void		write_results(ofstream *fptr);

	void		run();
};


#endif /* REACTION_H_ */
