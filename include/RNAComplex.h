/*
 * RNAComplex.h
 *
 *  Created on: Aug 29, 2010
 *      Author: kamal
 */

#ifndef RNACOMPLEX_H_
#define RNACOMPLEX_H_

#include "RNAObject.h"
#include "random.h"
//#include <cstring>

using namespace std;

#define		RNA_COMPLEX_SPACE   '+'
//#define		RNA_COMPLEX_SPACE   ' '


#define		OFF_STATE	0
#define		ON_STATE	1


class RNAComplex:public RNAObject{
private:
	RNAObject*			rna[MAX_RNA_NUMBER];
	int					rna_self_bound_sites[MAX_RNA_NUMBER];
	int					rna_total_bound_sites[MAX_RNA_NUMBER];
	int 				rna_rna_bound_sites[MAX_RNA_NUMBER][MAX_RNA_NUMBER];

	int 				rna_position_base[MAX_RNA_NUMBER];
	int					rna_position_end[MAX_RNA_NUMBER];
	int					rna_toehold_length[MAX_RNA_NUMBER];
	int					total_toehold_length;
	int					mean_toehold_length;

	double				max_rna_concentration;
	double				min_rna_concentration;


	string				toehold_structure;

	int					total_bound_sites;
	int					self_bound_sites;
	double				self_bound_fraction;

	double				kon_rate; // reaction rate Kon
	double				keq_rate; // equilibrium reaction rate

	int					rbs_free_target_site_count;
	int					rbs_bound_target_site_count;

//	int					max_mutation_locations;

	int					total_rna_number;
	int					**rna_combination;
	int					max_combination;
	bool				is_target_structure_updated;

	bool				is_it_extra_rna_complex;

	int 				mrna_position;

	int					mass_vector[MAX_RNA_NUMBER];
	int					mass_vector_length;

	bool				target_state;
	double				rbs_translation_rate;
	double				rbs_g_total;


public:

	RNAComplex();
	RNAComplex(RNAObject *rna1);
	RNAComplex(RNAObject *rna1,RNAObject *rna2);
	RNAComplex(RNAObject *rna1,RNAObject *rna2, RNAObject *rna3);
	RNAComplex(RNAObject **rna,int rna_number);

	virtual ~RNAComplex();

	RNAComplex& operator=(RNAComplex &rnac);

	void		create_rna_complex_from_truth_table(int *truth_table);

	void		add_rna(RNAObject *rna);
	void 		add_rna(RNAObject &new_rna);
//	void		add_rna(string sequence, string structure);

	void		update(int calculation_mode = EVOLVE_FOR_STRUCTURE);
	void		update_rna_objects(const string &seq, const string &str);
	void		update_rna_objects(int calculation_evolve_mode = EVOLVE_FOR_STRUCTURE);
	void		update_complex_parameters();
	void		update_target_parameters();
	void		update_rnacomplex_sequence_and_structure();
	void		update_rna_sequence_to_complex();

	bool		mutate(int wlen, int evolution_mode = EVOLVE_FOR_STRUCTURE);

	void		calculate_structure_distance(string &str,string &target,int *total_distance,int *rbs_distance);
	void		calculate_structure_distance(int *total_distance,int *rbs_distance);
	void 		calculate_binding_fractions();
	int 		calculate_mutation_rna(int *priority_position_vector,int *vector_len);
	int			calculate_toehold_length(bool force_any_rna_calculation = false);
	double		calculate_rbs_translation_rate();

	int			get_total_rna_number();
	RNAObject&	get_rna(int rna_id);
	double		get_energy_difference();
	string&		get_complex_target_structure();
	string& 	get_toehold_structure();
	double		get_score_change();
	int			get_mrna_position();
	double 		get_self_bound_fraction(int rna_id = -1); 					 // default rna_id = -1 represents mRNA id in the complex
	double		get_rna_rna_binding_fraction(int rna1_id, int rna2_id = -1); // default rna_id = -1 represents mRNA id in the complex
	int			get_total_toehold_length();
	int			get_rna_toehold_length(int rna_id);
	int			get_mean_toehold_length();
	int 		get_rna_starting_position_in_complex(int rna_id);
	int			get_mutation_locaion_number();
	int			get_rna_id_from_position(int pos);
	double		get_kon_rate();
	double		get_keq_rate();
	bool 		get_extra_rna_complex_flag();
	void		get_rna_free_energy(double *energy_vector);
	double 		get_min_rna_concentration();
	bool 		get_target_state();
	double		get_rbs_translation_rate();


	void		generate_mass_vector(int *mvector);

	void		set_mrna_position(int pos);
	void		save_parameters();
	void		restore_parameters();
	void		save_current_score();
	void		set_rbs_target(bool flag);
	void		set_mutation_locaion_number(int ml_number);
	void		set_extra_rna_complex_flag(bool cflag);
	void		set_target_state(bool flag);

	void		reset_component_rna_id();
	void		reset_mutation_flag();

	void		write_results();
};

#endif /* RNACOMPLEX_H_ */
