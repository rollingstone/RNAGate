/*
 * RNAComplex.cpp
 *
 *  Created on: Aug 26, 2010
 *      Author: kamal
 */


#include "RNAComplex.h"

#include <cstring>

using namespace std;


string	superfolded_gfp_seq = "CGUAAAGGCGAAGAGCUGUUCACUGGUGUCGUCCCUAUU"; // (BBa_I746916)
string 	ribosome_16s_seq 	= "ACCUCCUUA";
string	ribosome_30s_seq 	= "CAU";

/***************** Constructors **************************/

RNAComplex::RNAComplex():RNAObject(){
	total_rna_number = 0;
	is_target_structure_updated = false;
	is_it_extra_rna_complex	= false;

	max_rna_concentration = -9999;
	min_rna_concentration =  9999;

	for(int i = 0; i < MAX_RNA_NUMBER; i++)
		mass_vector[i] = 0;
}

RNAComplex::RNAComplex(RNAObject *rna1){
	RNAComplex();

	add_rna(*rna1);
	update_target_parameters();
	update();
}


RNAComplex::RNAComplex(RNAObject *rna1,RNAObject *rna2){
	RNAComplex();

	add_rna(*rna1);
	add_rna(*rna2);
	update_target_parameters();
	update();
}

RNAComplex::RNAComplex(RNAObject *rna1,RNAObject *rna2, RNAObject *rna3){
	RNAComplex();

	add_rna(*rna1);
	add_rna(*rna2);
	add_rna(*rna3);
	update_target_parameters();
	update();
}


RNAComplex::RNAComplex(RNAObject **rna,int rna_number){
	RNAComplex();

	for(int i = 0; i < rna_number; i++)
		add_rna(*rna[i]);

	update_target_parameters();
//	init_nupack_object(true);
	update();
}

RNAComplex::~RNAComplex(){
// i kill myself :P
}

/******************** Method implementations **************************/

void RNAComplex::add_rna(RNAObject &new_rna){
	string	sp;

	mass_vector[new_rna.get_id()]++;

	if(total_rna_number == 0){
		rna_position_base[total_rna_number] = 0;
		rna_position_end[total_rna_number]  = new_rna.get_sequence_length() - 1;
	}
	else{
		rna_position_base[total_rna_number] = rna_position_base[total_rna_number-1] + rna[total_rna_number-1]->get_sequence_length() + 1;
		rna_position_end[total_rna_number]  = rna_position_end[total_rna_number-1] + new_rna.get_sequence_length() + 1;
	}

	rna[total_rna_number++] = &new_rna;


	for(int i = 0; i < total_rna_number; i++){
		if(rna[i]->get_concentration() <= min_rna_concentration)
			min_rna_concentration = rna[i]->get_concentration();
	}

	update_target_parameters();
	update_rnacomplex_sequence_and_structure();
}


void RNAComplex::add_rna(RNAObject *rna){
	add_rna(*rna);
}


int RNAComplex::calculate_mutation_rna(int *priority_position_vector,int *vector_len){
/*
 *
 *
 *
 *
 */
	int		pos_idx_count = 0;
	int		len = (int) get_sequence_length();
	int 	pos;
	int		comp_pos;
	int		count = 0;

	int 	rna_id = -1;

	for(int i = 0; i < len; i++){
		if(count >= rbs_free_target_site_count)
			break;

		if(get_target_structure()[i] == RBS_FREE){
			count++;
			if(get_structure()[i] != '.'){
				pos = i;
				comp_pos = find_complementary_site(pos,get_structure());

				if(rna_id < -1){
					rna_id = get_rna_id_from_position(comp_pos);
					if(rna_id == get_rna_id_from_position(pos)){
						pos_idx_count = 0;
						break;
					}
				}

				if(rna_id == get_rna_id_from_position(comp_pos)){
					priority_position_vector[pos_idx_count++] = comp_pos - get_rna_starting_position_in_complex(rna_id);
				}
			}
		}
	}

	int	k;

	if(pos_idx_count == 0)
		rna_id = (int) floor(random_num()*total_rna_number);
	else{
		for(k = pos_idx_count; k < len; k++){
			priority_position_vector[k] = floor(random_num()*len);
		}
		pos_idx_count = k;
	}

	*vector_len = pos_idx_count;

	return rna_id;
}


bool RNAComplex::mutate(int wlen, int evolution_mode){
/*
 * This function mutates the RNAComplex object with a word or multiple random words given parameter wlen
 * at number of randomly chosen location given by RNAComplex object member "max_mutation_location" which
 * has default value of 1.
 *
 * Features:
 *
 * 1. The function mutates only the rna molecules that are not already mutated somewhere else.
 * 2. It tries to preserve the original structure of the component rna molecules.
 * 3. It returns TRUE only if the state of the newly calculated rna complex approaches the target structure or remains in their current
 *    distance from the target structures. Here, distance is the number of structural mismatch between current and target structure.
 *
 */
	char	new_structure[MAX_CHAR];
	bool	mutation_flag;
	char	word[20];
	int		dist;
	int 	rbs_dist;
	int 	max_iteration = 100;
	double 	ev_temp;
	int 	monomer_distance;
	string	str_tmp;

	string 	seq;
	string 	str;
	double 	en;
	double	en_per_base;
	int		len;
	int		pos;
	int		mcount = 0;

	int 		rna_mutation_count;
	double		old_prob = 0;
	DBL_TYPE 	old_zvalue = 0;
	double 		new_prob;
	DBL_TYPE	new_zvalue;
	int 		selected_rna;

	int		rna_rand_id[MAX_RNA_NUMBER];

	generate_unique_random_id(total_rna_number,rna_rand_id);

	for(int idx = 0; idx < total_rna_number; idx++){
		selected_rna = rna_rand_id[idx];//floor(random_num() * total_rna_number);

		if(rna[selected_rna]->get_mutation_flag() == true){
			mcount++;
			continue;
		}

		seq = rna[selected_rna]->get_sequence();
		str = rna[selected_rna]->get_structure();
		en 	= rna[selected_rna]->get_energy();
		old_prob = rna[selected_rna]->get_probability();
		old_zvalue = rna[selected_rna]->get_pfunc_z_value();
		en_per_base = rna[selected_rna]->get_energy_per_base();
		len = (int) rna[selected_rna]->get_sequence().length();

		int	mutation_locations = rna[selected_rna]->get_max_mutation_locations();
		rna_mutation_count = 0;

		for(int j = 0; j < max_iteration; j++){
			mutation_flag = false;
			for(int mc = 0; mc < mutation_locations; mc++){ // mutate in locations given by "max_mutation_locations"
				pos = (int) floor((len-wlen) * random_num());
				generate_random_word((char *)word,wlen);

				for(int k = 0; k < wlen; k++){
					bool sm_ret = rna[selected_rna]->single_mutation(word[k],pos,get_sequence_check_mode());
					mutation_flag |= sm_ret;
				}
			}

			if(mutation_flag == false)
				continue;

			string 		tmp_str;

			ev_temp = rna[selected_rna]->fold_mfe(rna[selected_rna]->get_sequence(),tmp_str,evolution_mode,new_zvalue);

			strcpy(new_structure,tmp_str.c_str());

			if(evolution_mode == EVOLVE_FOR_PROBABILITY)
				new_prob = EXP_FUNC(-ev_temp/(Kb*(37+273)))/new_zvalue;
			else
				new_prob = old_prob = 0;

			int		bp_len;
			double 	ev_bp_new  = calculate_energy_per_base(new_structure,len,bp_len,ev_temp);

			monomer_distance = (int) bp_distance(rna[selected_rna]->get_structure_constraint().c_str(), (const char *)new_structure);

			if(monomer_distance <=  rna[selected_rna]->get_distance_limit()){
				if(true || !RNA_PROB_CALC_ENFORCE_FLAG || (new_prob >= 0.5 || new_prob >= old_prob))
				{

					string		cx_seq = get_sequence();
					string		cx_str = get_structure();
					double  	cx_energy = get_energy();
					double		cx_energy_per_base = get_energy_per_base();
					int			rd = get_rbs_distance();
					int 		td = get_total_distance();
					DBL_TYPE 	cx_zvalue = get_pfunc_z_value();
					double   	cx_prob_val = get_probability();

					update(evolution_mode);
					calculate_structure_distance(&dist,&rbs_dist);

					if(rbs_dist <= get_rbs_distance()){
						set_rbs_distance(rbs_dist);
						set_total_distance(dist);
						set_mutation_flag(true);

						rna[selected_rna]->set_structure(new_structure);
						rna[selected_rna]->set_energy(ev_temp);
						rna[selected_rna]->set_energy_per_base(ev_bp_new);
						rna[selected_rna]->set_mutation_flag(true);

						rna[selected_rna]->set_probability(new_prob);
						rna[selected_rna]->set_pfunc_z_value(new_zvalue);

						if(total_rna_number == 1){
							set_probability(new_prob);
							set_pfunc_z_value(new_zvalue);
						}
						return true;
					}

					set_sequence(cx_seq);
					set_structure(cx_str);
					set_energy(cx_energy);
					set_energy_per_base(cx_energy_per_base);
					set_rbs_distance(rd);
					set_total_distance(td);
					set_pfunc_z_value(cx_zvalue);
					set_probability(cx_prob_val);
					set_mutation_flag(false);
				}
			}

			rna[selected_rna]->set_sequence(seq);
			rna[selected_rna]->set_structure(str);
			rna[selected_rna]->set_energy(en);
			rna[selected_rna]->set_energy_per_base(en_per_base);
			rna[selected_rna]->set_mutation_flag(false);
			rna[selected_rna]->set_pfunc_z_value(old_zvalue);
			rna[selected_rna]->set_probability(old_prob);
		}
	} /* end of for(selected_rna = 0....*/

	if(mcount > 0)
		return true;

	return false;
}



void RNAComplex::update_target_parameters(){


	for(int i = 0; i < total_rna_number; i++){
		if(i == 0){
			set_sequence_constraint(rna[i]->get_sequence_constraint());
			set_structure_constraint(rna[i]->get_structure_constraint());
			set_target_structure(rna[i]->get_target_structure());
		}
		else{
			set_sequence_constraint(get_sequence_constraint() + RNA_COMPLEX_SPACE + rna[i]->get_sequence_constraint());
			set_structure_constraint(get_structure_constraint() + RNA_COMPLEX_SPACE + rna[i]->get_structure_constraint());
			set_target_structure(get_target_structure() + RNA_COMPLEX_SPACE + rna[i]->get_target_structure());
		}
	}
}

void RNAComplex::update_rnacomplex_sequence_and_structure(){

	for(int i = 0; i < total_rna_number; i++){
		if(i == 0){
			set_sequence(rna[i]->get_sequence());
			set_structure(rna[i]->get_structure());
		}
		else{
			set_sequence(get_sequence() + RNA_COMPLEX_SPACE  + rna[i]->get_sequence());
			set_structure(get_structure() + RNA_COMPLEX_SPACE + rna[i]->get_structure());
		}
	}
}

void RNAComplex::update_rna_objects(int calculation_evolve_mode){
	int		len = (int) get_sequence_length();
	int		rna_index;
	int		i0,i;
	int		chunk = 8;


	if(total_rna_number == 1){
		rna[0]->set_sequence(get_sequence());
		rna[0]->update(calculation_evolve_mode);
		return;
	}
//#pragma omp parallel private(i,i0,rna_index)
	{
		i0 = 0;
		rna_index = 0;

//#pragma omp for schedule(static,chunk) nowait
		for(i = 0; i < len; i++){

			//		if(rna_index > total_rna_number){
			//			cerr << "ERROR: in update_rna_objects\nrna_count exceeded total_rna_number.\nExiting program" << endl;
			//			exit(-1);
			//		}

			cerr << "in update_rna_objects....::::" << omp_in_parallel() << endl;
			if(rna_index == total_rna_number - 1){
				rna[rna_index]->set_sequence(get_sequence().substr(i0, (len - i0)));
				rna[rna_index]->update(calculation_evolve_mode);
				i = 9999; // to go out of the loop
				continue;
//				return;
			}

			if(get_sequence()[i] == ' ' || get_sequence()[i] == '+'){
				rna[rna_index]->set_sequence(get_sequence().substr(i0, (i-i0)));
				rna[rna_index]->update(calculation_evolve_mode);
				i0 = i + 1;
				rna_index++;
			}
		}
	}
}


void RNAComplex::update(int calculation_mode){
/*
 * Calculates structure and free energy of the Complexes using MultiRNAFold package.
 *
 * calculation_mode = EVOLVE_FOR_STRUCTURE, energy and structure is calculated
 * calculation_mode = EVOLVE_FOR_PROBABILITY, energy, structure, partition function and probability is calculated
 *
 */

	for(int i = 0; i < total_rna_number; i++){
		if(i == 0){
			set_sequence(rna[i]->get_sequence());
			set_structure(rna[i]->get_structure());
		}
		else{
			set_sequence(get_sequence() + RNA_COMPLEX_SPACE  + rna[i]->get_sequence());
		}
	}

	if(total_rna_number == 1){
		set_energy(rna[0]->get_energy());
		set_probability(rna[0]->get_probability());
		set_pfunc_z_value(rna[0]->get_pfunc_z_value());
	}else{
		if(calculation_mode == EVOLVE_FOR_PROBABILITY){
			RNAObject::update(EVOLVE_FOR_PROBABILITY);
		}
		else if(is_it_extra_rna_complex == false || calculation_mode != EVOLVE_FOR_PROBABILITY){
				RNAObject::update(EVOLVE_FOR_STRUCTURE);
		}
	}

	calculate_toehold_length();
	calculate_energy_per_base();

	kon_rate = KON * exp(DG_EFF* get_total_toehold_length()/KT);
	keq_rate = exp(-get_energy_difference()/KT);
}


double RNAComplex::get_energy_difference(){
/*
 * Calculates reaction free energy.
 *
 */
	double en = get_energy();

	if(total_rna_number <= 1)
		return 0;

	for(int i=0; i < total_rna_number; i++){
		en -= rna[i]->get_energy();
	}

	return en;
}

void RNAComplex::write_results(){

	cerr << "RNA Complex ID:\t" << get_id() << endl;
	cerr << "RBS distance = " << get_rbs_distance() <<
			"\nTotal distance =  " << get_total_distance() << endl;

	cerr << "Sequence:\t" << get_sequence() <<
			"\nStructure:\t" << get_structure() <<
			"\nTraget:\t" << get_target_structure() << endl;
}

void RNAComplex::calculate_binding_fractions(){
/*
 * Calculates the binding among the RNA molecules in a complex
 *
 *
 */
	int			len;
	int			pos,comp_pos,i,j;
	bool		already_visited[MAX_CHAR];


	if(get_extra_rna_complex_flag() == true)
		return;

	len = get_sequence_length();

	if(len > MAX_CHAR){
		cerr << "ERROR: Increase the value of MAX_CHAR" << endl;
		exit(-1);
	}


	total_bound_sites = 0;
	self_bound_sites = 0;

	for(i = 0; i < total_rna_number; i++){
		rna_total_bound_sites[i] = 0;
		for(j = 0; j < total_rna_number; j++){
			rna_rna_bound_sites[i][j] = 0;
		}
	}

	for(pos = 0; pos < len; pos++)
		already_visited[pos] = false;


	for(pos = 0; pos < len; pos++){
		if(already_visited[pos] == true)
			continue;

		comp_pos = find_complementary_site(pos,get_structure());

		if(comp_pos == -1)
			continue;

		int id1 = get_rna_id_from_position(pos);
		int id2 = get_rna_id_from_position(comp_pos);

		if(id1 > id2){
			int tmp = id1;
			id1 = id2;
			id1 = tmp;
		}

		already_visited[pos] = true;
		already_visited[comp_pos] = true;

		rna_rna_bound_sites[id1][id2]++;
		total_bound_sites += 2;

		if(id1 == id2){
			self_bound_sites += 2;
			rna_total_bound_sites[id1] += 2;
		}
		else{
			rna_total_bound_sites[id1]++;
			rna_total_bound_sites[id2]++;
		}
	}

	self_bound_fraction = (double) self_bound_sites / total_bound_sites;
}

void RNAComplex::calculate_structure_distance(string &str,string &target,int *total_distance,int *rbs_distance){
/*
 * calculates rbs_distance: mismatch between the RBS regions of the target and current structure
 * and total_distance: total mismatch between the target and current structure
 *
 */
	if(get_extra_rna_complex_flag() == true)
		return;

	int	 dist = 0;
	int	 rbs_dist = 0;

	const char *tval = target.c_str();
	const char *sval = str.c_str();

	int len = (int) str.size();

	for(int i = 0; i < len; i++){

		if(sval[i] == '+' || sval[i] == ' ')
			continue;

		if(tval[i] == RBS_BOUND || tval[i] == RBS_FREE){
			if(((tval[i] == RBS_BOUND && (sval[i] == '(' || sval[i] == ')')) ||
			   (tval[i] == RBS_FREE && sval[i] == SITE_FREE)) )
				continue;

			rbs_dist++;
			continue;
		}

		if(!(tval[i] == ANYTHING ||
			tval[i] == sval[i]  ||
			(tval[i] == SITE_BOUND && (sval[i] == '(' || sval[i] == ')'))) )
			dist++;
	}

	*total_distance = dist + rbs_dist;
	*rbs_distance = rbs_dist;
}


void RNAComplex::calculate_structure_distance(int *total_distance,int *rbs_distance){
	calculate_structure_distance(get_structure(),get_target_structure(),total_distance,rbs_distance);
}


// ***************** Getters *******************

int RNAComplex::get_mrna_position(){
	return mrna_position;
}


double RNAComplex::get_self_bound_fraction(int rna_id){

	if(rna_id == -1)
		rna_id = mrna_position;

	return (double)rna_rna_bound_sites[rna_id][rna_id]/rna_total_bound_sites[rna_id];
}


double	RNAComplex::get_rna_rna_binding_fraction(int rna1_id,int rna2_id){
	if(rna2_id == -1)
		rna2_id = mrna_position;
	else if(rna1_id == rna2_id)
		return get_self_bound_fraction(rna1_id);

	if(rna1_id > rna2_id)
		return (double) rna_rna_bound_sites[rna1_id][rna2_id]/rna_total_bound_sites[rna2_id];
	else
		return (double) rna_rna_bound_sites[rna2_id][rna1_id]/rna_total_bound_sites[rna2_id];
}


int RNAComplex::get_rna_id_from_position(int pos){
	int 	k;

	for(k = 0; k < total_rna_number; k++){
		if(pos <= rna_position_end[k]){
			return k;
		}
	}

	cout << "\nERROR in get_rna_id_from_position with position =  " << pos << endl;
	cout << "total_rna_number " << total_rna_number << endl;
	cout << "length " << get_sequence_length() << endl;
	cout << "k = " << k << endl;

	exit(-1);
	return -1;
}

int RNAComplex::get_rna_starting_position_in_complex(int rna_id){
	return rna_position_base[rna_id];
}


int RNAComplex::calculate_toehold_length(bool force_any_rna_calculation){
/*
 * calculates toehold lengths. and generates the toehold site marked structure
 *
 */
	//  Brackets to identify the reaction initiation points....
	static const char toehold_pairs[6][2] = {{'{','}'},{'[',']'},{'<','>'},{'\\','/'},{'$','#'},{'&','%'}};

	if(!force_any_rna_calculation && get_extra_rna_complex_flag())
			return 0;

	int			pos;
	int 		comp_pos;
	int 		rna1_id;
	int			rna2_id;

	string 		&str = get_structure();
	int 		len =  get_sequence_length();
	bool 		already_visited[MAX_CHAR];


	total_toehold_length = 0;
	mean_toehold_length = 0;

	if(total_rna_number == 1)
		return 0;

	toehold_structure = str;

	for(int k = 0; k < total_rna_number; k++)
		rna_toehold_length[k] = 0;

	for(int i = 0; i < len; i++)
		already_visited[i] = false;

	for(int x = 0; x < len; x++){

		pos = 0 + x;

		if(already_visited[pos] == true)
			continue;

		if(str[pos] == '.' || str[pos] == '+' || str[pos] == ' ')
			continue;

		comp_pos = find_complementary_site(pos, str);
		rna1_id = get_rna_id_from_position(pos); //mrna_position;
		rna2_id = get_rna_id_from_position(comp_pos);

		already_visited[pos] = true;
		already_visited[comp_pos] = true;

		if(rna1_id == rna2_id)
			continue;

		int r1_base = get_rna_starting_position_in_complex(rna1_id);
		int r2_base = get_rna_starting_position_in_complex(rna2_id);

		if(rna[rna1_id]->get_structure()[pos - r1_base] == '.' &&
		   rna[rna2_id]->get_structure()[comp_pos - r2_base] == '.'){
			toehold_structure[pos] = toehold_pairs[rna1_id][0];
			toehold_structure[comp_pos] = toehold_pairs[rna1_id][1];

			rna_toehold_length[rna1_id]++;
			rna_toehold_length[rna2_id]++;
			total_toehold_length++;
		}
	}

	int rcount = 0;

	for(int i = 0; i < total_rna_number; i++){
		if(rna_toehold_length[i] != 0)
			rcount++;
	}
	mean_toehold_length = round((double)total_toehold_length / total_rna_number);

	return total_toehold_length;
}




string& RNAComplex::get_toehold_structure(){
	return toehold_structure;
}

int RNAComplex::get_total_rna_number(){
	return total_rna_number;
}

RNAObject& RNAComplex::get_rna(int rna_id){
	return *rna[rna_id];
}


int RNAComplex::get_total_toehold_length(){
	return total_toehold_length;
}

int RNAComplex::get_rna_toehold_length(int rna_id){
	return rna_toehold_length[rna_id];
}

int RNAComplex::get_mean_toehold_length(){
	return total_toehold_length;
}


double RNAComplex::get_min_rna_concentration(){
	return min_rna_concentration;
}

void RNAComplex::generate_mass_vector(int *mvector){
	int i;

	for(i = 0; i < MAX_RNA_NUMBER; i++)
		mvector[i] = 0;

	for(int i = 0; i < total_rna_number; i++){
		mvector[rna[i]->get_id()]++;
	}
}

double RNAComplex::get_kon_rate(){
	return kon_rate;
}

double RNAComplex::get_keq_rate(){
	return keq_rate;
}


bool RNAComplex::get_extra_rna_complex_flag(){
	return is_it_extra_rna_complex;
}

bool RNAComplex::get_target_state(){
	return target_state;
}


double RNAComplex::get_rbs_translation_rate(){
	return rbs_translation_rate;
}

// ***************** Setters *******************

void RNAComplex::set_mrna_position(int pos){
	mrna_position = pos;
}

void RNAComplex::set_rbs_target(bool flag){
/*
 * Generates RBS target sites according to the TRUTH table output given by the parameter flag.
 *
 * The RBS sites are given by symbols R, Q and P
 *
 * if flag = true,  R is replaced by F (RBS site Free)
 * if flag = false, R is replaced by B (RBS site Blocked)
 *
 * if flag = true,  Q is replaced by F (RBS site Free)
 * if flag = false, Q is replaced by @ (Anything)
 *
 * if flag = true,  P is replaced by @ (Anything)
 * if flag = false, P is replaced by B (RBS site Blocked)
 *
 */

	string str = get_target_structure();
	int	 len = (int) str.size();

	for(int k = 0; k < len; k++){
		if(str[k] == RBS_SITE){
			if(flag == true)
				str[k] = RBS_FREE;
			else
				str[k] = RBS_BOUND;
		} else if(str[k] == 'Q'){
			if(flag == true)
				str[k] = RBS_FREE;
			else
				str[k] = ANYTHING;
		} else if(str[k] == 'P'){
			if(flag == true)
				str[k] = ANYTHING;
			else
				str[k] = RBS_BOUND;
		}

	}
	set_target_structure(str);
}


void RNAComplex::set_target_state(bool flag){
	target_state = flag;
}


void RNAComplex::reset_mutation_flag(){
	for(int i = 0; i < total_rna_number; i++){
		rna[i]->set_mutation_flag(false);
	}
}


void RNAComplex::reset_component_rna_id(){
	for(int i = 0; i < total_rna_number; i++)
		rna[i]->set_id(i);
}

void RNAComplex::save_parameters(){
	for(int k = 0; k < total_rna_number; k++){
		rna[k]->save_parameters();
	}
	RNAObject::save_parameters();
}


void RNAComplex::restore_parameters(){
	RNAObject::restore_parameters();
	for(int k = 0; k < total_rna_number; k++){
		rna[k]->restore_parameters();
	}
	update_rnacomplex_sequence_and_structure();
}


void RNAComplex::set_extra_rna_complex_flag(bool cflag){
	is_it_extra_rna_complex = cflag;
}



