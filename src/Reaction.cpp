/*
 * Reaction.cpp
 *
 *  Created on: Aug 30, 2010
 *      Author: Marzuk M Kamal
 */

#include "Reaction.h"
#include <ctime>


using namespace std;

int truth_table_OR_gate[14] = {4, 3, // row and column
							0, 0, 0,
							0, 1, 1,
							1, 0, 1,
							1, 1, 1};

int truth_table_NOR_gate[14] = {4, 3, // row and column
							0, 0, 1,
							0, 1, 0,
							1, 0, 0,
							1, 1, 0};


int truth_table_AND_gate[14] = {4, 3,
							0, 0, 0,
							0, 1, 0,
							1, 0, 0,
							1, 1, 1};


int truth_table_NAND_gate[14] = {4, 3,
							0, 0, 1,
							0, 1, 1,
							1, 0, 1,
							1, 1, 0};


int truth_table_XOR_gate[14] = {4, 3,
							0, 0, 0,
							0, 1, 1,
							1, 0, 1,
							1, 1, 0};


int truth_table_XNOR_gate[14] = {4, 3,
							0, 0, 1,
							0, 1, 0,
							1, 0, 0,
							1, 1, 1};


int truth_table_YES_gate[6] = {2, 2,
							   0, 0,
							   1, 1};

int truth_table_NOT_gate[6] = {2, 2,
							   0, 1,
							   1, 0};



Reaction::Reaction(){
	id = 0;
	total_truth_table_rna_complex = 0;
	extra_rna_complex = 0;
	max_iteration = 1000;
	total_free_energy = 0;
	main_rna_complex = NULL;
	output_mode = 3;	// 1 | 2 = 3,  write to both stdout and file;
	weight = 1;
	calculation_mode = MONTE_CARLO_MODE;
	mc_steps_count = 0;
	total_mc_count = 0;
	total_cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

	omp_set_dynamic(false); // de-activate dynamics threading... will check it later
//	omp_set_dynamic(true); // activate dynamics threading... will check it later
	omp_set_nested(true);  // nested parallelism
	omp_set_num_threads(total_cpu_count*2);
//
}

Reaction::Reaction(RNAComplex rnac, RNAComplex target){
	Reaction();
}

Reaction::Reaction(RNAComplex &rnac){
	Reaction();
	add_rna_complex(rnac);
}


Reaction::~Reaction(){
}


void Reaction::create_rna_complex_from_truth_table(RNAComplex *full_rna_complex,int *truth_table,int mrna_pos,bool place_mrna_first){
/*
 * Generates RNA complexes according to the truth table.
 *
 * NOTE: Do not delete/free any of the input parameters outside of this function. These data are required during
 * 		 the running time of the RNAGate!
 *
 * full_rna_complex = this RNAComplex object contains all the RNAs that will take part in the reaction.
 * truth_table = 	  the target truth table.
 *
 * This function also makes some combinations of RNAs which is required for concentration calculation
 * For example,
 *
 * if we have 2 sRNA and 1 UTR and Truth table is AND gate
 * the RNA combination will be,
 *
 *  UTR
 *  UTR+sRNA1
 *  UTR+sRNA2
 *  UTR+sRNA1+sRNA2
 *
 *  and
 *
 *  sRNA1+sRNA1
 *  sRNA2+sRNA2
 *  sRNA1+sRNA2
 *
 * The combinations of sRNAs is required solely for concentration calculations.
 *
 */

	int 		row, col;
	RNAComplex *rnac   = NULL;


	int			base= 2;
	int			other_rna_idx[50];

	cout << "Creating RNA Complex from the Truth Table..." << endl;

	if(truth_table == NULL){
		cout << "truth_table == NULL" << endl;
		exit(-1);
	}

	if(full_rna_complex == NULL){
		cout << "full_rna_complex == NULL" << endl;
		exit(-1);
	}

	main_rna_complex = full_rna_complex;

	full_rna_complex->reset_component_rna_id();

	row = truth_table[0];
	col = truth_table[1];

	int k = 0;
	for(int i = 0; i < full_rna_complex->get_total_rna_number(); i++){
		if(i != mrna_pos){
			other_rna_idx[k++] = i;
		}
	}

	for(int i = 0; i < row; i ++){
		rnac = new RNAComplex;

		if(rnac == NULL){
			cerr << "ERROR: Cannot allocate memory for creating new RNA Complex." << endl;
			exit(-1);
		}

		if(place_mrna_first == true){
			rnac->add_rna(full_rna_complex->get_rna(mrna_pos));
			rnac->set_mrna_position(0);
		}

		for(int j = 0; j < col-1; j++){
			if(truth_table[base + i*col + j]){
				rnac->add_rna(full_rna_complex->get_rna(other_rna_idx[j]));
			}
		}

		if(place_mrna_first == false){
			rnac->add_rna(full_rna_complex->get_rna(mrna_pos));
			rnac->set_mrna_position(full_rna_complex->get_total_rna_number()-1);
		}

		rnac->set_id(i);
		rnac->update_target_parameters();
		rnac->update();

		rnac->set_rbs_target(truth_table[base + (i + 1) * col - 1]);

		int	rbs_dist, total_dist;

		rnac->calculate_structure_distance(&total_dist,&rbs_dist);
		rnac->set_rbs_distance(rbs_dist);
		rnac->set_total_distance(total_dist);

		add_rna_complex(*rnac);
	}

	calculate_shared_rna_count();

//********** This following part is necessary for concentration calculation.... ***************
	extra_rna_complex= 0;
	int complex_count = 0;

	int total_rna = full_rna_complex->get_total_rna_number();

	cout << "Creating additional RNA complexes for concentration calculation..." << endl;

	complex_count = 0;

	// add sRNA dimers...
	for(int ix = 0; ix < total_rna; ix++){
		if(ix == mrna_pos)
			continue;

		rnac = new RNAComplex;

		if(rnac == NULL){
			cerr << "ERROR: cannot allocate memory for rnac while generating Additional RNA Complex." << endl;
			exit(-1);
		}

		rnac->add_rna(full_rna_complex->get_rna(ix));
		rnac->add_rna(full_rna_complex->get_rna(ix));
		rnac->set_extra_rna_complex_flag(true);

		rnac->set_id(total_truth_table_rna_complex + complex_count++);
		rnac->update_target_parameters();
		rnac->update();
		add_extra_rna_complex(*rnac);
	}

	// add sRNA dimer combinations... e.g. sRNA1+sRNA2
	for(int ix = 0; ix < total_rna; ix++){
		if(ix == mrna_pos)
			continue;

		for(int iy = ix+1; iy < total_rna; iy++){
			if(iy == mrna_pos)
				continue;

			rnac = new RNAComplex;

			if(rnac == NULL){
				cerr << "ERROR: cannot allocate memory for rnac while generating Additional RNA Complex." << endl;
				exit(-1);
			}

			rnac->add_rna(full_rna_complex->get_rna(ix));
			rnac->add_rna(full_rna_complex->get_rna(iy));
			rnac->set_extra_rna_complex_flag(true);

			rnac->set_id(total_truth_table_rna_complex + complex_count++);
			rnac->update_target_parameters();
			rnac->update();
			add_extra_rna_complex(*rnac);
		}
	}


	// add the RNA monomers except the 5'UTR....
	for(int ix = 0; ix < total_rna; ix++){
		if(ix == mrna_pos)
			continue;

		rnac = new RNAComplex;

		if(rnac == NULL){
			cerr << "ERROR: cannot allocate memory for rnac while generating Additional RNA Complex." << endl;
			exit(-1);
		}

		rnac->add_rna(full_rna_complex->get_rna(ix));
		rnac->set_extra_rna_complex_flag(true);

		rnac->set_id(total_truth_table_rna_complex + complex_count++);
		rnac->update_target_parameters();
		rnac->update();
		add_extra_rna_complex(*rnac);
	}

	for(int i = 0;i < total_truth_table_rna_complex+extra_rna_complex;i++){
		cout << "RNAComplex: " << rna_complex[i]->get_id() << " (" << rna_complex[i]->get_total_rna_number() << ")";

		for(int j = 0; j < rna_complex[i]->get_total_rna_number(); j++){
			cout << "\t" << rna_complex[i]->get_rna(j).get_id();
		}
		cout << endl;
	}

	cout << "Total RNA:\t" << total_truth_table_rna_complex << endl <<
			"Extra RNA:\t" << extra_rna_complex << endl;

	int 	amatrix_complexes[100];
	int		amatrix_columns;
	int		amatrix_rows;
	int		index_vector[100];

	int		**Am;
	int		*cpos_vector;
	double	*rna_conc;
	double	*complex_conc;
	int 	main_cpos;

	for(int i = 0; i < 100; i++)
		index_vector[i] = 0;


	for(int i = 0; i < total_truth_table_rna_complex; i++){
		if(!rna_complex[i]){
			cout << "rna_complex[" << i << "] is NULL" << endl;
			continue;
		}

		amatrix_rows = rna_complex[i]->get_total_rna_number();
		if(amatrix_rows <= 1)
			continue;

		Am = new int*[amatrix_rows];

		get_combinatorial_complexes(i,amatrix_complexes,amatrix_columns,main_cpos);
		cpos_vector = new int[amatrix_columns];
		memcpy(cpos_vector,amatrix_complexes,amatrix_columns*sizeof(int));

		cout << "i:\t" << i << "\nrows: " << amatrix_rows << "col: " << amatrix_columns << endl;

		cout << "cpos_vector: ";
		for(int p = 0; p < amatrix_columns;p++){
			cout << cpos_vector[p] << " ";
		}
		cout << endl;

		cout << "combinatorial complexes: ";
		for(int n = 0; n < amatrix_columns; n++)
			cout << amatrix_complexes[n] << "\t";
		cout << endl;

		for(int j = 0; j < amatrix_rows; j++){
			Am[j] = new int[amatrix_columns];
			memset(Am[j], 0, amatrix_columns*sizeof(int));
			index_vector[rna_complex[i]->get_rna(j).get_id()] = j;
		}

		for(int k = 0; k < amatrix_columns; k++){
			RNAComplex *c1 = rna_complex[amatrix_complexes[k]];

			for(int p = 0; p < c1->get_total_rna_number(); p++){
				Am[index_vector[c1->get_rna(p).get_id()]][k]++;
			}
		}

		cout << "Created Mass matrix for RNAComplex " << i << ":"<< endl;

		rna_conc = new double[amatrix_rows];
		for(int row = 0; row < amatrix_rows; row++){
			rna_conc[row] = rna_complex[main_cpos]->get_rna(row).get_concentration();

			for(int col = 0; col < amatrix_columns; col++){
				cout << Am[row][col] << "\t";
			}
			cout << endl;
		}

		for(int j = 0; j < amatrix_rows; j++){
			cout << "row: " << j << " idx: " << index_vector[rna_complex[main_cpos]->get_rna(j).get_id()] << " conc: " << rna_conc[j] << endl;
		}

		cout << "Main position: " << main_cpos << endl;

		concentration_obj[i] = new Concentration;
		concentration_obj[i]->init_parameters(amatrix_rows,amatrix_columns,rna_conc,Am,cpos_vector,main_cpos);
		concentration_obj[i]->print_A_matrix();
	}

	for(int ix = 0; ix < total_truth_table_rna_complex; ix++){
		cout << "Complex: " << ix <<
				"\nSequence:\t"  << rna_complex[ix]->get_sequence() <<
				"\nStructure:\t" << rna_complex[ix]->get_structure() <<
				 "\nTarget:\t\t" << rna_complex[ix]->get_target_structure() << endl;
	}
}

void Reaction::calculate_shared_rna_count(){
/*
 * counts distinct number of RNA molecules involved in the RNAGate calculation.
 *
 */

	bool	rflag[MAX_RNA_NUMBER];
	int 	rna_count = 0;

	for(int i = 0; i < MAX_RNA_NUMBER; i++)
		rflag[i] = false;


	for(int i = 0; i < total_truth_table_rna_complex; i++){
		for(int j = 0; j < rna_complex[i]->get_total_rna_number(); j++){
			int idx = rna_complex[i]->get_rna(j).get_id();
			if(rflag[idx] == false){
				rflag[idx] = true;
				rna_count++;
			}
		}
	}

	total_rna = rna_count;
}


void Reaction::add_rna_complex(RNAComplex &rnac){
	rna_complex[total_truth_table_rna_complex++] = &rnac;
}


void Reaction::add_extra_rna_complex(RNAComplex &rnac){
	rna_complex[total_truth_table_rna_complex + extra_rna_complex++] = &rnac;
}


void Reaction::rna_complex_random_select(int rnac_num, int *rnac_indeces){
	bool 	ri[MAX_RNA_NUMBER];
	int		i;
	int		val;

	for(i = 0; i < total_truth_table_rna_complex; i++){
		ri[i] = false;
	}

	i = 0;
	while(true){
		if(i == rnac_num)
			break;

		val = floor(random_num() * total_truth_table_rna_complex);

		if(ri[val] == false){
			ri[val] = true;
			rnac_indeces[i++] = val;
		}
	}
}

double	Reaction::calc_score_kin(){
/*
 * Calculates the kinetic score function that was used in RNAdes... this scoring function is NOT used in this program
 * for testing purpose only...
 */

	double 	dg;
	double 	sc_kin = 0;
	int		tl;

	for(int i = 0; i < total_truth_table_rna_complex; i++){
		if(rna_complex[i]->get_total_rna_number() > 1){
			dg = rna_complex[i]->get_energy_difference();
			tl = rna_complex[i]->get_total_toehold_length();

			if(abs(dg) < DG_PARAM)
				sc_kin += DG_PARAM - abs(dg);

			if(tl < TOEHOLD)
				sc_kin += (TOEHOLD - tl) * DG_BP;
		}
		else{
			sc_kin += abs(dg);

			sc_kin += tl * DG_BP;
		}
	}

	return sc_kin/total_truth_table_rna_complex;
}

double Reaction::calc_score_struct(){
/*
* Calculates the structural score function that was used in RNAdes... this scoring function is NOT used in this program
* for testing purpose only...
*/
	double sc_str = 0;

	for(int i = 0; i < total_truth_table_rna_complex; i++){
		sc_str += rna_complex[i]->get_total_distance();
	}

	return sc_str/total_truth_table_rna_complex;
}

double Reaction::calculate_total_rnagate_score(double &prob_score, double &conc_score){
/*
 * Calculates the the scoring function Gs of the complexes in the truth table.
 *
 * free_energy = total free energy gain in reaction - mean structural probability - mean concentration of the rna complexes.
 *
 */

	double	free_energy = 0;
	double  pscore = 0;
	double  cscore = 0;
	double 	rna_pscore = 0;

	int 	rc_count = 0;
	for(int i = 0; i < total_truth_table_rna_complex; i++){
		free_energy += 0*rna_complex[i]->get_energy_difference() - rna_complex[i]->get_total_toehold_length() * DG_EFF;

		if(rna_complex[i]->get_total_rna_number() > 1){
			pscore += rna_complex[i]->get_probability();
			rc_count++;
		}

		cscore += get_concentration_fraction(i);
	}

	double l1_rna = 100;
	double l1_rna_complex = 100;
	double l2 = 100;

	int rna_count = main_rna_complex->get_total_rna_number();

	for(int j = 0; j < rna_count; j++){
		rna_pscore += main_rna_complex->get_rna(j).get_probability();
	}

//	int 	total_rna_n_complex = rc_count + rna_count;
//	double 	total_pscore = (pscore + rna_pscore) / total_rna_n_complex;

	//	int 	total_rna_n_complex = rc_count + rna_count;
	free_energy = free_energy/rc_count/KT;
	pscore = l1_rna_complex * pscore/rc_count;
	rna_pscore = l1_rna * rna_pscore/rna_count;
	cscore = l2 * cscore /rc_count;


	double total_val = (-free_energy) + pscore + rna_pscore + cscore;

//	double 	total_pscore = l1_rna_complex * pscore/rc_count + l1_rna * rna_pscore/rna_count;



//	int 	total_rna_n_complex = rc_count + rna_count;
	double 	total_pscore = l1_rna_complex * pscore/rc_count + l1_rna * rna_pscore/rna_count;

	if(score_weighting_mode == ADAPTIVE_SCORE_WEIGHTING){
		double w_f =  1 + free_energy/total_val;
		double w_p_rc = 1 - pscore/total_val;
		double w_p_rna = 1 - rna_pscore/total_val;
		double w_conc = 1 - cscore/total_val;

		total_free_energy = w_f * free_energy - w_p_rc * pscore - w_p_rna * rna_pscore - w_conc * cscore;
	}
	else{
		total_free_energy = free_energy -  pscore -  rna_pscore - cscore;
	}

	prob_score = total_pscore;
	conc_score = cscore;

	return total_free_energy;
}


bool Reaction::evolve(int word_len){
/*
 * This is where almost everything happens!
 *
 * evolve randomly mutates the component RNA molecules and calculates their and corresponding structures
 * of the rna complexes.
 *
 * It primarily accepts the newly mutated complexes only if
 *
 * 1. the distances (or mismatches) between the RBS sites of ALL the new complexes and the target structures either
 *    remain the same or decrease. The optimal situation would be RBS distance = 0.
 * 2. the distance of the remaining parts of the complexes from their target is within a given range. (not implemented yet!)
 *
 */
//	double		p;
	bool		do_evolve = true;
	int			max_iter = 1000;
	double		old_energy = 0;
	string		str;
	string		seq;
	int			rbs_dist,total_dist;
	int			i,k;
	int	 		rnac_number = total_truth_table_rna_complex;
	int  		rnac_index[30];

	int			chunk = 4*total_cpu_count;

	double  	new_free_energy = 0;
	double		old_free_energy = 0;


	int count  = 0;

	double old_score = 9999;


	double 		new_complex_mean_pvalue = 0;
	double		old_complex_mean_pvalue = 0;

	double 		new_total_probability = 0;
	double 		old_total_probability;

	double		old_pscore;
	double		old_cscore;

	bool 		mode_set = false;

	while(count++ < max_iter){
		double	rfree;
		double  rbound;

		calculate_total_rbs_distance(rfree,rbound);

		if(get_calculation_evolution_mode() == EVOLVE_FOR_STRUCTURE
				&& rfree <= RBS_FREE_THRESHOLD
				&& rbound <= RBS_BOUND_THRESHOLD
		){

			cout << "Calculation mode changed to EVOLVE_FOR_PROBABILITY\n";

			set_calculation_evolution_mode(EVOLVE_FOR_PROBABILITY);
			int i;

			for(int i = 0; i < main_rna_complex->get_total_rna_number();i++){
				main_rna_complex->set_probability(0);
				main_rna_complex->set_pfunc_z_value(1);
				main_rna_complex->set_concentration(1);
			}

			for(int i = 0; i < total_truth_table_rna_complex+extra_rna_complex;i++){
				rna_complex[i]->set_probability(0);
				rna_complex[i]->set_pfunc_z_value(1);
				rna_complex[i]->set_concentration(0);
			}
		}

		old_complex_mean_pvalue = calculate_complex_mean_pvalue();
		old_free_energy = calculate_total_rnagate_score(old_pscore,old_cscore);

		total_mc_count++;

		old_total_probability = 0;
		RNAObject::generate_unique_random_id(rnac_number,rnac_index);

		int 	all_accepted = true;
		int		rcount = 0;
		int		j;
		double 	new_pvalue;
		double 	pvalue;

		chunk = 1;

		for(i = 0; i < rnac_number; i++){
			k = rnac_index[i];									// chose a complex randomly

			for(int n = 0; n < total_truth_table_rna_complex; n++){
				rna_complex[n]->reset_mutation_flag();
				rna_complex[n]->save_parameters();
			}

			old_energy = rna_complex[k]->get_energy();

			do_evolve = rna_complex[k]->mutate(word_len,get_calculation_evolution_mode());		///***** mutate rna_complex[k] *******/

			if(do_evolve == false){
				for(int n = 0; n < total_truth_table_rna_complex; n++){
					rna_complex[n]->restore_parameters();
					rna_complex[n]->reset_mutation_flag();
				}
				continue;
			}

			all_accepted = true;
			rcount = 0;
			pvalue = 0;

			for(j = 0; j < total_truth_table_rna_complex; j++){			// check rest of the rna_complexes
				if(j == k)
					continue;

				if(check_complexes_with_shared_rna(k,j) == false) // skip the complexes whose RNA is not modified
					continue;
				else
					rna_complex[j]->set_mutation_flag(true);

				pvalue = rna_complex[j]->get_probability();

				rna_complex[j]->update(get_calculation_evolution_mode());
				rna_complex[j]->calculate_structure_distance(&total_dist,&rbs_dist);

				if(get_calculation_evolution_mode() == EVOLVE_FOR_PROBABILITY)
					new_complex_mean_pvalue = calculate_complex_mean_pvalue();
				else
					new_complex_mean_pvalue = old_complex_mean_pvalue;


				if(rbs_dist <= rna_complex[j]->get_rbs_distance()
						//				   && (new_pvalue >= 0.5 || new_pvalue >= pvalue)
						//							&& (new_complex_mean_pvalue >= 0.5*total_truth_table_rna_complex || new_complex_mean_pvalue >= old_complex_mean_pvalue)
				){
					// are the RBS distance remains same or approaching zero

					old_complex_mean_pvalue = new_complex_mean_pvalue;

					all_accepted &= true;
					rcount++;

					rna_complex[j]->set_rbs_distance(rbs_dist);
					rna_complex[j]->set_total_distance(total_dist);
				}
				else{
					//					rna_complex[j]->set_probability(pvalue);
					all_accepted &= false;
					break;
				}
			}



			if(all_accepted == true){								// are all the rna_complexes converging?
				if(get_calculation_evolution_mode() == EVOLVE_FOR_PROBABILITY)
					calculate_all_complex_concentration();

				double	new_pscore;
				double	new_cscore;

				new_free_energy =  calculate_total_rnagate_score(new_pscore,new_cscore);
				double de = new_free_energy - old_free_energy;


				if(de <= 0){
						return true;
				}
				else{
					double p = exp(-de * get_weight());

					if(p > random_num()){
						cout << "\nde = " << de <<  ", free energy = " << new_free_energy <<  " old_free_energy = " << old_free_energy << endl;
						cout << "ACCEPTED with p = " << p << endl;

						return true;
					}
				}
			}

			for(int i = 0; i < total_truth_table_rna_complex; i++){
				rna_complex[i]->restore_parameters();
				rna_complex[i]->reset_mutation_flag();
			}
		} // for(i = 0....)
	} // while(....
	return false;
}


void Reaction::run(){
/*
 * Starting point of the simulation
 *
 * The calculation mode is either  Monte Carlo or Simulated Annealing depending on the mode set by
 * set_calculation_mode(cmode);
 *
 * also the function writes to either stdout or file, or both depending on the command line option.
 *
 */
	int			word_len = 3;
	bool		do_evolve = false;
	int			count = 0;
	string		energy_out_file_name;
	string		stmp = output_file_name;
	int			data_count = 0;
	double		sc_kin;
	double		sc_str;
	double		sw;


	char	host_computer_name[200];
	time_t	prog_time;
	const struct tm *ltime;
	time(&prog_time);
	ltime = localtime(&prog_time);

	gethostname(host_computer_name,sizeof(host_computer_name));

	int pos = stmp.find('.');

	energy_out_file_name = stmp.substr(0,pos) + "_energy.txt";

	output.open(output_file_name);

	if(output.is_open() == false){
		cerr << "ERROR: Cannot create output file " << output_file_name <<  endl;
		exit(-1);
	}

	output << "RNAGate v11\n(c) 2011 Marzuk Kamal\nCalculation with command line instructions:" << endl;
	for(int i = 0; i < ARGC_VALUE; i++){
		output << ARGV_VALUE[i] << " ";
	}
	output << "\nProgram run in system: " << host_computer_name << "\nDate: " << asctime(ltime) << endl;
	output << "*******************************************************\n\n\n";


	energy_output.open(energy_out_file_name.c_str());

	if(energy_output.is_open() == false){
		cerr << "ERROR: Cannot create output file " << energy_out_file_name <<  endl;
		exit(-1);
	}

	cout << "Starting Monte Carlo optimization..." << endl;
//	write_results();

	energy_output << "Time\t" << "N\t" << "MC Step\t";
	for(int i = 0; i < main_rna_complex->get_total_rna_number(); i++){
		energy_output << "RNA_" << i
					  << "\tPr_" << i
				      << "\t";
	}
	for(int i = 0; i < total_truth_table_rna_complex; i++){
		energy_output << "Complex_" << i
					  <<"\tDGc_" << i
					  << "\tToehold_" << i
					  << "\tDistance_" << i
					  << "\tConc_" << i
					  << "\tKon_" << i
					  << "\tKeq_" << i
					  << "\tCxPr_" << i
					  << "\t";
	}
	energy_output << "Total_DE" << "\tScore" << "\tSc_Kin" << "\tSc_dist" <<endl;

	time(&start_time);

	double	dw = 1.0/(1000.0);
	double	weight_value = 0;

	sw = 0;
	int start_writing = max_iteration*80/100;

	score_kin  = calc_score_kin();
	score_dist = calc_score_struct();

	score_weight = 0.5;//score_kin/(score_kin + score_dist);
	sscore = score_weight * score_kin + (1 - score_weight) * score_dist;

	time_t	current_time;

	set_calculation_evolution_mode(EVOLVE_FOR_STRUCTURE);

	while(count++ < max_iteration){

		if(get_calculation_mode() == SIMULATED_ANNEALING_MODE){
			if(weight_value < 1.0)
				weight_value += dw;
			else
				weight_value = 1.0;

			set_weight(weight_value);
		}

		do_evolve = evolve(word_len); /**************************/

		if(do_evolve == true){
			mc_steps_count = count;
			cout << "Count =\t" << count << endl;

			if(get_calculation_evolution_mode() == EVOLVE_FOR_PROBABILITY){
				if(count > 0*start_writing){
					if(output_mode & WRITE_TO_FILE)
						write_results(&output);

					if(output_mode & WRITE_TO_STDOUT)
						write_results(NULL);
				}
			}
			else if(output_mode & WRITE_TO_STDOUT){
						write_results(NULL);
			}

			time(&current_time);
			double	time_val = difftime(current_time,start_time)/60.0;

			energy_output <<  setprecision(2) << time_val << "\t" << data_count++ << "\t" << count << "\t";
			for(int i = 0; i < main_rna_complex->get_total_rna_number(); i++){
				energy_output << main_rna_complex->get_rna(i).get_energy() << "\t"
							  << main_rna_complex->get_rna(i).get_probability()
							  << "\t";
			}

			for(int i = 0; i < total_truth_table_rna_complex; i++){
				energy_output << rna_complex[i]->get_energy() << "\t"
							  << rna_complex[i]->get_energy_difference() << "\t"
							  << rna_complex[i]->get_total_toehold_length() << "\t"
							  << rna_complex[i]->get_rbs_distance() << "\t"
							  << rna_complex[i]->get_concentration() << "\t"
							  << rna_complex[i]->get_kon_rate() << "\t"
							  << rna_complex[i]->get_keq_rate() << "\t"
							  << rna_complex[i]->get_probability() << "\t";
			}
			energy_output << get_total_free_energy() << "\t"
							 << sscore << "\t"
							 << score_kin << "\t"
							 << score_dist
							 << endl;

		}
	}

	time(&end_time);

	total_time_taken = difftime(end_time,start_time)/ 60; // minutes

	if(output_mode & 1)
		output << "Total time taken:\t" << total_time_taken << endl;

	if(output_mode & 2)
		cout << "Total time taken:\t" << total_time_taken << endl;

	energy_output.close();
	output.close();
}


bool Reaction::check_complexes_with_shared_rna(int complex1, int complex2){
/*
 * Checks whether complex1 and complex2 two are sharing common RNA
 * if an RNA in complex1 is mutated and also present in complex2
 * then the function returns true
 *
 * otherwise, it's false.
 *
 */
	RNAComplex *c1 = rna_complex[complex1];
	RNAComplex *c2 = rna_complex[complex2];

	for(int i = 0; i < c1->get_total_rna_number(); i++){
		if(c1->get_rna(i).get_mutation_flag() == true){
			int rna_id = c1->get_rna(i).get_id();

			for(int j = 0; j < c2->get_total_rna_number(); j++){
				if(rna_id == c2->get_rna(j).get_id())
					return true;
			}
		}
	}
	return false;
}


bool Reaction::check_combinatorial_complexes(int complex1, int complex2){
/*
 * Check whether complex2 is a contains common RNAs from complex1
 *
 * if complex2 contains RNAs that are not present in complex1
 * function returns false
 *
 * otherwise functions returns true;
 *
 */
	bool		rna_does_not_exist;

	RNAComplex *c1 = rna_complex[complex1];
	RNAComplex *c2 = rna_complex[complex2];

	for(int i = 0; i < c2->get_total_rna_number(); i++){
		int rna_id = c2->get_rna(i).get_id();
		rna_does_not_exist = true;
		for(int j = 0; j < c1->get_total_rna_number(); j++){
			rna_does_not_exist &= (rna_id != c1->get_rna(j).get_id());
		}

		if(rna_does_not_exist == true)
			return false;
	}

	return true;
}


void Reaction::get_combinatorial_complexes(int complex1, int *idx_data, int &idx_len,int &main_cpos){
	int count = 0;

	for(int i = 0; i < total_truth_table_rna_complex+extra_rna_complex; i++){
		if(check_combinatorial_complexes(complex1,i) == true){
			if(complex1 == i)
				main_cpos = count;
			idx_data[count++] = rna_complex[i]->get_id();
		}
	}
	idx_len = count;
}

void Reaction::create_vector_from_complex(int complex1,int *complex_vec){
/*
 *
 *
 */
	for(int i = 0; i < rna_complex[complex1]->get_total_rna_number(); i++){
		complex_vec[rna_complex[complex1]->get_rna(i).get_id()] = 0;
	}

	for(int i = 0; i < rna_complex[complex1]->get_total_rna_number(); i++){
		complex_vec[rna_complex[complex1]->get_rna(i).get_id()]++;
	}
}

void Reaction::calculate_all_complex_concentration(bool force_all_calculation){
	double		en_buff[100];
	int 		k;
	int			x;
	int 		i;
	int 		chunk = total_cpu_count;
	int			complex_mut_flag[100];

#pragma omp parallel private( maxGapIndex,k,x)
	{
#pragma omp for schedule(static,chunk) nowait
		for(k = 0; k < total_truth_table_rna_complex+extra_rna_complex; k++){
			if(k >= total_truth_table_rna_complex){
				for(x = 0; x < rna_complex[k]->get_total_rna_number(); x++){
					if(force_all_calculation || rna_complex[k]->get_rna(x).get_mutation_flag() == true){
						rna_complex[k]->update(EVOLVE_FOR_STRUCTURE);
						break;
					}
				}
			}
			en_buff[k] = rna_complex[k]->get_energy()/KT;
			en_buff[k] = (rna_complex[k]->get_energy() + rna_complex[k]->calculate_toehold_length(true)*DG_EFF)/KT;

		}
	}

#pragma omp parallel private(i)
	{
#pragma omp for schedule(static,chunk) nowait
		for(i = 0; i < total_truth_table_rna_complex; i++){
			if(concentration_obj[i] ){
				if(force_all_calculation || rna_complex[i]->get_mutation_flag() == true)
				{
					concentration_obj[i]->update_complex_free_energy_vector(en_buff);
					concentration_obj[i]->calculate_concentration();
					rna_complex[i]->set_concentration(get_concentration(i));
				}
			}
		}
	}
}


void Reaction::calculate_total_rbs_distance(double &rbs_free_dist,double &rbs_bound_distance){
/*
 *
 *
 */
	int	free_val = 0;
	int	bound_val = 0;

	int 	rna_free_count = 0;
	int 	rna_bound_count = 0;


	for(int i = 0; i < total_truth_table_rna_complex; i++){
		int dst, tot_site;
		rna_complex[i]->get_total_free_distance(dst,tot_site);

		free_val += dst;
		rna_free_count += tot_site;

		rna_complex[i]->get_total_bound_distance(dst,tot_site);

		bound_val += dst;
		rna_bound_count += tot_site;
	}

	if(rna_free_count == 0)
		rna_free_count = 1;
	if(rna_bound_count == 0)
		rna_bound_count = 1;

	rbs_free_dist = (double)free_val/rna_free_count;
	rbs_bound_distance = (double)bound_val/rna_bound_count;
}


int Reaction::sequence_distance(string &seq1,string &seq2){

	return 0;
}


double Reaction::sequence_alignment_score(string &seq1,string &seq2){
	double	score;

	return score;
}

int Reaction::calculate_2D_sequence_map(string *seq, int total_seq){
	return 0;
}

/********************** Output **************************/

void Reaction::write_output(){
	for(int k = 0; k < total_truth_table_rna_complex; k++){
		rna_complex[k]->write_results();
	}
}

void Reaction::write_results(ofstream *fptr){
/*
 * Display results and write to file!
 *
 */
	  streambuf *psbuf, *backup;

	  if(fptr == NULL){
		  write_results();
		  return;
	  }

	  backup = cout.rdbuf();
	  psbuf = fptr->rdbuf();
	  cout.rdbuf(psbuf);
	  write_results();
	  cout.rdbuf(backup);
}


void Reaction::write_results(){
	string	rna_current_structure;
	string	str_mono;

	for(int i = 0; i < total_truth_table_rna_complex; i++){
		for(int j = 0; j < rna_complex[i]->get_total_rna_number(); j++){
			if(j == 0){
				str_mono = rna_complex[i]->get_rna(j).get_structure();
			}else{
				str_mono += " " + rna_complex[i]->get_rna(j).get_structure();
			}
		}

		cout << "MC Step:\t" << mc_steps_count << "(" << total_mc_count << ")" <<
				"\nRNA Complex:\t" << i <<
				"\nSequence:\t" << rna_complex[i]->get_sequence() <<
				"\nStrConstr:\t" << rna_complex[i]->get_structure_constraint() <<
				"\nMonoStructure:\t" << str_mono <<
				"\nStructure:\t" << rna_complex[i]->get_structure() <<
				"\nToeStruct:\t" << rna_complex[i]->get_toehold_structure() <<
				"\nTarget:\t\t" << rna_complex[i]->get_target_structure() <<
				"\nRBS Distance:\t" << rna_complex[i]->get_rbs_distance() <<
				"\nTotDistance:\t" << rna_complex[i]->get_total_distance() <<
				"\nKon:\t" << rna_complex[i]->get_kon_rate() <<
				"\nKeq:\t" << rna_complex[i]->get_keq_rate() <<
				"\nConcentration:\t" << get_concentration(i) <<
				"\nProbability:\t" << rna_complex[i]->get_probability() << endl;



		cout << "Energy: ";
		for(int k = 0; k < rna_complex[i]->get_total_rna_number(); k++){
			cout << "RNA " << k << " = " << rna_complex[i]->get_rna(k).get_energy() << "(pb: "
					<< rna_complex[i]->get_rna(k).get_energy_per_base() << ", "<<  rna_complex[i]->get_rna(k).get_probability()<<  "); ";
		}

		cout <<	" Complex: " << rna_complex[i]->get_energy() << "(pb: " << rna_complex[i]->get_energy_per_base()<< ", "
				<<  rna_complex[i]->get_probability() << ")" << endl;
		cout << "Toehold: Total = " << rna_complex[i]->get_total_toehold_length() << "; Average = " << rna_complex[i]->get_mean_toehold_length() << endl;
	}
}


double calculate_concentration(int complex_id){
	double conc_value;


	return conc_value;
}

/************************* Getters **************************/

int	Reaction::get_id(){
	return id;
}

double Reaction::get_energy(){
	return energy;
}

double Reaction::get_concentration(){
	return concentration;
}

double Reaction::get_concentration(int complex_idx){
	if(complex_idx >= total_truth_table_rna_complex){
		cout << "ERROR: In Reaction::get_concentration, complex_idx is greater than total_truth_table_rna_complex." << endl;
		exit(-1);
	}

	if(rna_complex[complex_idx]->get_total_rna_number() <= 1)
		return rna_complex[complex_idx]->get_rna(0).get_concentration();

	return concentration_obj[complex_idx]->get_concentration();
}


double Reaction::get_concentration_fraction(int complex_idx){

	if(complex_idx >= total_truth_table_rna_complex){
		cout << "ERROR: In Reaction::get_concentration, complex_idx is greater than total_truth_table_rna_complex." << endl;
		exit(-1);
	}

	if(rna_complex[complex_idx]->get_total_rna_number() <= 1)
		return (rna_complex[complex_idx]->get_rna(0).get_concentration());

	return (concentration_obj[complex_idx]->get_concentration()/rna_complex[complex_idx]->get_min_rna_concentration());
}


double Reaction::get_total_free_energy(){
	return total_free_energy;
}


double Reaction::get_weight(){
	return weight;
}


int	Reaction::get_mc_steps_count(){
	return mc_steps_count;
}

void Reaction::get_rna_and_complex_vectors(double *rna_energy_vector, double *complex_energy_vector){
	int i;

	for(i = 0; i < main_rna_complex->get_total_rna_number(); i++){
		rna_energy_vector[i] = main_rna_complex->get_rna(i).get_energy();
	}

	for(i = 0; i < total_truth_table_rna_complex + extra_rna_complex; i++){
		complex_energy_vector[i] = rna_complex[i]->get_energy();
	}
}


bool Reaction::get_calculation_mode(){
	return calculation_mode;
}

double Reaction::get_sscore(){
	return sscore;
}


int Reaction::get_calculation_evolution_mode(){
	return calculation_evolution_mode;
}


double Reaction::calculate_complex_mean_pvalue(){
	int 	i;
	double 	pval = 0;

	for(i = 0; i < total_truth_table_rna_complex; i++)
		pval += rna_complex[i]->get_probability();

	return pval;// /total_truth_table_rna_complex;
}


bool Reaction::get_score_weighting_mode(){
	return score_weighting_mode;
}

/********************* Setters **************************/


void Reaction::set_max_iteration(int max_itr){
	max_iteration = max_itr;
}

double Reaction::calculate_energy(){
	return 0;
}

void Reaction::set_input_file(const char *fname){
	strcpy(input_file_name,fname);
}

void Reaction::set_output_file(const char *fname){
	strcpy(output_file_name,fname);
}


void Reaction::set_output_mode(int mode){
	output_mode = mode;
}


void Reaction::set_weight(double new_weight){
	weight = new_weight;
}

void Reaction::set_calculation_mode(bool cmode){
	calculation_mode = cmode;
}

void Reaction::set_sequence_check_mode(int smode){
	for(int i = 0; i < main_rna_complex->get_total_rna_number(); i++)
		main_rna_complex->set_sequence_check_mode(smode);

	for(int i = 0; i < total_truth_table_rna_complex+extra_rna_complex; i++)
		rna_complex[i]->set_sequence_check_mode(smode);
}


void Reaction::set_calculation_evolution_mode(int cmode){
	calculation_evolution_mode = cmode;
}

void Reaction::set_score_weighting_mode(bool mode){
	score_weighting_mode = mode;
}


