/*
 * RNAObject.cpp
 *
 *  Created on: Aug 26, 2010
 *      Author: kamal
 */


#include "RNAObject.h"

using namespace std;

int g_rna_object_count = 0;	// Used to count total number of RNAObject created and to assign unique ID for each RNAObject


#define MZK_OMP_DYNAMIC_LIBRARY_ENV		"DYLD_LIBRARY_PATH"
#define	MZK_OMP_DYNAMIC_LIBRARY_PATH	":./omp_lib:/usr/local/lib"


RNAObject::RNAObject(){
	char env_data[1000];

	base_repeat_limit = 3; //REPEATS;
	distance_limit = DIST;
	temperature = KT;
	stack_height = 0;

	is_mutated = false;
	rbs_distance = 99999;
	total_distance = 99999;

//	id = 0;//g_rna_object_count++;
	energy = 0;
	temperature = KT;
	concentration = 1; // 1 uM

	max_mutation_locations = 1;

	is_rbs_site_unrestricted = false;
	sequence_check_mode = EXTENDED_SEQUENCE_CHECK;
	DNARNACOUNT = 2; // rna37
	DANGLETYPE = 2; // 0:none, 1:some, 2:all

	mfe_structure.validStructs = NULL;
	mfe_structure.nAlloc = 0;
	mfe_structure.nStructs = 0;
	mfe_structure.seqlength = 0;
	mfe_structure.minError = NAD_INFINITY;

	do_symetry_correction = 1;
	do_calculate_pairs = 0;
	subopt_energy_gap = 1; // 1 should give good results
	complexity_num = 3;
	perm_symmetry = 1;
	temp_val = 37+ZERO_C_IN_KELVIN;
	TEMP_K = temp_val;

	probability_value = 0;
	pfunc_z_value = 1;
	concentration = 0;

	save_probability_value = 0;
	save_pfunc_z_value = 1;
	save_concentration = 0;

	rna_structure_mutation_mode = STEM_PRESERVE_ON;

	nu_obj_initialized = false;

	char *ev = getenv(MZK_OMP_DYNAMIC_LIBRARY_ENV);

	if(ev != NULL){
		strcpy(env_data,ev);
		strcat(env_data,MZK_OMP_DYNAMIC_LIBRARY_PATH);
		setenv(MZK_OMP_DYNAMIC_LIBRARY_ENV,env_data,1);
	}
	else{
		strcpy(env_data,MZK_OMP_DYNAMIC_LIBRARY_PATH);
		setenv(MZK_OMP_DYNAMIC_LIBRARY_ENV,env_data,0);
	}

	DNARNACOUNT = 2; // rna37
	DANGLETYPE = 1;// some
	TEMP_K = 37 + ZERO_C_IN_KELVIN;
	LoadEnergies(); // Load nupack energy
	NUPACK_ENERGY_LOADED = 1;

}

RNAObject::RNAObject(string &seqin, string &seq_constraint_in, string &strin, string &target, int id_number){
/*
 * Part of this function is taken from RNAdesigner program
 *
 */
	char 	seq_temp[MAX_CHAR];
	char	str_temp[MAX_CHAR];
	double	energy_gap;
	bool 	fold_possible = false;

	string	seq = seqin;
	string  seq_constr = seq_constraint_in;
	string	str = strin;
	int		seq_len = (int) seq.size();

	RNAObject();

	for(unsigned i = 0; i < seq.length(); i++)
		if(toupper(seq[i]) == 'T')
			seq[i] = 'U';

	for(unsigned i = 0; i < seq_constr.length(); i++)
		if(toupper(seq_constr[i]) == 'T')
			seq_constr[i] = 'U';


	sequence_constraint = seq_constr;
	seq_len = (int) seq_constr.length();

	if(seq[0] == '#'){ // is a line starts with #, no initial sequence is provided, generate randon sequence that conforms with the sequence and structure constraints.
		seq = seq_constr;
	}
	else{
		for(int i = 0; i < seq.length(); i++){
			if(seq_constr[i] != 'N'){
				if(seq_constr[i] != seq[i]){
					cout << "ERROR: seqquence does not match woth the constraint!!" << endl;
					cout << "Check position " << i << " of the sequence and sequence constraint!" << endl;
					exit(-1);
				}
			}
		}
	}

	id = id_number;

	if(seq.length() != str.length() || seq.length() != target.length() || str.length() != target.length()){
		cerr << "ERROR: The length of sequence, structure and target must be equal." << endl;
		exit(-4);
	}

	if(check_structure(str) == false){
		cerr << "ERROR in structure: " << str << endl;
		cerr << "Unbalanced parenthesis." << endl;
		exit(-5);
	}


	int k = 0;

	structure_constraint = str;
	target_structure = target;

	k = 0;
	while((seq[k++] != 'N') && k < (int) seq.length());

	if(k == (int) str.length()){
		sequence = seq;
		structure = str; // only for memory allocation of "structure"! it's overwritten in update()
		update(EVOLVE_FOR_PROBABILITY);
		structure_constraint = structure;

		cout << endl;
		cout << "sequence: " <<  sequence << endl;
		cout << "structure: " << structure << endl;
		cout << "energy: " << energy << endl;
		cout << "probability: " << probability_value << endl;
		return;
	}

	cout << endl;
	for(int i=0; i < 300; i++){
		strcpy(seq_temp, seq.c_str());
		generate_random_sequence(seq_temp, seq_len);
		cout << "Random sequence ("<< i <<"): " << seq_temp << endl;
		energy_gap = (double) inverse_fold(seq_temp, str.c_str());

		string st = seq_temp;
		if(energy_gap < 0.5 || check_valid_sequences(st,-1,BASIC_REPEAT_CHECK)){
			cout << "\nenergy_gap = " << energy_gap << endl;
			fold_possible = true;
			break;
		}
		cout << "Sequence: " << seq_temp << endl << "Structure: " << str << endl;
	}

	for(unsigned idx = 0; idx < strlen(seq_temp); idx++){
		seq_temp[idx] = toupper(seq_temp[idx]);
	}

	sequence = seq_temp;
	structure = str;

	cout << endl;
	cout << "sequence: " <<  sequence << endl;
	cout << "structure: " << structure << endl;

	nu_obj_initialized = false;
	update(EVOLVE_FOR_PROBABILITY); // update energy, structure, and the probability of the of the sequence.
}

RNAObject::~RNAObject(){
	clearDnaStructures(&mfe_structure);
}

void RNAObject::update_sequence(const string &seq){
	sequence = seq;
}

void RNAObject::update_structure(const string &str){
	structure = str;
}

RNAObject& RNAObject::operator=(RNAObject &rna){
	id = rna.id;
	sequence = rna.sequence;
	structure = rna.structure;
	sequence_constraint = rna.sequence_constraint;
	structure_constraint = rna.structure_constraint;
	energy = rna.energy;
	temperature = rna.temperature;
	concentration = rna.concentration;
	base_repeat_limit = rna.base_repeat_limit;
	distance_limit = rna.distance_limit;
	is_mutated = rna.is_mutated;

	return *this;
}

bool RNAObject::check_repeats(const char *seq, int len){
/*
 * This function is taken from RNAdesigner code.
 *
 */

	char 	nt,nt2;
	int 	count;
	int		limit = base_repeat_limit;


	if(len == -1){
		len = 0;
		char *st  = (char *) seq;
		while(*st++)
			len++;
	}


	if(len == -1)
		len = strlen(seq);


	for (int i = 0; i < len-limit; i++){
		nt = toupper(seq[i]);
		if(nt=='N')
			continue;
		count = 0;
		for (int j = 1; j <= limit; j++){
			nt2 = seq[i+j];
			if(nt == toupper(nt2))
				count++;
			else
				break;
		}
		if(count >= limit)
			return false;
	}
	return true;
}


bool RNAObject::check_local_repeats(const string &seqin,int pos){
/*
 * Checks repeats of bases within a region [pos-max_base_repeat_limit, pos+max_base_repeat_limit];
 *
 */

//	const char	*ts = seqin.c_str();
//	return check_repeats(ts,(int) seqin.length());

	return check_valid_sequences(seqin,pos);


	char 		nt;
	char		nt_new;
	int 		count;
	const char	*seq = seqin.c_str();

	int	limit = 3;

	int	len = (int) seqin.size();

	int pos_min;// = pos;
	int	pos_max;// = pos;

	if(pos >= 0){
		pos_min = pos - limit;
		pos_max = pos;

		if(pos_min < 0)
			pos_min = 0;

		if(pos_max + limit >= len-1)
			pos_max = len - limit - 1;
	}else{
		pos_min = 0;
		pos_max = len - limit -1;
	}

	for (int i = pos_min; i < pos_max; i++){
		nt = seq[i];

		if(nt=='N')
			continue;

		count = 0;
		for(int j = 1; j <= limit; j++){
			if(nt == toupper(seq[i+j]))
				count++;
			else
				break;
		}

		if(count >= limit){
			return false;
		}
	}

	return true;
}


bool RNAObject::check_valid_sequences(const string &seqin, int pos, int check_mode){
/*
 * Checks repeats of bases within a region [pos-max_base_repeat_limit, pos+max_base_repeat_limit];
 *
 * The function returns false if any of the following sequences is encountered around position pos
 * the sequence to be avoided: AAAA,CCCC,GGGG,UUUU,KKKKKK, MMMMMM,RRRRRR,SSSSSS,WWWWWW,YYYYYY
 *
 *  Nupack, base wildcards...
 *  N = A C G U
 *  R = A G
 *  Y = C U
 *  M = A C
 *  K = G U
 *  S = C G
 *  W = A U
 *  V = A C G
 *  H = A C U
 *  B = C G U
 *  D = A G U
 *
 */

		char 		nt;
		char 		nt_old;
		const char	*seq = seqin.c_str();
		int			limit = 5;
		int			len = (int) seqin.length();

		int			k6_count = 0;
		int			m6_count = 0;
		int			r6_count = 0;
		int			s6_count = 0;
		int			w6_count = 0;
		int			y6_count = 0;
		int			repeat4_count = 0;

		int 		pos_min;// = pos;
		int			pos_max;// = pos;

		if(pos >= 0){
			pos_min = pos - limit;
			pos_max = pos + limit;

			if(pos_min < 0)
				pos_min = 0;

			if(pos_max > len-1)
				pos_max = len - 1;
		}
		else{	// is pos < 0... check entire sequence
			pos_min = 0;
			pos_max = len - 1;
		}

		nt_old = seq[pos_min];
		for (int i = pos_min+1; i < pos_max; i++){
			nt = seq[i];

			if(nt_old == nt){
				if(++repeat4_count >= 3)
					return false;
			}
			else{
				repeat4_count = 0;
			}

			if(check_mode == EXTENDED_SEQUENCE_CHECK){
				if((nt_old == 'A' || nt_old == 'G') && (nt == 'A' || nt == 'G')){
					if(++r6_count >= 5)
						return false;
				}
				else{
					r6_count = 0;
				}

				if((nt_old == 'C' || nt_old == 'U') && (nt == 'C' || nt == 'U')){
					if(++y6_count >= 5)
						return false;
				}
				else{
					y6_count = 0;
				}

				if((nt_old == 'A' || nt_old == 'C') && (nt == 'A' || nt == 'C')){
					if(++m6_count >= 5)
						return false;
				}
				else{
					m6_count = 0;
				}

				if((nt_old == 'G' || nt_old == 'U') && (nt == 'G' || nt == 'U')){
					if(++k6_count >= 5)
						return false;
				}
				else{
					k6_count = 0;
				}

				if((nt_old == 'C' || nt_old == 'G') && (nt == 'C' || nt == 'G')){
					if(++s6_count >= 5)
						return false;
				}
				else{
					s6_count = 0;
				}

				if((nt_old == 'A' || nt_old == 'U') && (nt == 'A' || nt == 'U')){
					if(++w6_count >= 5)
						return false;
				}
				else{
					w6_count = 0;
				}
			}
			nt_old = nt;
		}
	return true;
}


void RNAObject::calculate_base_fractions(){
	int len = sequence.length();

	for(int i = 0; i < len; i++){
		if(sequence[i] == 'A')
			A_count++;
		else if(sequence[i] == 'U')
			U_count++;
		else if(sequence[i] == 'G')
			G_count++;
		else if(sequence[i] == 'C')
			C_count++;
	}
	total_bases = A_count + U_count + G_count + C_count;
}

void RNAObject::generate_random_sequence(char *seq, int length){
	double	p;

	for(int i=0; i<length; i++){
		if(seq[i] == 'N'){
			p = random_num();
			if(p < 0.25) 		seq[i]='A';
			else if(p < 0.5) 	seq[i]='U';
			else if(p < 0.75) 	seq[i]='G';
			else 				seq[i]='C';
		}
		else{
			if(seq[i]=='A') 	 seq[i]='a';
			else if(seq[i]=='U') seq[i]='u';
			else if(seq[i]=='T') seq[i]='u';
			else if(seq[i]=='G') seq[i]='g';
			else if(seq[i]=='C') seq[i]='c';
		}
	}
}


bool RNAObject::check_structure(const string &str) {
	int unpaired = 0;

	for(int i = 0; i < (int) str.length(); i++) {
		if(str[i] == '(')
			unpaired++;
		else if(str[i] == ')'){
			unpaired--;
			if(unpaired < 0)
				return false;
		}
	}

	if (unpaired == 0)
		return true;

	return false;
}

void RNAObject::string_upper(string &seq){
	for(unsigned k = 0; k < seq.length(); k++)
		seq[k] = toupper(seq[k]);
}


void RNAObject::generate_unique_random_id(int max_num, int *vec_indeces){
	bool 	ri[3*MAX_RNA_NUMBER];
	int		i;
	int		val;

	for(i = 0; i < max_num; i++){
		ri[i] = false;
	}

	i = 0;
	while(true){
		if(i == max_num)
			break;

		val = floor(random_num() * max_num);

		if(ri[val] == false){
			ri[val] = true;
			vec_indeces[i] = val;
			i++;
		}
	}
}

void RNAObject::generate_str_weight(){

}



int RNAObject::find_complementary_site(int pos, const string &strin){
/*
 * Calculates the complementary site position of a given position "pos"
 *
 *
 */
	int count;
	int	comp_count;

	if(strin[pos] == '('){
		count = 1;
		comp_count = 1;
		while(1){
			if(strin[pos+count] == '(')
				comp_count++;
			else if(strin[pos+count] == ')')
				comp_count--;

			if(comp_count == 0)
				return (pos+count);

			count++;
		}
	}
	else if(strin[pos] == ')'){
		count = 0;
		comp_count = 0;
		while(1){
			if(strin[pos-count] == ')')
				comp_count++;
			else if(strin[pos-count] == '(')
				comp_count--;

			if(comp_count == 0)
				return (pos-count);

			count++;
		}
	}

	return -1;
}

int RNAObject::calculate_structure_mismatch(const char *str1, const char *str2, int str_len){
	int len = str_len;
	int	dist_count = 0;

	for(int i = 0; i < len; i++){
		if((str1[i] == '.' && str2[i] != '.') ||
		   (str1[i] != '.' && str2[i] == '.'))
			dist_count++;
	}
	return dist_count;
}


void RNAObject::calculate_structure_weight(){
	int len = (int) structure.length();

	for(int i = 0; i < len; i++){
		str_weight[i] = 0;
		str_coords[i] = -1;
	}

	for(int i = 0; i < len; i++){
		if(str_weight[i] > 0)
			continue;

		if(target_structure[i] != RBS_FREE)
			continue;

		if(structure[i] == SITE_FREE)
			continue;

		if(structure_constraint[i] == ANY_NUCLEOTIDE){
			str_weight[i] = 1;
			str_coords[i] = i;
			int comp_pos = find_complementary_site(i,structure);
			str_weight[comp_pos] = 1;
		}
	}
}


double RNAObject::calculate_energy_per_base(){
	return (energy_per_base = calculate_energy_per_base((char *)structure.c_str(),structure.length(),base_pair_length,energy));
}

double RNAObject::calculate_energy_per_base(char *str,int len, int &bp_length,double en){
	int 	bp_count = 0;

	for(int i = 0; i < len; i++)
		if(str[i] == '(')
			bp_count++;

	bp_length = bp_count;
	return (en/bp_count);
}

void RNAObject::save_parameters(){
	save_is_mutated = is_mutated;
	save_sequence = sequence;
	save_structure = structure;
	save_energy = energy;
	save_energy_per_base = energy_per_base;
	save_total_distance = total_distance;
	save_rbs_distance = rbs_distance;
	save_concentration = concentration;
	save_probability_value = probability_value;
	save_pfunc_z_value = pfunc_z_value;

}

void RNAObject::restore_parameters(){
	is_mutated = save_is_mutated;
	sequence = save_sequence;
	structure = save_structure;
	energy = save_energy;
	energy_per_base = save_energy_per_base;
	total_distance = save_total_distance;
	rbs_distance = save_rbs_distance;
	concentration = save_concentration;
	probability_value = save_probability_value;
	pfunc_z_value = save_pfunc_z_value;
}


void RNAObject::generate_random_word(char *word,int word_len){

	for(int i = 0; i < word_len; i++){
		double p = random_num();

		if(p < 0.25) 	  word[i] = 'G';
		else if(p < 0.50) word[i] = 'C';
		else if(p < 0.75) word[i] = 'U';
		else 			  word[i] = 'A';
	}
}


bool RNAObject::single_mutation(char w, int pos, int check_mode,int mutation_mode){

	char 	wc;
	string	seq;
	int		comp_pos;

	seq = sequence;


	if(pos >= (int) seq.length()){
		cout << "ERROR: pos > sequence_constraint.size... inside single_mutation()\n";
		return false;
	}

	if(seq[pos] == w)
		return false;


	if(sequence_constraint[pos] == ' ' || sequence_constraint[pos] == '+'){
		cout << "Found a space or +....inside single_mutation(\n";
		return false;
	}

	comp_pos = find_complementary_site(pos,structure);

	wc = calc_comp_base(w);

	if(comp_pos < 0 || get_rna_structure_mutation_mode() == STEM_PRESERVE_OFF){ // Free base or mutation without preserving base pairing
		if(sequence_constraint[pos] == 'N'){
			seq[pos] = w;
			if(check_valid_sequences(seq,pos,check_mode)){
				sequence = seq;
				return true;
			}
		}
	}
	else{ // paired base
		if(sequence_constraint[pos] == 'N'){
			if(sequence_constraint[comp_pos] == 'N'){
				seq[pos]= w;
				seq[comp_pos]= wc;
				if(check_valid_sequences(seq,pos,check_mode) && check_valid_sequences(seq,comp_pos,check_mode)){
					sequence = seq;
					return true;
				}
			}
		}
	}

	return false;
}


void RNAObject::update(int calc_evolve_mode){
/*
 * Calculate the MFE structure and energy of the sequence
 * using MultiRNAFold package.
 *
 * if do_fold is false, calculates energy only, given that structure is already calculated;
 * otherwise, calculates both energy and structure.
 *
 */
	DBL_TYPE   zval;

	energy = fold_mfe(sequence,structure,calc_evolve_mode,zval);

	if(calc_evolve_mode == EVOLVE_FOR_STRUCTURE){
		probability_value = 0;
		pfunc_z_value = 1;
	}
	else{
		probability_value = EXP_FUNC(-energy/(Kb*TEMP_K))/zval;
		pfunc_z_value = zval;
	}

	calculate_energy_per_base();

	return;

	if(calc_evolve_mode == EVOLVE_FOR_STRUCTURE){
		if(is_mutated)
			energy = fold_mfe(sequence,structure,calc_evolve_mode,zval);
		pfunc_z_value = 1;
	}
	else{
		energy = fold_mfe(sequence,structure,calc_evolve_mode,zval);
		pfunc_z_value = zval;
		probability_value = EXP_FUNC(-energy/(Kb*TEMP_K))/zval;
	}
	calculate_energy_per_base();
}

double RNAObject::fold_mfe(string &seq, string &str,int calc_evolve_mode,DBL_TYPE &pf_z_value, int calc_one_mfe){
	double 	en;
	int		nseq[MAX_NUPACK_SEQLEN];
	char	char_str[MAX_NUPACK_SEQLEN];

	DANGLETYPE = 2;
	ONLY_ONE_MFE = (calc_one_mfe == 1);
	TEMP_K = 37 + ZERO_C_IN_KELVIN;

	extern int use_cache;

	int dv = use_cache;

//	convertSeq((char *)seq.c_str(),nseq,seq.length());

	if(calc_evolve_mode == EVOLVE_FOR_STRUCTURE){
//		en = mzk_mfeFullWithSym_SubOpt(nseq,seq.length(),&mfe_structure,3,DNARNACOUNT,DANGLETYPE,37,1,subopt_energy_gap,
//				ONLY_ONE_MFE,SODIUM_CONC,MAGNESIUM_CONC,USE_LONG_HELIX_FOR_SALT_CORRECTION);
		en = mfeFullWithSym_SubOpt((char *)seq.c_str(),&mfe_structure,3,DNARNACOUNT,DANGLETYPE,37,1,subopt_energy_gap,
				ONLY_ONE_MFE);

//		DBL_TYPE mfeFullWithSym_SubOpt( char inputSeq[], dnaStructures *mfeStructures,
//				int complexity, int naType,
//				int dangles, DBL_TYPE temperature,
//				int symmetry, DBL_TYPE range, int onlyOne);

		pf_z_value = 1;
	}
	else if(calc_evolve_mode == EVOLVE_FOR_PROBABILITY){
//#pragma omp parallel shared(seqHash,use_cache) private(sizeTerm,pairPr,maxGapIndex,EXTERN_Q,EXTERN_QB)
#pragma omp parallel //shared(seqHash,use_cache) private(sizeTerm,pairPr,maxGapIndex)
		{
#pragma omp sections nowait
			{
#pragma omp section
//				en = mzk_mfeFullWithSym_SubOpt(nseq,seq.length(),&mfe_structure,3,DNARNACOUNT,DANGLETYPE,
//						TEMP_K-ZERO_C_IN_KELVIN,1,subopt_energy_gap,
//						ONLY_ONE_MFE,SODIUM_CONC,MAGNESIUM_CONC,USE_LONG_HELIX_FOR_SALT_CORRECTION);

				en = mfeFullWithSym_SubOpt((char *)seq.c_str(),&mfe_structure,3,DNARNACOUNT,DANGLETYPE,37,1,subopt_energy_gap, ONLY_ONE_MFE);

#pragma omp section
				pf_z_value = pfuncFullWithSym((char *)seq.c_str(),3,DNARNACOUNT,DANGLETYPE, 37, 0, 1);

//				DBL_TYPE pfuncFullWithSym(  char inputSeq[], int complexity, int naType,
//						int dangles, DBL_TYPE temperature,
//						int calcPairs, int permSym);

			}
		}
	}

	make_rna_structure(mfe_structure,seq.c_str(),char_str);
	str = char_str;
	return en;
}

double RNAObject::subsequence_free_probability(string &seq, DBL_TYPE &pf, DBL_TYPE &pf_subopt_value, int &total_str,int &total_str_free, DBL_TYPE subopt_gap, int start_pos,int end_pos, int natype, int dangle_type){
/*
 * Calculates the probability if a subsequence (defined by start_pos and end_pos) of a sequence, seq
 * in a given suboptimal energy gap.
 *
 * by default natype = 2, rna37
 * and dangletype = 2, all dangles
 *
 * Returns
 * subsequence free probability
 *
 * Also returns,

 * pf, the total partition function
 * pf_subopt_value, the suboptimal partition function the sequnece and the partition function in within the energy gap
 */

	int				nseq[MAX_NUPACK_SEQLEN];
	char			char_str[MAX_NUPACK_SEQLEN];
	DBL_TYPE 		zvalue;
	double 			en;
	dnaStructures	ds = {NULL,0,0,0,NAD_INFINITY};

	DANGLETYPE = dangle_type;
	DNARNACOUNT = natype; //2 =  rna37
	TEMP_K = 37+ZERO_C_IN_KELVIN;
	ONLY_ONE_MFE = 0;

//	convertSeq((char *)seq.c_str(),nseq,seq.length());

	zvalue = pfuncFullWithSym((char *)seq.c_str(),3,natype,dangle_type,
			TEMP_K - ZERO_C_IN_KELVIN,0,1);


//	en = mzk_mfeFullWithSym_SubOpt(nseq,seq.length(),&ds,3,natype,dangle_type,
//			TEMP_K-ZERO_C_IN_KELVIN,1,subopt_gap,
//			ONLY_ONE_MFE,SODIUM_CONC,MAGNESIUM_CONC,USE_LONG_HELIX_FOR_SALT_CORRECTION);

	en = mfeFullWithSym_SubOpt((char *)seq.c_str(),&ds,3,natype,dangle_type,
			TEMP_K-ZERO_C_IN_KELVIN,1,subopt_gap,
			ONLY_ONE_MFE);

	make_rna_structure(ds,seq.c_str(),char_str);

	DBL_TYPE pval = 0;
	DBL_TYPE zval_subopt = 0;
	DBL_TYPE en_structure;

	total_str = ds.nStructs;
	int scount = 0;

	bool 	 free_flag;
	for(int i = 0; i < ds.nStructs; i++){
		free_flag = true;
		make_rna_structure(ds,seq.c_str(),char_str,i);

		en_structure = naEnergyPairsOrParensFullWithSym(NULL, char_str, (char *)seq.c_str(), natype, dangle_type, TEMP_K-ZERO_C_IN_KELVIN,1);


		DBL_TYPE ptmp = expl(-en_structure/KT)/zvalue;

		zval_subopt += ptmp;

		for(int k = start_pos-1; k < end_pos; k++){
			if(char_str[k] != '.'){
				free_flag = false;
				break;
			}
		}
		if(free_flag == true){
			pval += ptmp;
			scount++;
		}
	}
	clearDnaStructures(&ds);

	pf = zvalue;
	total_str_free = scount;
	pval = pval/zval_subopt;
	pf_subopt_value = zval_subopt;

	return (double) pval;
}


void RNAObject::calculate_subseq_free_probability(ifstream &infile,string &output_file_name){
	DBL_TYPE 	subopt_gap;
	int			start_pos;
	int			end_pos;
	DBL_TYPE	zval;
	DBL_TYPE	zval_subopt;
	ofstream	outfile;
	int 		total_str,total_str_free;
	string 		seq;
	string		buffer;

	infile >> subopt_gap >> start_pos >> end_pos;

	outfile.open(output_file_name.c_str());

	if(!outfile.is_open()){
		cout << "ERROR: Cannot open output file: " << output_file_name << endl;
		exit(-4);
	}

	streambuf *backup, *fbuffer;

	fbuffer = outfile.rdbuf();
	backup = cout.rdbuf();

	for(int j = 0; j < 2; j++){
		cout << "subopt gap:\t" << subopt_gap << endl;
		cout << "Start pos:\t" << start_pos << endl;
		cout << "End pos:\t" << end_pos << endl;
		cout << "\nNo.\tProb\tZ\tZ_sopt\tTotal\tFree_Str\tSequence\n";


		cout.rdbuf(fbuffer);
	}
	cout.rdbuf(backup);

	int cval = 0;
	bool skip_flag = false;
	double 	subseq_free_prob;
	int		ptmp;

	while(!infile.eof()){
		getline(infile,seq);

		if((ptmp = seq.find("%<")) != string::npos && ptmp == 0) // the block enclosed by SKIP_START SKIP_END will be skipped
			skip_flag = true;

		if(skip_flag == true){
			if((ptmp = seq.find("%>")) != string::npos && ptmp == 0){
				skip_flag = false;
			}
			continue;
		}

		if(seq[0] == '%' || seq.length() == 0){
			for(int j = 0; j < 2; j++){
				cout << seq << endl;
				cout.rdbuf(fbuffer);
			}
			cout.rdbuf(backup);
			continue;
		}

		if(seq.find("Sequence:") != string::npos){
			int pval = seq.find_first_of('\t');
			seq = seq.substr(pval+1,seq.length());
		}else if(seq.find_first_of(":=.()") != string::npos){
			continue;
		}

		subseq_free_prob = subsequence_free_probability(seq, zval, zval_subopt, total_str, total_str_free, subopt_gap, start_pos, end_pos);

		cval++;

		for(int j = 0; j < 2; j++){
			cout << cval << ":\t" << setprecision(2) << subseq_free_prob << '\t' << zval << '\t' << zval_subopt << '\t' << total_str << '\t' << total_str_free << '\t' << seq << endl;
			cout.rdbuf(fbuffer);
		}
		cout.rdbuf(backup);
	}

	cout << "\nCalculation finished!" << endl;
	cout << "The results are saved in " << output_file_name << endl;
	outfile.close();
}


DBL_TYPE RNAObject::calculate_pfunc(string &seq){
	int nseq[MAX_NUPACK_SEQLEN];

//	convertSeq((char *)seq.c_str(),nseq,seq.length());
//	return nu_obj->pfunc(nseq);;
	return pfuncFullWithSym((char *)seq.c_str(),3,DNARNACOUNT,DANGLETYPE,
			37,0,1);
}

DBL_TYPE RNAObject::calculate_pfunc(){
	int nseq[MAX_NUPACK_SEQLEN];

//	convertSeq((char *)sequence.c_str(),nseq,sequence.length());
	pfunc_z_value = pfuncFullWithSym((char *)sequence.c_str(),3,DNARNACOUNT,DANGLETYPE,
			37,0,1);
	return pfunc_z_value;
}

double RNAObject::calculate_probability(DBL_TYPE z_value){
	int nseq[MAX_NUPACK_SEQLEN];

	if(z_value > 0){
		pfunc_z_value = z_value;
		probability_value = EXP_FUNC(-energy/KT)/z_value;
		return probability_value;
	}

//	convertSeq((char *)sequence.c_str(),nseq,sequence.length());
	pfunc_z_value = pfuncFullWithSym((char *)sequence.c_str(),3,DNARNACOUNT,DANGLETYPE,
			37,0,0);
	probability_value = expl(-energy/KT)/pfunc_z_value;
	return probability_value;
}



void RNAObject::make_rna_structure(const dnaStructures &ds,const char *seq,char *str,int index){
	  int 	pos = 0;

	  for(int j = 0; j < ds.seqlength; j++) {
		  if(seq[pos] == '+' || seq[pos] == ' '){
			  str[pos++] = '+';
		  }

		  if((ds.validStructs)[index].theStruct[j] > j){
			  str[pos++]  = '(';
		  }
		  else if(  (ds.validStructs)[index].theStruct[j] == -1 ){
			  str[pos++] = '.';
		  }
		  else{
			  str[pos++] = ')';
		  }
	  }
	  str[pos] = 0;
}


char **RNAObject::make_multiple_rna_structures(const dnaStructures &ds,const char *seq,int &number_of_structures){
	  char	**strout;

	  strout = (char **)malloc(sizeof(char*)*ds.nStructs);
	  number_of_structures = ds.nStructs;

	  cout << "All the closely spaced structures..." << endl;

	  for(int count = 0; count < ds.nStructs; count++){
		  strout[count] = (char *)malloc(ds.seqlength*sizeof(char));
		  make_rna_structure(ds,seq,strout[count],count);
		  cout << strout[count] << endl;

	  }
	  return strout;
}



// ************** Getters *********************//


int RNAObject::get_id(){
	return id;
}

string &RNAObject::get_rna_name(){
	return rna_name;
}

int RNAObject::get_sequence_length(){
	return (int) sequence.length();
}

int RNAObject::get_base_pair_length(){
	return base_pair_length;
}

string& RNAObject::get_sequence(){
	return sequence;
}

string& RNAObject::get_structure(){
	return structure;
}

string& RNAObject::get_structure_constraint(){
	return structure_constraint;
}

string& RNAObject::get_sequence_constraint(){
	return sequence_constraint;
}

string& RNAObject::get_target_structure(){
	return target_structure;
}

double RNAObject::get_energy(){
	return energy;
}

double RNAObject::get_energy_per_base(){
	return energy_per_base;
}

double RNAObject::get_save_energy(){
	return save_energy;
}
double RNAObject::get_concentration(){
	return concentration;
}

double RNAObject::get_temperature(){
	return temperature;
}

bool RNAObject::get_mutation_flag(){
	return is_mutated;
}

int RNAObject::get_rbs_distance(){
	return rbs_distance;
}

int RNAObject::get_total_distance(){
	return total_distance;
}

int RNAObject::get_distance_limit(){
	return distance_limit;
}

int RNAObject::get_max_mutation_locations(){
	return max_mutation_locations;
}

int	RNAObject::get_repeat_limit(){
	return base_repeat_limit;
}


double RNAObject::get_GC_fraction(){
	return (double) (G_count+C_count)/total_bases;
}

double RNAObject::get_AU_fraction(){
	return (double) (A_count+U_count)/total_bases;
}

void RNAObject::get_base_content_fraction(double &A_frac,double &U_frac,double &G_frac,double &C_frac){
	int 	len = total_bases;

	A_frac = (double) A_count/len;
	U_frac = (double) U_count/len;
	G_frac = (double) G_count/len;
	C_frac = (double) C_count/len;
}
//
//NUPACKObject& RNAObject::get_nupack_object(){
//	return *nu_obj;
//}



int RNAObject::get_sequence_check_mode(){
	return sequence_check_mode;
}

void RNAObject::set_sequence_check_mode(int cmode){
	sequence_check_mode = cmode;
}


double RNAObject::get_probability(){
	return probability_value;
}

DBL_TYPE RNAObject::get_pfunc_z_value(){
	return pfunc_z_value;
}

int RNAObject::get_evolve_mode(){
	return evolve_mode;
}


void RNAObject::get_total_free_distance(int &distance,int &free_site_count){
	int len = (int) sequence.length();
	int dist_count = 0;
	int total_site_count = 0;

	for(int i = 0; i < len; i++){
		if(target_structure[i] == RBS_FREE){
			total_site_count++;
			if(structure[i] != '.')
				dist_count++;
		}
	}

	free_site_count = total_site_count;
	distance = dist_count;
}


void RNAObject::get_total_bound_distance(int &distance,int &bound_site_count){
	int len = (int) sequence.length();
	int dist_count = 0;
	int total_site_count = 0;

	for(int i = 0; i < len; i++){
		if(target_structure[i] == RBS_BOUND){
			total_site_count++;
			if(structure[i] == '.')
				dist_count++;
		}
	}

	bound_site_count = total_site_count;
	distance = dist_count;
}

double RNAObject::getcputime()
{
	double t;
	struct timeval tim;
	struct rusage ru;

	getrusage(RUSAGE_SELF, &ru);
	tim = ru.ru_stime;
	t = (double)tim.tv_sec * 1000000.0 + (double)tim.tv_usec;
	return t;
}


int	RNAObject::get_subopt_energy_gap(){
	return subopt_energy_gap;
}


char RNAObject::calc_non_comp_base(char nt){
/*
 * Allow woble pairing with probability given in WOBBLE_PROBABILITY.
 *
 *
 */
	if(nt == 'G') return random_num() >= WOBBLE_PROBABILITY? 'G':'U';
	if(nt == 'U') return random_num() >= WOBBLE_PROBABILITY? 'U':'G';

	return nt;
}

char RNAObject::calc_comp_base(char nt){

	if(nt == 'A') return 'U';
	if(nt == 'C') return 'G';
	if(nt == 'U') return (random_num() >= WOBBLE_PROBABILITY? 'A':'G');
	if(nt == 'G') return (random_num() >= WOBBLE_PROBABILITY? 'C':'U');

	return nt;
}

int RNAObject::	get_rna_structure_mutation_mode(){
	return rna_structure_mutation_mode;
}


// ************** Setters *********************//

void RNAObject::set_id(int id_value){
	id = id_value;
}

void RNAObject::set_rna_name(string &name){
	rna_name = name;
}

void RNAObject::set_energy(double energy_value){
	energy = energy_value;
}

void RNAObject::set_energy_per_base(double energy_pb_value){
	energy_per_base = energy_pb_value;
}

void RNAObject::set_mutation_flag(bool flag){
	is_mutated = flag;
}

void RNAObject::set_sequence(const string &seq){
	sequence = seq;
}

void RNAObject::set_structure(const string &str){
	structure = str;
}

void RNAObject::set_sequence_constraint(const string &seq){
	sequence_constraint = seq;
}

void RNAObject::set_structure_constraint(const string &str){
	structure_constraint = str;
}

void RNAObject::set_target_structure(const string &str){
	target_structure = str;
}

void RNAObject::set_base_pair_length(int bp_length){
	base_pair_length = bp_length;
}

void RNAObject::set_rbs_distance(int rbs_dist){
	rbs_distance = rbs_dist;
}

void RNAObject::set_total_distance(int total_dist){
	total_distance = total_dist;
}


void RNAObject::set_distance_limit(int dist_limit){
	distance_limit = dist_limit;
}


void RNAObject::set_concentration(double conc_value){
	concentration = conc_value;
}

void RNAObject::set_repeat_limit(int rlimit){
	base_repeat_limit = rlimit;
}


void RNAObject::set_max_mutation_locations(int mutation_lications){
	max_mutation_locations = mutation_lications;
}


void RNAObject::set_probability(double prob_value){
	probability_value = prob_value;
}

void RNAObject::set_pfunc_z_value(DBL_TYPE zvalue){
	pfunc_z_value = zvalue;
}


void RNAObject::set_evolve_mode(int e_mode){
	evolve_mode = e_mode;
}

void RNAObject::set_subopt_energy_gap(int egap){
	subopt_energy_gap = egap;
}


void RNAObject::set_rna_structure_mutation_mode(int mode){
	rna_structure_mutation_mode = mode;
}

