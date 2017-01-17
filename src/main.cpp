//============================================================================
// Name        : RNAgate.cpp
// Author      : Marzuk M Kamal
// Version     :
// Copyright   : (c) Marzuk M Kamal, 2010
// Description : RNAgate v1.1
//============================================================================

/*
 * RNAGate v 1.1
 *
 * Usage:
 * rnagate inputfilename.txt outputfilename.txt -OPTIONS
 *
 * OPTIONS must be either or both of the following:
 *
 * f = write to file outputfilename.txt
 * s = write to standard output, i.e. display the results on screen
 *
 * r = calculate the subsequence free probability....
 *
 * To write both to file and display results on screen
 *
   rnagate inputfilename.txt outputfilename.txt -sf

 * command line for subsequence free probability calculation
 *
   rnagate inputfilename.txt outputfilename.txt -r
 *
 *
 * Note that the input file format is different for -r option. See below.
 *
 */

/*********************************************************************************
 * In order to compile the code, please use Eclipse C++ development tools
 *
 * You need to install gcc 4.5.0+ to compile and run the code
 *
 * I used gcc 4.6.0 release verison to compile the code. it's required for openmp library
 *
 *********************************************************************************/


/*
 * input file format for OPTION = -s or -f or -sf
 %
 % comments
 % 5'UTR sequence must be the first one
 %
 % Please follow exact gaps between paragraphs
 % All inputs must start from the first column
 %

 N (number of RNAs involved including UTR)
 .....((((....))).... (monomer structure constraint)
 AGCAGCGCAGGCTUAGAAUG (initial sequence, if the initial sequnece is not present, put only a # in this line)
 NNNNNNNNNNGCTUNNNAUG (monomer sequence constraint, the constraint bases must match with the initial sequence (if provided) in the above line)
 @@@@@@@@@@RRRRRRRRRR (constraint when the RNAs make complexes, R means RBS or any associated base that affects translation)
 X			(number of mutations sites for the UTR)
 D			(maximum allowed distance of this monomer)
 Y			(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken)

 .....((((....))).... (monomer structure constraint for sRNA1)
 #
 NNNNNNNNNNNNNNNNNNNN (monomer sequence constraint)
 @@@@@@@@@@@@@@@@@@@@ (constraint when the RNAs make complexes)
 X			(number of mutations sites for the RNA)
 D			(maximum allowed distance of this monomer)
 Y			(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken)

 .....((((....))).... (monomer structure constraint sRNA2)
 #
 NNNNNNNNNNNNNNNNNNNN (monomer sequence constraint)
 @@@@@@@@@@@@@@@@@@@@ (constraint when the RNAs make complexes)
 X			(number of mutations sites for the sRNA2)
 D			(maximum allowed distance of this monomer)
 Y			(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken)

 AND	(name of the predefined logic gate that we would like to achieve, the options are YES,NOT,AND, OR,XOR,NAND,XNOR,NOR)
 BASIC_REAPEAT_CHECK (or EXTENTED_SEQUENCE_CHECK)  (Please see the documentation of check_valid_sequence)
 MONTE_CARLO_METHOD  (or SIMULATED_ANNEALING_METHOD)  (choose the method of calculation)
  */


/*
 * For calculating RBSfree probability:
 *
 * command line:
 *
 * rnagate free_in.txt rfree_out.txt -r
 *
 * the input file format.
 *
 * rfee.txt
 *
2 % kcal/mol, sub-optimal energy gap, comments are allowed in the input file.
4 10 % start and end positions of the subsequence
% each line of the input file contains an RNA or RNA-complex
%
% RNA sequence one
ACGUGUAGUGAUAGUGAGA
% another RNA-complex
ACGUGUAGUGAUAGUGAGA+GAUGAUAGUAGAUAGCAGC
ACGUGUAGUGAUAGUGAGA+GAUGAUAGUAGAUAGCAGC+CGCGGCGAGAGAG

*/


/*
 * Before running RNAGate, the environmental variable RNAGATEHOME, must be set to the home folder of RNAGate which contains rnagate binary
 * and the "parameters" folder
 *
 * Path tree:
 *
 * ~/
 *  RNAGate/
 *  		parameters/
 *
 *
 *
 * command:
 *
 * export RNAGATEHOME=~/RNAGate
 * export PATH=$PATH:$RNAGATEHOME
 *
 * if not set, the default value of RNAGATEHOME will be
 *
 * RNAGATEHOME=~/RNAGate
 *
 *
 */


/*
 *
 * All the function and variable names are self-explanatory, in case of any question,
 * please write to marzuk9999@gmail.com
 *
 */


#include "RNAObject.h"
#include "RNAComplex.h"
#include "Reaction.h"
#include "random.h"


using namespace std;

extern double DG_PARAM;

int main(int argc, char **argv) {
	double temperature = 37;
	time_t start_time, end_time;
	string input_name;
	string output_name;
	int *selected_table = NULL;
	int calc_mode;
	int seq_check_mode = EXTENDED_SEQUENCE_CHECK;

	ifstream infile;
	string buffer;
	int val;
	int output_mode;

	string seq_val;
	string seq_constraint_val;
	string str_val;
	string target_val;

	ARGC_VALUE = argc;
	ARGV_VALUE = argv;

	SCORING_FUNCTION_MODE = 1; // 1 for Standard free energy one, 0 for rnades scoring function
	RNA_PROB_CALC_ENFORCE_FLAG = 1;	// force probability improvement in RNA structure

//****************** set up environmental variables if not set already **************

	string rnagate_path =  getenv("HOME") + (string)"/RNAGate";
	string export_path =  (string)"export PATH=$PATH:" + rnagate_path;

	setenv("RNAGATEHOME",rnagate_path.c_str(),0);
	system(export_path.c_str());

//***********************************************************************************

	init_random(); // initialize random number generator

	cout << "\nRNAGate v 1.0\n(c) 2010-2011 Marzuk M Kamal\nInstitute of Systems and Synthetic Biology, Genopole-CNRS, Evry, France";


	if (argc < 3) {
		cout    << "Usage: rnagate inputfilename.txt outputfilename.txt -OPTIONS\n"
				<< endl;
		cout << "OPTIONS:" << endl;
		cout << "-f: write output to file." << endl;
		cout << "-s: write output to the standard output." << endl;
		cout << "-sf: write output to both file and standard output." << endl << endl;
		cout << "-r: calculate probability of a subsequence of an RNA/RNA-complex to be free within a given energy gap." << endl;

		cout << "\nIMPORTANT:The folder, \"parameters\" must be located at " << getenv("RNAGATEHOME") << ",\n\t or in a folder defined by RNAGATEHOME"
			 << endl;

		exit(-1);
	}

	input_name = argv[1];
	output_name = argv[2];

	output_mode = 0;

	int 	param_count = 4;

	if (argc >= param_count) {
		if (argv[param_count-1][0] != '-') {
			cout << "OPTION must be preceded by -" << endl;
			exit(-1);
		}
		if(find_char(argv[param_count-1], 'f'))
			output_mode |= WRITE_TO_FILE;
		if(find_char(argv[param_count-1], 's'))
			output_mode |= WRITE_TO_STDOUT;

		if(find_char(argv[param_count-1],'r')){
			output_mode = CALCULATE_SUBSEQUENCE_FREE_PROBABILITY;
		}

	} else
		output_mode = WRITE_TO_FILE; // write to file only


	if(argc >= ++param_count){
		int v = atoi(argv[param_count-1]);

		if(v >= 0)
			RNA_PROB_CALC_ENFORCE_FLAG = v;
	}

	if(argc >= ++param_count){
		double v =  atof(argv[param_count-1]);
		if(v >= 0)
			DG_PARAM = v;
	}

	if(argc >= ++param_count)
		SCORING_FUNCTION_MODE = atoi(argv[param_count-1]);

	cout << "Input file name:\t"  << input_name << endl
		 << "Output file name:\t" << output_name << endl;


	infile.open(input_name.c_str());

	if (!infile.is_open()) {
		cout << "ERROR: cannot open file: " << input_name.c_str() << endl;
		exit(-1);
	}

	if(output_mode & CALCULATE_SUBSEQUENCE_FREE_PROBABILITY){
			RNAObject::calculate_subseq_free_probability(infile,output_name);
			infile.close();
			return 0;
	}


	time(&start_time);

	RNAObject 	*rna_object;
	RNAComplex 	rnacomplex_object;
	Reaction 	reaction_object;

	reaction_object.set_score_weighting_mode(ADAPTIVE_SCORE_WEIGHTING);


	while (!infile.eof()) {
		scanline(infile, buffer);
		cout << buffer << endl;
		if (!(buffer[0] == '%' || buffer[0] == ' ' || buffer.length() == 0))
			break;
	}

	val = atoi(buffer.c_str());

	int idx = 0;

	cout << "\nReading input file " << input_name << endl;

	int	rna_count = 0;

	while(1){
		scanline(infile, buffer);

		if (buffer.length() == 0) {
			cout << endl;
			continue;
		}

		if (idx++ >= val)
			break;

		str_val = buffer;
		scanline(infile, seq_val);
		scanline(infile, seq_constraint_val);
		scanline(infile, target_val);
		scanline(infile, buffer);
		int mut_sites = atoi(buffer.c_str());
		scanline(infile, buffer);
		int dist_param = atoi(buffer.c_str());
		scanline(infile, buffer);
		int preserve_stem = atoi(buffer.c_str());

		cout << str_val << endl << seq_constraint_val << endl << seq_val << endl << target_val << endl;


		rna_object = new RNAObject(seq_val, seq_constraint_val, str_val, target_val, rna_count++);
		rna_object->set_distance_limit(dist_param); // set distance between the
		rna_object->set_max_mutation_locations(mut_sites);
		rna_object->set_concentration(1.0); // will be changed
		rna_object->set_rna_structure_mutation_mode(preserve_stem);
		rnacomplex_object.add_rna(*rna_object);
	}

	while(buffer.length() == 0 && infile.eof() == false)
		scanline(infile, buffer);

	cout << buffer << endl;

	if (buffer == "YES") {
		if (val != 2) {
			cout
					<< "ERROR: only two RNAs must be specified for a YES and NOT gate."
					<< endl;
			exit(-1);
		}
		selected_table = truth_table_YES_gate;
	} else if (buffer == "NOT") {
		if (val != 2) {
			cout
					<< "ERROR: only two RNAs must be specified for a YES and NOT gate."
					<< endl;
			exit(-1);
		}
		selected_table = truth_table_NOT_gate;
	} else if (buffer == "AND") {
		if (val != 3) {
			cout
					<< "ERROR: only three RNAs must be specified for a AND, OR and XOR gates."
					<< endl;
			exit(-1);
		}
		selected_table = truth_table_AND_gate;
	} else if (buffer == "OR") {
		if (val != 3) {
			cout
					<< "ERROR: only three RNAs must be specified for a AND, OR and XOR gates."
					<< endl;
			exit(-1);
		}
		selected_table = truth_table_OR_gate;
	} else if (buffer == "XOR") {
		if (val != 3) {
			cout    << "ERROR: only three RNAs must be specified for a AND, OR and XOR gates."
					<< endl;
			exit(-1);
		}
		selected_table = truth_table_XOR_gate;
	}

	scanline(infile, buffer);
	int max_iter = atoi(buffer.c_str());
	scanline(infile, buffer);
	if (compare_string_cci(buffer.c_str(), "BASIC_REPEAT_CHECK"))
		seq_check_mode = BASIC_REPEAT_CHECK;
	else
		seq_check_mode = EXTENDED_SEQUENCE_CHECK;

	scanline(infile, buffer);
	if (compare_string_cci(buffer.c_str(), "MONTE_CARLO_MODE"))
		calc_mode = MONTE_CARLO_MODE;
	else
		calc_mode = SIMULATED_ANNEALING_MODE;

	cout << "Maximum iteration:\t" << max_iter << endl;

	reaction_object.create_rna_complex_from_truth_table(&rnacomplex_object, selected_table);

	reaction_object.set_max_iteration(max_iter); // set maximum iteration
	reaction_object.set_output_mode(output_mode);
	reaction_object.set_sequence_check_mode(seq_check_mode);
	reaction_object.set_calculation_mode(calc_mode);
	reaction_object.set_output_file(output_name.c_str());
	reaction_object.run();

	time(&end_time);
	double dt = difftime(end_time, start_time) / 60.0;
	cout << "Total time of simulation = " << dt << "minutes." << endl;
	cout << "Phew! that was hard work!" << endl;

	double	cpu_time  = RNAObject::getcputime(); // not verified

	cout << "Total CPU time taken " << cpu_time << endl;

	delete_random();
	return 0;
}

