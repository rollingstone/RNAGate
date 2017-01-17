/*
 * ReadInputFile.h
 *
 *  Created on: Nov 8, 2010
 *      Author: kamal
 */

#ifndef READINPUTFILE_H_
#define READINPUTFILE_H_

#include <iostream>
#include <fstream>
#include <string>


#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

enum RNA_INPUT_FILE_IDS{
	BX = 0,//0xCAFE, // :)
	UNKNOWN_ID,
	NUMBER_ID,
	UTR_ID,
	YES_ID,
	NO_ID,
	MAX_ITER_ID,
	RNA_DATA_HEADER_ID,
	ID_ID,
	SEQUENCE_ID,
	STRUCTURE_ID,
	TARGET_ID,
	GATE_TYPE_ID,
	GATE_NAME_ID,
	ROW_ID,
	COL_ID,
	BEGIN_ID,
	END_ID,
	NAME_ID,
	CONCENTRATION_ID,
	MUTATION_WORD_SIZE_ID,
	MUTATION_LOCATIONS_ID,
	DISTANCE_LIMIT_ID,
	REPEAT_LIMIT_ID,
	TRUTH_TABLE_DATA_HEADER_ID,
	PREDEFINED_TRUTH_TABLE_ID,
	TABLE_SIZE_ID
};

//#define		NUMBER_ID		1+BX
//#define		UTR_ID			2+BX
//#define		YES_ID			3+BX
//#define		NO_ID			4+BX
//#define		MAX_ITER_ID		5+BX
//#define		ID_ID			6+BX
//#define		SEQUENCE_ID		7+BX
//#define		STRUCTURE_ID	8+BX
//#define		TARGET_ID		9+BX
//#define		GATE_NAME_ID	10+BX
//#define		ROW_ID			11+BX
//#define		COL_ID			12+BX
//#define		BEGIN_ID		13+BX
//#define		END_ID			14+BX
//#define		NAME_ID			15+BX
//#define		CONCENTRATION_ID 16+BX
//#define		MUTATION_WORD_SIZE_ID	17+BX
//#define		MUTATION_LOCATIONS_ID	18+BX
//#define		DISTANCE_LIMIT_ID		19+BX
//#define		REPEAT_LIMIT_ID			20+BX
//#define		UNKNOWN_ID		-1




string 	remove_comment(string str);
void  	kill_leading_and_trailing_spaces_v0(string &str);
string& kill_leading_and_trailing_spaces(string &str);
string& kill_all_spaces(string &str);
bool 	compare_string(string &str, char *char_str,bool is_case_sensitive = false);
bool 	compare_string_obj(string &str1, string &str2, bool is_case_sensitive = false);
bool 	compare_string_char(const char *str1, const char *str2, bool is_case_sensitive = false);
bool 	compare_string_cci(const char *str1,const char *str2);

//bool 	compare_string(string &str1, string &str2);
void 	get_numbers_from_string(string &str, int *val, int &n);
bool 	split_string(string &str, string &str_left, string &str_right);
bool 	read_line(ifstream &file, string &str);
void 	scanline(ifstream &file,string &str);
void 	read_numbers_from_string(string &str,int *num, int &num_count);
bool 	make_even_spaced_string(string &str);
int 	get_param_id(string &str);
int 	read_params(string &str, string &str_value, int &value,string &strout);

bool 	find_char(char * str, char ch);


#ifdef __cplusplus
}
#endif
#endif /* READINPUTFILE_H_ */
