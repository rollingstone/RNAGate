/*
 * read_write.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: kamal
 */

#include "ReadWrite.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define		CONFIGURATION_HEADER		"[RNA_Configuration]"
#define		TRUTH_TABLE_HEADER			"[Truth_Table]"
#define		SEQUENCE_LABEL				"Sequence"
#define		SEQUENCE_CONSTARINT_LABEL	"Sequence_Constraint"
#define		STRUCTURE_LABEL				"Structure"
#define		TARGET_STRUCTURE_LABEL		"Target_Structure"
#define		ROW_LABEL					"Row"
#define		COLUMN_LABEL				"Column"
#define		MAX_ITERATION_LABEL			"Maximum_Iteration"


#define		TMP_STRLEN		2000


char allowed_chars[150] = "1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz±!@#$%^&*()_+=-¤~`,.:;'\"<>?|{}[]";

string remove_comment(string str){
	int 	pos = 0;
	int		len = (int) str.length();
	string	tmp_str;

	if(len == 0){
		return str;
	}

	while(pos < len && str[pos++] != '%');

	if(pos == len)
		return str;

	tmp_str = str.substr(0,pos-1);
//	str = tmp_str;

	return tmp_str;
}

string remove_init_spaces(string str){
	unsigned 	pos = 0;
	string		str1 = str;

	while(pos++ < str1.length())
		if(str1[pos] != ' ')
			break;

	str1 = str1.substr(pos,str1.length()-1);
	return str1;
}

string remove_end_spaces(string str){
	string			str1 = str;
	unsigned int 	pos = str1.length()-1;

	while(pos--)
		if(str1[pos] != ' ')
			break;

	str1 = str1.substr(0,pos);
	return str1;
}

void make_single_spaced_string(string &str){
	int		len;
	int 	i;
	int		idx = 0;
	string	sv;

	kill_leading_and_trailing_spaces_v0(str);

	len = (int) str.length();

	i = 0;
	while(i < len){
		if(str[i] == ' '){
			sv[idx++] = str[i];
			while(i < len && str[i++] == ' ');
		}else{
			sv[idx++] = str[i++];
		}
	}
	str = sv;
}


void get_numbers_from_string(string &str, int *val, int &n){

	int 	len;
	string 	sv = str;
	int 	i;
	int 	sc = 0;
	int 	idx = 0;
	bool	sflag = false;

	make_single_spaced_string(str);

	len = (int) str.length();
	int 	i_old = 0;
	for(i = 0; i < len; i++){
		if(str[i] == ' ' || i == len-1){
			sv = str.substr(i_old, i-i_old+1);
			val[idx++] = atoi(sv.c_str());
			i_old = i+1;
		}
	}
	n = idx;
}


void kill_leading_and_trailing_spaces_v0(string &str){
	int 			pos = 0;
	int 			len = (int) str.length();
	string			tmp_str;

	if(!str.size())
		return;

	if(str[0] == ' '){
		while(pos < len && str[pos++] == ' ');
		tmp_str = str.substr(pos-1,str.length()-pos+1);
		str = tmp_str;
	}

	len = str.length();
	pos = len-1;

	if(str[pos] == ' '){
		while(pos >=0 && str[pos--] == ' ');
		tmp_str = str.substr(0,pos+2);
		str = tmp_str;
	}
}


string& kill_leading_and_trailing_spaces(string &str){
	int 			pos = 0;
	int 			len = (int) str.length();
	string			tmp_str = str;
	string			str1;

//	cout << "str = " << str << endl;
//	cout << "TMP_STR  = " << tmp_str << endl;

	if(!str.size())
		return str;

	if(str[0] == ' '){
		while(pos < len && str[pos++] == ' ');
		tmp_str = str.substr(pos-1);
	}

//	cout << "After space removal Left =\t" << tmp_str << endl;

	str1 = tmp_str;
	len = str1.length();
	pos = len-1;

	if(str1[pos] == ' '){
		while(pos >=0 && str1[pos--] == ' ');
		tmp_str = str1.substr(0,pos+2);
	}

//	cout << "After space removal Right =\t" << tmp_str << endl;

//	str = tmp_str;
	return tmp_str;
}

void lower_string(string &str){
	int	pos = 0;
	while(pos++ < str.length())
		str[pos] = tolower(str[pos]);
}


bool split_string(string &str, string &str_left, string &str_right){
	size_t 	pos;
	string 	stemp;

	string 	str1 = str;
	string 	str2;

	pos = str.find('=');

	if(pos == string::npos){
		if(str.length() == 0)
			return false;

		str_left = str1;
		str_right = "";
		return true;
	}

	stemp = str.substr(0,pos-1);
	str1 = kill_leading_and_trailing_spaces(stemp);
	stemp = str.substr(pos+1,str.length()-1);
	str2 = kill_leading_and_trailing_spaces(stemp);

	str_left  = str1;
	str_right = str2;

//	cout << "left: " << str_left << endl << "right: " << str_right << endl;

	return true;
}

bool compare_string_cci(const char *str1,const char *str2){

	if(strlen(str1) != strlen(str2))
		return false;

	while(*str1){
		if(toupper(*str1++) != toupper(*str2++))
			return false;
	}

	return true;
}


void scanline(ifstream &file,string &str){
	string 	buffer;
	getline(file,buffer);
	if(buffer.length() == 0){
		str.clear();
		return;
	}

	unsigned pos = buffer.find_first_of('%');
	if(pos != string::npos){
		if(pos == 0){
			str.clear();
			return;
		}
		buffer = buffer.substr(0,pos);
	}

	pos = buffer.find_first_of(allowed_chars);
	if(pos != string::npos){
		buffer = buffer.substr(pos,buffer.length()-pos);
	}

	pos = buffer.find_last_of(allowed_chars);
	if(pos != string::npos){
		buffer = buffer.substr(0,pos+1);
	}

	str = buffer;
}

void read_numbers_from_string(string &str,int *num, int &num_count){
	char 	*strbuf = (char *)str.c_str();
	char 	*cend;
	int 	count = 0;

	cend = &strbuf[strlen(strbuf)];

	while(strbuf != cend){
		num[count++] = (int) strtol(strbuf,&strbuf,10);
	}
	num_count = count;
}


bool find_char(char * str, char ch){
	int len = strlen(str);
	int	i = 0;

	while(i++ < len){
		if(toupper(str[i]) == toupper(ch))
			return true;
	}

	return false;
}


bool read_line(ifstream &file, string &str){
	string	str_buffer;

	if(file.eof())
		return false;

	getline(file,str_buffer);

//	cout << "STR_buffer = " << str_buffer << endl;

	str_buffer = remove_comment(str_buffer);

//	cout << "string read = " << str << endl;

	if(!str_buffer.length())
		return false;

	str_buffer = kill_leading_and_trailing_spaces(str_buffer);

//	str_buffer = kill_all_spaces(str_buffer);

//	cout << "after killing space, str = " << str << endl;


	if(str_buffer.length() == 0)
		return false;

	str = str_buffer;
	return true;
}

bool make_even_spaced_string(string &str){
	string	newstr;
	int		len = (int) newstr.length();
	int 	pos = 0;

	int i = 0;
	while(i < str.length()){
		if(str[i] == ' '){
			newstr.append(str.substr(pos,i));
			int k =0;

			while(str[i + ++k] == ' ');

			pos = i + k +1;
			i = pos;
		}
		else{
			i++;
		}
	}

//	str

	return true;
}


string& kill_all_spaces(string &str){
	string	str_buffer;
	int		len = (int) str.length();
	int		idx = 0;

	for(int i = 0; i < len; i++){
		if(str[i] != ' ')
			str_buffer[idx++] = str[i];
	}

	return str_buffer;
}

bool compare_string(string &str1, char *str2, bool is_case_sensitive){
	return compare_string_char(str1.c_str(), str2, is_case_sensitive);
}

bool compare_string_obj(string &str1, string &str2, bool is_case_sensitive){
	return compare_string_char(str1.c_str(), str2.c_str(), is_case_sensitive);
}

bool compare_string_char(const char *str1, const char *str2, bool is_case_sensitive){
	int 	len = strlen(str1);
	int 	i = 0;
	char	str1_buffer[TMP_STRLEN];
	char	str2_buffer[TMP_STRLEN];

	strcpy(str1_buffer,str1);
	strcpy(str2_buffer,str2);

	if(len != (int) strlen(str2_buffer)){
//		cout << str1 << " and " << str2 << " are of different length."<< endl;
		return false;
	}

	while(i++ < len){
		if(is_case_sensitive == false){
			if(toupper(str1_buffer[i]) != toupper(str2_buffer[i]))
				return false;
		}
		else if(is_case_sensitive == true){
			if(str1[i] != str2[i])
				return false;
		}
	}

	return true;
}


int get_param_id(string &str){
/*
 *
 * Returns an ID corresponding to a given string.
 *
 */
	if(compare_string(str,"Number"))
		return NUMBER_ID;

	if(compare_string(str,"is5UTRFirst"))
		return UTR_ID;

	if(compare_string(str,"YES"))
		return YES_ID;
	else
		return NO_ID;

	if(compare_string(str,"Maximum_Iteration"))
		return MAX_ITER_ID;

	if(compare_string(str,"ID"))
		return ID_ID;

	if(compare_string(str,"Sequence"))
		return SEQUENCE_ID;

	if(compare_string(str,"Structure"))
		return STRUCTURE_ID;

	if(compare_string(str,"Structure_Target"))
		return TARGET_ID;

	if(compare_string(str,"Gate_Name")){
		cerr << "compared Gate_name true" << endl;
		return GATE_NAME_ID;
	}

	if(compare_string(str,"Concentration"))
		return CONCENTRATION_ID;

	if(compare_string(str,"Rows"))
		return ROW_ID;

	if(compare_string(str,"Columns"))
		return COL_ID;

	if(compare_string(str,"Begin"))
		return BEGIN_ID;

	if(compare_string(str,"End"))
		return END_ID;

	if(compare_string(str,"Name"))
		return NAME_ID;

	if(compare_string(str,"Truth_Table_Data"))
		return TRUTH_TABLE_DATA_HEADER_ID;

	if(compare_string(str,"RNA_Data"))
		return RNA_DATA_HEADER_ID;

	cout << "UNKNOWN_ID" << endl;
	return UNKNOWN_ID;
}


int read_params(string &str, string &str_value, int &value,string &strout){
	char	str_buffer[TMP_STRLEN];
	char	str_const_buffer[TMP_STRLEN];


	if(compare_string(str,"Number")){
		value = atoi(str_value.c_str());
		strout = "";
		return NUMBER_ID;
	}

	if(compare_string(str,"Is5UTRFirst")){
		value = compare_string(str_value,"YES") ? 1:0;
		strout = "";
		return UTR_ID;
	}

	if(compare_string(str,"Maximum_Iteration")){
		value = atoi(str_value.c_str());
		strout = "";
		return MAX_ITER_ID;
	}

	if(compare_string(str,"RNA_Data")){
		return RNA_DATA_HEADER_ID;
	}

	if(compare_string(str,"ID")){
		value = atoi(str_value.c_str());
		strout = "";
		return ID_ID;
	}

	if(compare_string(str,"Sequence")){
		value = -1;
		strout = str_value;
		return SEQUENCE_ID;
	}

	if(compare_string(str,"Structure")){
		value = -1;
		strout = str_value;
		return STRUCTURE_ID;
	}

	if(compare_string(str,"Structure_Target")){
		value = -1;
		strout = str_value;
		return TARGET_ID;
	}

	if(compare_string(str,"Gate_Name")){
		cerr << "compared Gate_name true" << endl;
		value = -1;
		strout = str_value;
		return GATE_NAME_ID;
	}

	if(compare_string(str,"Concentration")){
		value = (int) atof(str_value.c_str()) * 1000;
		strout = "";
		return CONCENTRATION_ID;
	}

	if(compare_string(str,"Name")){
		value = -1;
		strout = str_value;
		return NAME_ID;
	}

	if(compare_string(str,"Truth_Table_Data")){
		value = -1;
		strout = "";
		return TRUTH_TABLE_DATA_HEADER_ID;
	}

	if(compare_string(str,"Table_Size")){
		value = -1;
		strout = kill_all_spaces(str_value);
		return TABLE_SIZE_ID;
	}

	cout << "UNKNOWN_ID" << endl;
	return UNKNOWN_ID;
}



