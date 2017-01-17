# RNAGate
Non-coding small RNA reaction and structure calculation


============================================================================
 Name        : RNAgate.cpp
 Author      : Marzuk M Kamal
 Version     :
 Copyright   : (c) Marzuk M Kamal, 2010
 Description : RNAgate v1.1
============================================================================


  # RNAGate v 1.1
 
  Usage:
  rnagate inputfilename.txt outputfilename.txt -OPTIONS
 
  OPTIONS must be either or both of the following:
 
  f = write to file outputfilename.txt
  s = write to standard output, i.e. display the results on screen
 
  r = calculate the subsequence free probability....
 
  To write both to file and display results on screen
 
   rnagate inputfilename.txt outputfilename.txt -sf
  command line for subsequence free probability calculation
 
   rnagate inputfilename.txt outputfilename.txt -r
 
 
  Note that the input file format is different for -r option. See below.
 
 

 *********************************************************************************
 In order to compile the code, please use Eclipse C++ development tools

 You need to install gcc 4.5.0+ to compile and run the code
 
 I used gcc 4.6.0 release verison to compile the code. it's required for openmp library
 



 * input file format for OPTION = -s or -f or -sf
 %
 % comments
 % 5'UTR sequence must be the first one
 %
 % Please follow exact gaps between paragraphs
 % All inputs must start from the first column
 %
 
 N (number of RNAs involved including UTR)
 
 .....((((....))).... %(monomer structure constraint)
 AGCAGCGCAGGCTUAGAAUG %(initial sequence, if the initial sequnece is not present, put only a # in this line)
 NNNNNNNNNNGCTUNNNAUG %(monomer sequence constraint, the constraint bases must match with the initial sequence (if provided) in the above line)
 @@@@@@@@@@RRRRRRRRRR %(constraint when the RNAs make complexes, R means RBS or any associated base that affects translation) <br/>
 X			%(number of mutations sites for the UTR) <br/>
 D			%(maximum allowed distance of this monomer)<br/>
 Y			%(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken) <br/>
 
 .....((((....))).... (monomer structure constraint for sRNA1) <br/>
 \# <br/>
 NNNNNNNNNNNNNNNNNNNN %(monomer sequence constraint) <br/>
 @@@@@@@@@@@@@@@@@@@@ %(constraint when the RNAs make complexes) <br/>
 X			%(number of mutations sites for the RNA) <br/>
 D			%(maximum allowed distance of this monomer) <br/>
 Y			%(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken) <br/>
 .....((((....))).... %(monomer structure constraint sRNA2)<br/>
 \#<br/>
 NNNNNNNNNNNNNNNNNNNN (monomer sequence constraint)<br/>
 @@@@@@@@@@@@@@@@@@@@ (constraint when the RNAs make complexes)<br/>
 X			%(number of mutations sites for the sRNA2)<br/>
 D			%(maximum allowed distance of this monomer)<br/>
 Y			%(Y = 1 is the computation attempts to preserve, during mutations, the stems of the RNA structures, if Y = 0, no preservation steps taken)<br/>
 AND	(name of the predefined logic gate that we would like to achieve, the options are YES,NOT,AND, OR,XOR,NAND,XNOR,NOR)<br/>
 BASIC_REAPEAT_CHECK (or EXTENTED_SEQUENCE_CHECK)  %(Please see the documentation of check_valid_sequence)<br/>
 MONTE_CARLO_METHOD  (or SIMULATED_ANNEALING_METHOD)  %(choose the method of calculation)<br/>
  


 For calculating RBSfree probability:<br/>
 <br/>
 command line:<br/>
 
 rnagate free_in.txt rfree_out.txt -r<br/>
 <br/>
 the input file format.<br/>
 <br/>
 rfee.txt<br/>
 
2 % kcal/mol, sub-optimal energy gap, comments are allowed in the input file.<br/>
4 10 % start and end positions of the subsequence<br/>
% each line of the input file contains an RNA or RNA-complex<br/>
<br/>
% RNA sequence one<br/>
ACGUGUAGUGAUAGUGAGA<br/>
% another RNA-complex<br/>
ACGUGUAGUGAUAGUGAGA+GAUGAUAGUAGAUAGCAGC<br/>
ACGUGUAGUGAUAGUGAGA+GAUGAUAGUAGAUAGCAGC+CGCGGCGAGAGAG<br/>



 Before running RNAGate, the environmental variable RNAGATEHOME, must be set to the home folder of RNAGate which contains rnagate binary and the "parameters" folder<br/>

 Path tree:<br/>
 
 ~/
   RNAGate/
    		parameters/

command:

export RNAGATEHOME=~/RNAGate
export PATH=$PATH:$RNAGATEHOME

if not set, the default value of RNAGATEHOME will be

RNAGATEHOME=~/RNAGate


