/*
 * nupack_mfe_pfunc.h
 *
 *  Created on: Feb 28, 2011
 *      Author: kamal
 */

#ifndef NUPACK_MFE_PFUNC_H_
#define NUPACK_MFE_PFUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "pfuncUtilsHeader.h"

/* Declare Global Variables. */
#include"DNAExternals.h"

#include "global_variables.h"

//int containsPk;
/* End Global Variables */


extern DBL_TYPE mzk_mfeFullWithSym_SubOpt( int inputSeq[], int seqLen,
		dnaStructures *mfeStructures, int complexity, int naType,
		int dangles, DBL_TYPE temperature, int symmetry, DBL_TYPE fixedSubOptRange,
		int onlyOne, DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
		int uselongsalt);
extern DBL_TYPE mzk_mfe( int seq[], int seqLen, int *thepairs);
extern DBL_TYPE mzk_mfeFullWithSym( int inputSeq[], int seqLen,
		dnaStructures *mfeStructures, int complexity, int naType,
		int dangles, DBL_TYPE temperature, int symmetry, int onlyOne,
		DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);
extern void mzk_initMfeStructures( dnaStructures *mfeStructures, int seqlength);
extern DBL_TYPE mzk_mfeFull( int inputSeq[], int seqLen, int *thepairs, int complexity,
		int naType, int dangles, DBL_TYPE temperature,
		DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

extern int mzk_compareDnaStructs( const void *p1, const void *p2);
extern int mzk_compareDnaStructsOutput(const void *p1, const void * p2);



extern DBL_TYPE mzk_pfuncFullWithSymHelper( int inputSeq[], int seqlength, int nStrands,
                                int complexity, int naType, int dangles, DBL_TYPE temperature,
                                int calcPairs, int permSymmetry, DBL_TYPE sodiumconc,
                                DBL_TYPE magnesiumconc, int uselongsalt,DBL_TYPE **pair_probability = NULL);

double calc_pair_probability(int row,int col, DBL_TYPE *pair_pr,int seq_len);

extern DBL_TYPE mzk_pfuncFullWithSym( int inputSeq[], int complexity, int naType,
			   int dangles, DBL_TYPE temperature, int calcPairs, int permSymmetry,
			   DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt,DBL_TYPE **pair_probability = NULL);

extern DBL_TYPE mzk_pfunc( int seq[]);
extern DBL_TYPE mzk_pfuncFull( int inputSeq[], int complexity, int naType, int dangles,
                    DBL_TYPE temperature, int calcPairs,
		    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt,DBL_TYPE **pair_probability = NULL);

extern DBL_TYPE mzk_ExplDangle( int i, int j, int seq[], int seqlength, unsigned int seqHash);
void mzk_InitEtaN( int **etaN, const int *nicks, int seqlength);


#endif /* NUPACK_MFE_PFUNC_H_ */
