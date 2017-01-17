
#ifndef __NUPACK_OBJECT_H
#define __NUPACK_OBJECT_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <iostream>

#include "pfuncUtilsHeader.h"
#include"DNAExternals.h"
//#include "hash.h"



#define		MAX_NUPACK_ARRAY_LENGTH		3000


//#define 	fold	fold_nupack
//#define  	DNA		DNA_1
//#define 	RNA		RNA_1
//#define		RNA37	RNA37_1
//#define		COUNT 	COUNT_1

//using namespace std;


void process_multi_sequence( char *inputSeq, int seqlength, int nStrands, char *seq, int *nicks);
int  compare_dna_structs(const void *p1, const void *p2);


class NUPACKObject{
private:

	int 	complexity_val;
	int 	naType_val;
	int 	dangle_val;
	DBL_TYPE temperature_val;
	int 	symmetry_val;
	DBL_TYPE fixedSubOptRange_val;
	DBL_TYPE *sizeTerm_ptr;



	short 	*possiblePairs;
	int 	**etaN;
	int 	arraySize;
	int 	nStrands;
	//	char 	*seq;
	int 	*seq;

	DBL_TYPE salt_correction;

	//	char 	*foldparens;
	int 	*foldparens;
	//	int		*nicks;

	char 	*sequence_nupack;
	char	*structure_nupack;
	double	 energy_nupack;
	dnaStructures 	mfeStructs_nupack;

	DBL_TYPE mfeEpsilon;
	DBL_TYPE *minILoopEnergyBySize;
	int 	*maxILoopSize;

	DBL_TYPE localEnergy;
	int 	oldp1;
	int 	symmetryOfStruct; // Must be initialized
	int 	*thepairs;
	int		seqlength;

	long int	maxIndex;

	int 	*nicks;//[ MAXSTRANDS];  //the entries must be strictly increasing
	//nicks[i] = N means a strand ends with base N, and a new one starts at N+1

	int 	calcPairs;
	int 	permSymmetry;

	DBL_TYPE *F;
	DBL_TYPE *Fb;
	DBL_TYPE *Fm;

	//N^3
	DBL_TYPE *Fx;
	DBL_TYPE *Fx_1;
	DBL_TYPE *Fx_2;
	DBL_TYPE *Fs;
	DBL_TYPE *Fms;

	//PKNOTS
	DBL_TYPE *Fp;
	DBL_TYPE *Fz;  //O(N^2)
	DBL_TYPE *Fg; //O(N^4)

	//N^5
	DBL_TYPE *FgIx;
	DBL_TYPE *FgIx_1;
	DBL_TYPE *FgIx_2;
	DBL_TYPE *Fgls;
	DBL_TYPE *Fgrs;

	DBL_TYPE *Fgl;
	DBL_TYPE *Fgr; //O(N^4) space



	DBL_TYPE *Q;
	DBL_TYPE *Qb;
	DBL_TYPE *Qm; //O(N^2)

	//Multiple strand arrays


	//N^3 arrays
	DBL_TYPE *Qx;
	DBL_TYPE *Qx_1;
	DBL_TYPE *Qx_2;
	DBL_TYPE *Qs;
	DBL_TYPE *Qms;

	//Pseudoknot arrays
	DBL_TYPE *Qp;
	DBL_TYPE *Qz; //O(N^2)
	DBL_TYPE *Qg; //O(N^4) space

	DBL_TYPE *QgIx;
	DBL_TYPE *QgIx_1;
	DBL_TYPE *QgIx_2;
	DBL_TYPE *Qgls;
	DBL_TYPE *Qgrs; //O(N^4)
	DBL_TYPE *Qgl;
	DBL_TYPE *Qgr; //O(N^4) space

	//extern DBL_TYPE *sizeTerm;

	//Pair probabilities
	DBL_TYPE *Pb;
	DBL_TYPE *P;
	DBL_TYPE *Pm;
	DBL_TYPE *Pms;
	DBL_TYPE *Ps;

	//pseudoknots
	DBL_TYPE *Pz;
	DBL_TYPE *Pp;
	DBL_TYPE *Pg;
	DBL_TYPE *Pbg;

	//N^5
	DBL_TYPE *Pgl;
	DBL_TYPE *Pgr;
	DBL_TYPE *Pgls;
	DBL_TYPE *Pgrs;


	int 	isPairPrExtern;

	int		max_gap_index_value;
public:

	NUPACKObject();
	NUPACKObject(int *inputSeq, int seqLen,
			int complexity, int naType,
			int dangles, DBL_TYPE temperature, int symmetry,DBL_TYPE fixedSubOptRange,int calc_Pairs, int perm_Symmetry,
			DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
			int uselongsalt);
	~NUPACKObject();
	//	bool init_nupack_data();
	bool init_nupack_data(int *inputSeq, int seqLen,
			int complexity, int naType,
			int dangles, DBL_TYPE temperature, int symmetry,DBL_TYPE fixedSubOptRange,int calc_Pairs, int perm_Symmetry,
			DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
			int uselongsalt);

	DBL_TYPE mfe_sym_subopt( int *inputSeq, dnaStructures *mfeStructures,DBL_TYPE fixedSubOptRange,int onlyOne);
	DBL_TYPE mfe_sym_subopt_with_pfunc( int *inputSeq, dnaStructures *mfeStructures,DBL_TYPE fixedSubOptRange = 1,
									    int onlyOne = 1, DBL_TYPE *pfunc_val = NULL);

	DBL_TYPE mfe_sym_subopt_with_pfunc2( int *inputSeq, dnaStructures *mfeStructures,DBL_TYPE fixedSubOptRange = 1,
									    int onlyOne = 1, DBL_TYPE *pfunc_val = NULL);


	DBL_TYPE mfe_sym_suboptmfe_sym_subopt(int *inputSeq, dnaStructures *mfeStructures,DBL_TYPE fixedSubOptRange,int onlyOne);

	DBL_TYPE pfunc(int inputSeq[]);

	DBL_TYPE probability(int *inputSeq);

	//	DBL_TYPE mfeFullWithSym_SubOpt_Progress( char inputSeq[], dnaStructures *mfeStructures, int complexity, int naType,
	//			int dangles, DBL_TYPE temperature, int symmetry, DBL_TYPE fixedSubOptRange,
	//			double prog, char *progName, int onlyOne);
	//
	//	DBL_TYPE mfe( char *seq, int *thepairs);
	//	DBL_TYPE mfeFull( char inputSeq[], int *thepairs, int complexity, int naType,
	//			int dangles, DBL_TYPE temperature, int onlyOne);
	//	DBL_TYPE mfeFullWithSym( char inputSeq[], dnaStructures *mfeStructures, int complexity, int naType,
	//			int dangles, DBL_TYPE temperature, int symmetry, int onlyOne);
	//	DBL_TYPE mfeFullWithSym_Progress( char inputSeq[], dnaStructures *mfeStructures, int complexity, int naType,
	//			int dangles, DBL_TYPE temperature, int symmetry, double prog, char  *progName,
	//			int onlyOne);
	void initMfeStructures( dnaStructures *mfeStructures, int seqlength);
	void clearDnaStructures( dnaStructures *ds);

	void InitLDoublesMatrix( DBL_TYPE **Q, int size, char name[]);
	void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]);
	bool delete_vector(DBL_TYPE **data_ptr);

	void nonZeroInit( DBL_TYPE Q[], int seq[], int seqlength);
	void manageQgIx( DBL_TYPE **QgIx, DBL_TYPE **QgIx_1, DBL_TYPE **QgIx_2, int d, int seqlength);
	void manageFgIx( DBL_TYPE **FgIx, DBL_TYPE **FgIx_1, DBL_TYPE **FgIx_2, int d, int seqlength);

	int compareDnaStructs( const void *p1, const void *p2);
	void make_rna_structure(const dnaStructures &ds,const char *seq,char *str, int index = 0);
	char **make_multiple_rna_structures(const dnaStructures &ds,const char *seq,int &number_of_structures);
	void process_multi_sequence( char *inputSeq, int seqlength, int nStrands, char *seq, int *nicks);
	int  convertSeq(const char *seqchar, int *seqnum, int seqlength);
	void initPF( int seqlength);
	void initMfe( int seqlength);

	void MakeFs_Fms_and_Qs_Qms( int i, int j);
	void MakeF_Fm_N3_and_Q_Qm_N3( int i, int j);
	void NickedEmptyF_n_Q( int i, int j,DBL_TYPE &F_pf, DBL_TYPE &Q_pf);



//	extern long int maxGapIndex;

	/* ***************************** */
//	void PrecomputeValuesN5( int seqlength);
//	void PrecomputeValuesN5f( int seqlength);
//	void manageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1, DBL_TYPE **Qx_2, int len, int seqlength);
//	void manageFx( DBL_TYPE **Fx, DBL_TYPE **Fx_1, DBL_TYPE **Fx_2, int len, int seqlength);

};


#endif
