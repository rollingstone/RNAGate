/*
  complexesHeader.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 8/2006 and Justin Bois 1/2007

  Header file for use with complexes.c.
 */

#include "pfuncUtilsHeader.h"
#include "complexesStructs.h"
#include "physical_constants.h"
#include "runtime_constants.h"
#include "DNAExternals.h"



#define 	fold	fold_nupack
#define  	DNA		DNA_1
#define 	RNA		RNA_1
#define		RNA37	RNA37_1
#define		COUNT 	COUNT_1
#ifndef FUNCTIONSCX_H
#define FUNCTIONSCX_H


// Functions from complexesUtils.c
int getPermutation( int, int, int*);
void resetNicks( int, int*);
void nextMultiset( int, int*, int*, int*, int); //generate a multiset
int isCyclicP(int, int*, int*); //check if a cyclic permutation
void symmetryCheck( multiset*, int, permutation*); //check if symmetry
void printPerms( FILE*, int, int, multiset*); 
void printMfesToFile( const dnaStructures *ds, FILE *fp, 
		      const int *nicks);
int compareMultisets( const void*, const void*);


//Functions for permBG.c
void neg( int t, int n, int k);
void gen( int t, int n, int k);
void PrintIt( void);
void BigGen( int);
void swap( int, int, int);
void initializeMP( int, int*);
void freeMP( void);
int nextPerm( void);
void setPerm( int*);

// Functions from ReadCommandLine.c
int ReadCommandLine( int, char**); 
int ReadInputFileComplexes(  char *filePrefix, int *nStrands, 
			     char ***seqs, int **seqlength,
			     int *maxLength, int *maxComplexSize);
void printHeader( int nStrands, char **seqs, int maxComplexSize, 
		  int totalOrders, int nSets, int nNewComplexes, 
		  FILE *F_cx, int nargs, char **argv, int isPairs);

// Functions from utils.c in the shared directory
int gcd(int a, int b);
double factorial(int n);
long double factorial_long(int n);
long GetRandSeed(long s);

#endif
