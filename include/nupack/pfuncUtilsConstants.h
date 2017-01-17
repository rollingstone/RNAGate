/* 
  pfuncUtilsConstants.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2001, Justin Bois 1/2007

  Useful constants for running partition function applications.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "physical_constants.h"
#include "runtime_constants.h"

//sets the type of floating point variables
#define DBL_TYPE long double

//Minimum difference between energies before two are considered identical
#define ENERGY_TOLERANCE 0.0001

//max error in the bits of precision.  Used during pair probability
//calculations (where subtraction occurs) Setting this to zero can
//significantly slow down pair probability calculations.
#define MAXPRECERR 24 //max error in bits of precision

//Maximum seqeuence length
#define MAXSEQLENGTH 10000

//maximum # of strands in a complex
#define MAXSTRANDS 2000
 
//MATCH_PF will make the energy model used in energy calculations
//match the one used in mfe and partition function calculations.
//Otherwise, the energy of multiloops scales with the log of the size,
//rather than linearly.
//Other refinements, such as coaxial stacking, could also be included.
#define MATCH_PF

//Including NOGU will disallow all wobble pairs
//#define NOGU

//Including STRUCTURE_WARNINGS will produce error messages whenever
//a disconnected or illegal structure is evaluated, rather than just
//returning an "infinite" energy.
//#define STRUCTURE_WARNINGS

//Including DEBUG will cause various intermediate values to be printed during backtracking (backtrack.c)
//#define DEBUG

//Including FILEOUTPUT will cause pair probabilities to be spit out to a file named Pb_N#.txt, where
//the # is replaced with the time-complexity of the algorithm.

//#define FILEOUTPUT

//Including PRINTRESULTSONLY will limit output to the screen.  This needs to be updated.
//#define PRINTRESULTSONLY

//Including NODANGLES will set all dangle energies to zero.  Good for debugging purposes.
//#define NODANGLES

/* No changes below this line! */

//some physical constants

#define BASE_A 1  //do not change!
#define BASE_C 2  // Watson crick bp add to 5
#define BASE_G 3  // Wobble bp add to 7
#define BASE_T 4

//Multi base restrictions
#define BASE_AG 5 //R
#define BASE_CT 10 //Y
#define BASE_AC 6 //M
#define BASE_GT 9 //K
#define BASE_CG 7 //S
#define BASE_AT 8 //W

#define BASE_ACG 11 //V
#define BASE_ACT 12 //H
#define BASE_AGT 13 //D
#define BASE_CGT 14 //B

#endif
