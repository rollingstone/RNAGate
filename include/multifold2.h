#ifndef MULTIFOLD2_H
#define MULTIFOLD2_H

#define MAXNUMSEQ      50
#define MAXSLEN         5000 
#define MAXSUBSTR       100


double multifold_ordered (int num_sequences, char sequences[MAXNUMSEQ][MAXSLEN], char structure[MAXSUBSTR*MAXSLEN]);


#endif


