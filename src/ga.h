/*===============================================================================
ga : Genetic Algorithm
(C) 2003-2006 Jens Kleinjung
Read the COPYING file for license information.
*==============================================================================*/

#if !defined GA_H
#define GA_H

/*____________________________________________________________________________*/
/* includes */
/*
#define MPI "/usr/local/apps/mpich/mpich-1.2.2.3/include/mpi.h"
#include MPI
*/

#include <alloca.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
/*#include "aliss.h"*/

/*____________________________________________________________________________*/
/* define debug mode*/
/* #define DEBUG */

/* define alignment mode*/
//#define ALIGNMODE

/* define compress score */
/* #define COMPRESS_SCORE */

/* define bitfield data structure for genome */
/*#define BIT*/

/* define stuctural code */
#define SA_ALPHABET

/* define generation of background distribution for most fitted */
#define GENERATE_BACKGROUND

/* define benchmark database */
/*#define HOMSTRAD*/
/*#define CATH*/

/* parameters */
#define POPSIZE 100 /* population size */
#define FITMATE  5 /* number of fittest individuals to mate */
#define PARNUM 7/* number of parameters to optimize */
#define LOOP 2500 /* number of optimisation cycles */
#define UPLIM 1 /* upper limit of parameters */
#define LOWLIM 0 /* lower limit of parameters */

#define BREEDMODE 0 /* breeding mode: crossover(0) or equilibrium (1) */
#define JACKKNIFE 1 /* splitting of database in 'JACKKNIFE' parts (1 = no jackknife) */
#define REPGA 1 /* repeat entire GA 'REPGA' times (1 = no repeat) */

#define NUMRAN_POOL 2500 /* number of randomly generated pools for background distribution */
/* max macro */
#define max(a,b)  (((a) > (b)) ? (a) : (b))

/*____________________________________________________________________________*/
typedef struct
{
    char *protein_name; /* protein filename */
	char *aa_description; /* description in header of fastafile */
    char *aaseq; /* aaseq from fastafile */
	float aaseq_entropy; /* aa entropy */
	float aaseq_score; /* aa score */
	char *SA_description; /* description in header of fastafile */
    char *SAseq; /* SAseq from fastafile */
	float SAseq_entropy; /* sa entropy */
	float SAseq_score; /* sa score */
} ProteinEntry;

typedef struct
{
	ProteinEntry *protein_entries; /* protein data */
	int protnum; /*number of proteins */
	float aa_entropy_sum; /* sum of single protein entropy in aa code */
	float SA_entropy_sum; /* sum of single protein entropy in SA code */
} Prots;

/*____________________________________________________________________________*/
/* JK: replaced by 'freq' routines
typedef struct
{
	char aa;
	float freq;
} AminoacidFreq; */

/*____________________________________________________________________________*/
/* bit field for binary gene */
#ifdef BIT
typedef struct
{
	unsigned int state: 1; /* 1 bit integer */
} Bitgene;
#endif

typedef struct
{
#ifdef BIT
	Bitgene *bitgenome; /* bitgenome is an array of binary parameters=bitgenes */
#else
    int *genome; /* genome is an array of parameters=genes */
#endif
    float fitness; /* fitness of genome */
} Pool;

typedef struct
{
    float value;
    char name[64];
} Parameter;

/*____________________________________________________________________________*/
typedef struct
{
    char *family[2000]; /* max. 2000 family names */
    int nfam; /* number of families */
} Fam;

/*____________________________________________________________________________*/
FILE *safe_open(const char *name, const char *mode);
extern void *safe_malloc(size_t), *safe_realloc(void *, size_t);
float aliss(char *pdbfilename0, char *pdbfilename1, Parameter *parameter);

#endif
