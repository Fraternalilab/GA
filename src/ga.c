/*===============================================================================
ga : Genetic Algorithm
(C) 2003-2006 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <assert.h>
#include <math.h>

#include "freqs.h"
#include "ga.h"
#include "getfreqs.h"
#include "getseqs.h"
#include "huffman.h"
#include "lz.h"
#include "suffix_tree.h"
#include "safe.h"
#include "time.h"

const int popsize = POPSIZE;
const int fitmate = FITMATE;
const int parnum = PARNUM;
const int loop = LOOP;

const int uplim = UPLIM;
const int lowlim = LOWLIM;

const int breedmode = BREEDMODE;
const int jackknife = JACKKNIFE;
const int repga = REPGA;

const int numran_pool = NUMRAN_POOL;

int nfam = 0;
int genenum; /* number of gene per chromosome */
int k_word_len;
int freqset = 0;

/*-----------------------------------------------------------------------------*/
/* generate random integer number within range  0 to 'uplim' */
int get_rand(int limit)
{
    double dranum; /* random number in double format */
    int iranum; /* random number in integer format */
    dranum = rand()/(double)INT_MAX; /* generate rand and scale to range 0-1 */
    iranum = (int)((limit+1) * dranum); /* scale to integer of range 0-limit */
    if (iranum == (limit+1)) -- iranum; /* exclude upper limit */
    return iranum;
}

/*-----------------------------------------------------------------------------*/
/* fitness comparison routine used by 'qsort'  */
int cmp_fitness(const void *pp1, const void *pp2)
{
    Pool *p1 = (Pool*)pp1;
    Pool *p2 = (Pool*)pp2;
#ifdef ALIGNMODE
    return (p1->fitness < p2->fitness ? -1 : 1);  /* MINIMIZE */
#endif
#ifndef ALIGNMODE 
    return (p1->fitness > p2->fitness ? -1 : 1); /* MAXIMIZE */
#endif
}

/*-----------------------------------------------------------------------------*/
/* sort pool by fitness */
void sort_fitness(Pool *pool)
{
    qsort(&pool[0], popsize, sizeof(Pool), cmp_fitness);
}

/*-----------------------------------------------------------------------------*/
/* bitfield genome: set gene state to either 0 or 1, depending on passed integer 'i' */
#ifdef BIT
void set_bitgene(Bitgene *bitgen, int i) 
{
	bitgen->state = i > 0 ? 1 : 0;
}
#endif

/*-----------------------------------------------------------------------------*/
/* bitfield genome: return integer 0 or 1, depending on gene bit 'state' */
#ifdef BIT
int get_bitgene(Bitgene *bitgen) 
{
	return (int)bitgen->state > 0 ? 1 : 0;
}
#endif

/*-----------------------------------------------------------------------------*/
/* vary value */
int vary_value(int val, int maxvar)
{
    int newval; /* new value */
    int mp; /* minus/plus */

    mp = get_rand(2);

    if (mp == 1)
        while ((newval = val - get_rand(maxvar)) < lowlim);
    else
        while ((newval = val + get_rand(maxvar)) > uplim);
    return newval;
}

/*-----------------------------------------------------------------------------*/
/* initialise random pool and fitness*/
static void init_random_pool(Pool *pool)
{
    unsigned int i, j;

    for (i = 0; i < popsize; ++i) /* genomes (individuals) */
    {
        pool[i].fitness = 0; /* initialise fitness */
        for (j = 0; j < genenum; ++j) /* genes */
		{
#ifdef BIT
            set_bitgene(&pool[i].bitgenome[j], get_rand(uplim)); /* random value */
#else
            pool[i].genome[j] = get_rand(uplim); /* random value */
#endif
		}
    }
}

/*-----------------------------------------------------------------------------*/
/* initialise seeded pool and fitness */
void init_seeded_pool(Pool *pool)
{
    unsigned int i, j;
    int maxvar; /* maximal variation around equilibrium parameter */

    for (i = 0; i < popsize; ++i) /* genomes (individuals) */
    {
		pool[i].fitness = 0; /* initialise fitness */
		/* assign defined value to each gene */
		/*
		pool[i].genome[4] = 4; //gapo
		pool[i].genome[5] = 4; //gape
		*/
    }
    
    maxvar = 4; /* maxmum variation around seed value */

    /* vary values, not first genome */
    for (i = 1; i < popsize; ++i)
    {
		for (j = 0; j < genenum; ++j)
		{
#ifdef BIT
			set_bitgene(&pool[i].bitgenome[j], vary_value(get_bitgene(&pool[i].bitgenome[j]), maxvar));
#else
			pool[i].genome[j] = vary_value(pool[i].genome[j], maxvar);
#endif
		}
    }
}

/*-----------------------------------------------------------------------------*/
/* equilibrium: calc average parameter values in fittest genomes */
void equilibrium(Pool *pool, int *average)
{
    unsigned int i, j;

    for (j = 0; j < genenum; ++j)
    {
		average[j] = 0;
		for (i = 0; i < fitmate; ++i) /* sum up over the fittest genomes */
		{
#ifdef BIT
			average[j] += get_bitgene(&pool[i].bitgenome[j]);
#else
			average[j] += pool[i].genome[j];
#endif
		}
		average[j]  = (int)(average[j]/fitmate); /* normalize */
    }
}

/*-----------------------------------------------------------------------------*/
/* breed_equilibrium: generate equilibrium genomes */ 
void breed_equilibrium(Pool *pool, int *average, int ix)
{
    unsigned int j;
    int maxvar; /* maximal variation around equilibrium parameter */

    maxvar = uplim/loop; /* decay to zero variation */

    for (j = 0; j < genenum; ++j)
	{
#ifdef BIT
		set_bitgene(&pool[ix].bitgenome[j], vary_value(get_bitgene(&pool[ix].bitgenome[j]), maxvar));
#else
		pool[ix].genome[j] = vary_value(pool[ix].genome[j], maxvar);
#endif
	}
    pool[ix].fitness = 0; /* initialise fitness */
}

/*-----------------------------------------------------------------------------*/
/* calculate sum of genome: n of selected items */
int gene_sum(Pool *pool, int ix)
{
        int gsum = 0;
        int i;

        for(i = 0; i < genenum; i++)
		{
#ifdef BIT
			gsum += get_bitgene(&pool[ix].bitgenome[i]);
#else
            gsum += pool[ix].genome[i];
#endif
		}

        return gsum;
}

/*-----------------------------------------------------------------------------*/
/* constrain individual genome to minset_perc */
void constrain_genome(Pool *pool, int ix, int nselected)
{
    int j;

    while(gene_sum(pool, ix) > nselected) /* check that number of selected does not exceed requested */
    {
        /* deselect a random protein*/
#ifdef BIT
        do {
			j = (int) (rand()/(double)INT_MAX * genenum);
        } while(get_bitgene(&pool[ix].bitgenome[j]) == 0);
        set_bitgene(&pool[ix].bitgenome[j], 0); 
#else
        do {
			j = (int) (rand()/(double)INT_MAX * genenum);
        } while(pool[ix].genome[j] == 0);
        pool[ix].genome[j] = 0; 
#endif
    }
}

/*-----------------------------------------------------------------------------*/
/* breed_crossover: crossover between fitmate pairs and create genome 'ix' */ 
void breed_crossover(Pool *pool, int iy)
{
    int ia, ib, j, yn;

    ia = -1;
    ib = -1;
    for (j = 0; j < genenum; ++j)
    {
		ia = get_rand(fitmate); /* principal parent (1)*/

		do
			ib = get_rand(fitmate); /* crossover parent (2) */
		while (ib == ia);

		yn = get_rand(2); /* choose crossover yes/no */
#ifdef BIT
		if (yn == 1)
			set_bitgene(&pool[iy].bitgenome[j], get_bitgene(&pool[ib].bitgenome[j]));
		else
			set_bitgene(&pool[iy].bitgenome[j], get_bitgene(&pool[ia].bitgenome[j]));
#else
		if (yn == 1)
			pool[iy].genome[j] = pool[ib].genome[j];
		else
			pool[iy].genome[j] = pool[ia].genome[j];
#endif
    }

    pool[iy].fitness = 0; /* initialise fitness */
}

/*-----------------------------------------------------------------------------*/
/* transform gene values to run parameters */
void set_parameters(Pool *pool, int ix, Parameter *parameter, FILE *gaoutfile)
{
    unsigned int i = 0;

    /* parameter[0] -> gap open */
#ifdef BIT
    parameter[i].value = (float)((0.5 * get_bitgene(&pool[ix].bitgenome[i])) + 8.);
#else
    parameter[i].value = (float)((0.5 * pool[ix].genome[i]) + 8.);
#endif
    strcpy(parameter[i].name, "go");

    fprintf(stdout, "\t%d ", ix);
    fprintf(gaoutfile, "\t%d ", ix);
    for (i = 0; i < parnum; ++ i)
    {
#ifdef BIT
		fprintf(stdout, "%s:%d->%5.3f ", parameter[i].name, get_bitgene(&pool[ix].bitgenome[i]), parameter[i].value);
		fprintf(gaoutfile, "%s:%d->%5.3f ", parameter[i].name, get_bitgene(&pool[ix].bitgenome[i]), parameter[i].value);
#else
		fprintf(stdout, "%s:%d->%5.3f ", parameter[i].name, pool[ix].genome[i], parameter[i].value);
		fprintf(gaoutfile, "%s:%d->%5.3f ", parameter[i].name, pool[ix].genome[i], parameter[i].value);
#endif
    }
    /*
    fprintf(stdout, "\n");
    fprintf(gaoutfile, "\n");
    */
    fflush(stdout);
    fflush(gaoutfile);
}

/*-----------------------------------------------------------------------------*/
/* mem_fitness: get fitness from identical genome that has already been calculated*/
int mem_fitness(Pool *pool, int ix)
{
    unsigned int i, j;
    unsigned int id;

    for (i = 0; i < ix; ++ i)
    {
		id = 1;

		/* test whether any gene of 'i' and 'ix' differs */
		for (j = 0; j < genenum; ++ j)
		{
#ifdef BIT
			if (get_bitgene(&pool[i].bitgenome[j]) != get_bitgene(&pool[ix].bitgenome[j]))
				id = 0;
#else
			if (pool[i].genome[j] != pool[ix].genome[j])
				id = 0;
#endif
		}
		/* if identical and already fitnessed: copy fitness */
		if ((id == 1) && (pool[i].fitness > 0))
		{
			pool[ix].fitness = pool[i].fitness;
			return 1;
		} 
    }

    return 0;
}
/*-----------------------------------------------------------------------------*/
/* print pool bin */
/*JK void print_pool_bin(Pool *pool, FILE *pooloutfile)
{
    unsigned int i, j;
    char ci;

    fprintf(stdout, "\n");

    for (i = 0; i < popsize; ++i)
    {
        fprintf(stdout, "%3d: ", i);
        fwrite(&i, sizeof(int), 1, pooloutfile);

        for (j = 0; j < genenum; ++j)
        {
#ifdef BIT
            ci = get_bitgene(&pool[i].bitgenome[j]) ? '1':'0';
            fprintf(stdout, "%1d ", get_bitgene(&pool[i].bitgenome[j]));
#else
            ci = pool[i].genome[j] ? '1':'0';
            fprintf(stdout, "%1d ", pool[i].genome[j]);
#endif
            fwrite(&ci, sizeof(char), 1, pooloutfile);
        }

        fprintf(stdout, "%6.4f\n", pool[i].fitness);
        fwrite(&pool[i].fitness, sizeof(float), 1, pooloutfile);

    }
    fprintf(stdout, "\n");
    fflush(stdout);
    fflush(pooloutfile);
}
*/

/*-----------------------------------------------------------------------------*/
/* print pool ascii */
/* JK
void print_pool_ascii(Pool *pool, FILE *pooloutfile)
{
    unsigned int i, j;

    fprintf(stdout, "\n");
    fprintf(pooloutfile, "\n");

    for (i = 0; i < popsize; ++i)
    {
		fprintf(stdout, "%3d: ", i);
		fprintf(pooloutfile, "%3d: ", i);

		for (j = 0; j < genenum; ++j)
		{
#ifdef BIT
			fprintf(stdout, "%1d ", get_bitgene(&pool[i].bitgenome[j]));
			fprintf(pooloutfile, "%1d ", get_bitgene(&pool[i].bitgenome[j]));
#else
			fprintf(stdout, "%1d ", pool[i].genome[j]);
			fprintf(pooloutfile, "%1d ", pool[i].genome[j]);
#endif
		}

		fprintf(stdout, "%6.4f\n", pool[i].fitness);
		fprintf(pooloutfile, "%6.4f\n", pool[i].fitness);

    }
    fprintf(stdout, "\n");
    fflush(stdout);
    fprintf(pooloutfile, "\n");
    fflush(pooloutfile);
}
*/

/*-----------------------------------------------------------------------------*/  
/* get code counts*/
void get_counts(char *seq, int *count, int alphabet_array_len)
{
    int i;

    for(i=0; i< alphabet_array_len; i++)
        count[i] = 0;

    /*while((*seq) != 0)*/
    while((*seq) != '\0') /* JK: string end? */
    {
		/*assert(((*seq)-'A') > 0 && ((*seq)-'A') < alphabet_array_len);*/
        ++ count[(*seq)-'A'];
        ++ seq;
    }
}

/*-----------------------------------------------------------------------------*/  
/* Shannon formula */
float shannon(float p)
{
    return (-1 * (p * log(p) / log(2)));
}

/*-----------------------------------------------------------------------------*/  
/* Shannon entropy of code */
float shannon_entropy(float *p, int length)
{
    unsigned int i;
    float H;

    for (i = 0, H = 0; i < length; ++ i)
    {
        if (p[i] == 0)
            continue;
        H += shannon(p[i]);
    }

    return H;
}

/*-----------------------------------------------------------------------------*/  
/* relative entropy of code (Kullback-Leibler distance) */
float relative_entropy(float *p, float *q, int length)
{
    unsigned int i;
    float D;

    for (i = 0, D = 0; i < length; ++ i)
    {
        if (p[i] == 0 || q[i] == 0)
            continue;
        D += p[i] * log(p[i] / q[i]) / log(2);
    }

    return D;
}

/*-----------------------------------------------------------------------------*/  
/* word entropy from suffix tree: */
/* the entropy of the subset based on a constant-length word alphabet */

float word_entropy(char *ref_string, char *search_string, int *ptr_nsymbols)
{
    SUFFIX_TREE* tree; /* the suffix tree */
    DBL_WORD ref_len = 0; /* length of reference string */
    DBL_WORD search_len = 0; /* length of search (sub)string */
	extern int k_word_len;
    DBL_WORD sub_len = (DBL_WORD)k_word_len; /* alphabet word length */
	/*DBL_WORD sub_len = 5;*/
	unsigned int i = 0;
    char *sub_pc; /* pointer to char for search substring */
    unsigned int n_all = 0; /* total number of substrings */
    DBL_WORD position; /* position of search (sub)string in reference string */
    float entropy;

    ref_len = (DBL_WORD)strlen(ref_string);
    search_len = (DBL_WORD)strlen(search_string);

    tree = ST_CreateTree(ref_string, ref_len); /* generate the suffix tree */
    ST_InitTreeHits(tree);

    /* for the length of the search string */
    while (i < search_len - (sub_len - 1))
    {
        /* move pointer along search string to create substrings */
        sub_pc = search_string + i;

        /* skip all substrings containing the inter-sequence delimiter '-' */
        if (sub_pc[sub_len - 1] == '-')
        {
            i += sub_len; /* increment pointer position by substring length */
            continue;
        }

        ++ n_all; /* count all subsequences */

        /* search substring in suffix tree */
        position = ST_FindSubstring(tree, sub_pc, sub_len);

        ++ i; /* increment pointer position */
    }

    assert(n_all == (tree->allhit + tree->allmiss));

    entropy = ST_TreeEntropy(tree);
    *ptr_nsymbols = tree->allsymbol; 

    ST_DeleteTree(tree);

    return entropy;
}

/*-----------------------------------------------------------------------------*/  
/* score the seq */
float score_seq(char *subseq, char* search_polyfasta, char *superseq, FreqSet *bg_freq, int alphabet_array_len, float *ptr_KL_distance)
{
    int i;
    int strlength, strlength1;
    float H = 0.; /* entropy */
/*  float E = 0.; *//* expected entropy */
    float D = 0.; /* relative entropy */
    float score = 0.; /* fitness score */
    int nsymbols;
    int *count;
    float *p_count;
    float *p_bg;

    count = safe_malloc(sizeof(int) * alphabet_array_len);
    p_count = safe_malloc(sizeof(float) * alphabet_array_len);
    p_bg = safe_malloc(sizeof(float) * alphabet_array_len);

    strlength = strlen(subseq); /* count before moving the seq pointer on*/
    strlength1 = strlen(superseq);  /* count before moving the seq pointer on*/

    get_counts(subseq, count, alphabet_array_len);
    
    for(i =0; i < alphabet_array_len; ++i)
        p_count[i] = count[i] / (float) strlength;

    for(i =0; i < alphabet_array_len; ++i) {
        p_bg[i] = 0; /* init before assigning, to set unused characters */
        p_bg[i] = bg_freq->freq[i];
	}

    /*H = shannon_entropy(p_count, alphabet_array_len);*/ /* single-character alphabet */
    
    H = word_entropy(subseq, search_polyfasta, &nsymbols); /* contant-length word alphabet */
    D = relative_entropy(p_count, p_bg, alphabet_array_len);
/*  E = shannon_entropy(p_bg, alphabet_array_len); 
    score = H/E * (1 - D); */  /* first version of scoring: scaled by expected code entropy */

    score = H * (1 - D);

    *ptr_KL_distance = D;

#ifdef DEBUG
    fprintf(stderr, "l: %d, H: %f, E: %f, D: %f, score: %f\n", 
            strlength, H, E, D, score);
#endif

    free(count);
    free(p_count);
    free(p_bg);

    return score;
}

/*-----------------------------------------------------------------------------*/  
/* score the seq by compression ratio */
float score_compress(char *instring, int insize, int total_len)
{
    unsigned char *out;
    unsigned int *work;
    unsigned int outsize; 
    int maxoutsize;
    float compression_ratio;
    float zipscore;
    
    /* maximum output size w/ lz : (257/256)*insize + 1. */
    maxoutsize = ceil(257.0/256.0*insize) + 1;

/*  work - Pointer to internal working buffer, which must be able to hold (insize+65536) unsigned integers.*/
    work = safe_malloc(sizeof(unsigned int) * (insize + 65536));
    out = safe_malloc(sizeof(unsigned char) * maxoutsize);

/*  intersize = Huffman_Compress(instring, out, insize); */
    outsize   = LZ_CompressFast((unsigned char *) instring, out, insize, work); 

    if (insize == 0)
    {
        free(work);
        free(out);
        return 0.0; /* exclude null chromosome cases */
    }

    compression_ratio = (float) (insize - outsize) / insize; 
    compression_ratio = max(compression_ratio, 0.1);

    zipscore = (log(1 + insize / total_len) / log(2)) / compression_ratio; 

#ifdef DEBUG
        fprintf(stderr, "--> %d\t%d\t %d\t%5.1f\t%5.1f\t", insize, outsize, total_len, compression_ratio, zipscore*100);
#endif

    free(work);
    free(out);
    return zipscore;
}

/*-----------------------------------------------------------------------------*/  
/* calculate fitness */
void calculate_fitness(Pool *pool, int ix, Prots *prots, FreqSet *bg_freq, int alphabet_array_len, char *alphabet, int *polyfasta_count, int total_len, char *setfasta, float *ptr_KL_distance)
{
    int i = 0;
    int allocated = 1;
    int search_allocated = 1;
    char *polyfasta;
    char *search_polyfasta;                      /* search string with single sequence separeted by '-' char */
    polyfasta = safe_malloc(allocated);
    strcpy(polyfasta, "");
    search_polyfasta = safe_malloc(search_allocated);
    strcpy(search_polyfasta, "");

    for(i = 0; i < genenum; i++)
    {
#ifdef BIT
        if(get_bitgene(&pool[ix].bitgenome[i]) == 1)
#else
        if(pool[ix].genome[i] == 1)
#endif
        {
            /* JK: replaced by generic lines below
			if (strcmp(alphabet, "MV2000") == 0)
            {
                allocated += strlen(prots->protein_entries[i].aaseq);
                polyfasta = safe_realloc(polyfasta, allocated * sizeof(char));
                strcat(polyfasta, prots->protein_entries[i].aaseq);  

                search_allocated += strlen(prots->protein_entries[i].aaseq) + 1;
                search_polyfasta = safe_realloc(search_polyfasta, search_allocated * sizeof(char));
                if (strcmp(search_polyfasta, "") != 0)
                    strcat(search_polyfasta, "-");
                strcat(search_polyfasta, prots->protein_entries[i].aaseq);
            }
            if (strcmp(alphabet, "Camproux2004") == 0)
            {
                allocated += strlen(prots->protein_entries[i].SAseq);
                polyfasta = safe_realloc(polyfasta, allocated * sizeof(char));
                strcat(polyfasta, prots->protein_entries[i].SAseq);  

                search_allocated += strlen(prots->protein_entries[i].SAseq) + 1;
                search_polyfasta = safe_realloc(search_polyfasta, search_allocated * sizeof(char));
                if (strcmp(search_polyfasta, "") != 0)
                    strcat(search_polyfasta, "-");
                strcat(search_polyfasta, prots->protein_entries[i].SAseq);
            }
			*/

            allocated += strlen(prots->protein_entries[i].SAseq);
            polyfasta = safe_realloc(polyfasta, allocated * sizeof(char));
            strcat(polyfasta, prots->protein_entries[i].SAseq);  

            search_allocated += strlen(prots->protein_entries[i].SAseq) + 1;
            search_polyfasta = safe_realloc(search_polyfasta, search_allocated * sizeof(char));
            if (strcmp(search_polyfasta, "") != 0)
                strcat(search_polyfasta, "-");
            strcat(search_polyfasta, prots->protein_entries[i].SAseq);

#ifdef DEBUG
            fprintf(stderr,"|(%3d)==>%s<==|\n", i, polyfasta);
#endif
        }
    }
#ifndef COMPRESS_SCORE
    if (strcmp(polyfasta, "") == 0)     /* if empty selection the score is set to 0.0 to avoid segfault */
        pool[ix].fitness = 0.0;
    else
        pool[ix].fitness = score_seq(polyfasta, search_polyfasta, setfasta, bg_freq, alphabet_array_len, ptr_KL_distance); 
#endif
#ifdef COMPRESS_SCORE
    pool[ix].fitness = score_compress(polyfasta, strlen(polyfasta), total_len);
#endif

    get_counts(polyfasta, polyfasta_count, alphabet_array_len);

    free(polyfasta);
    free(search_polyfasta);
}

/*-----------------------------------------------------------------------------*/
/* print selected protein list */
void print_selected_proteins_lastrun(Pool *pool, Prots *prots, FILE *proteinoutfile, FreqSet *bg_freq, int alphabet_array_len, char *alphabet, int *polyfasta_count, int total_len, char
*selectedfilename, char *setfasta)
{
    unsigned int i, j, k;
    FILE *selectedout;
    float KL_distance = 1.0;

    selectedout = safe_open(selectedfilename, "w");

    fprintf(stdout, "\n");
    fprintf(proteinoutfile, "\n");

    fprintf(stdout, "     ");
    fprintf(proteinoutfile, "     ");

    for (j = 0; j < genenum; ++j)
    {
        fprintf(stdout, "%10s ", prots->protein_entries[j].protein_name);
        fprintf(proteinoutfile, "%10s ", prots->protein_entries[j].protein_name);
    }
    fprintf(proteinoutfile, "\tfitness KLdist\t");
    for(i = 0; i < alphabet_array_len; i++)
    {
        if(bg_freq->freq[i] == 0.0)
            continue;
        fprintf(proteinoutfile, "%4c", bg_freq->codeOrder[i]);
    }
    fprintf(stdout, "\n");
    fprintf(proteinoutfile, "\n");

    for (i = 0; i < popsize; ++i)
    {
        if (mem_fitness(&pool[0], i) == 0)
            calculate_fitness(&pool[0], i, prots, bg_freq, alphabet_array_len, alphabet, polyfasta_count, total_len, setfasta, &KL_distance);

        fprintf(stdout, "%4d: ", i);
        fprintf(proteinoutfile, "%4d: ", i);

        fprintf(selectedout, "%4d: ", i);

        for (j = 0; j < genenum; ++j)
        {
#ifdef BIT
            fprintf(stdout, "%10d ", get_bitgene(&pool[i].bitgenome[j]));
            fprintf(proteinoutfile, "%10d ", get_bitgene(&pool[i].bitgenome[j]));
#else
            fprintf(stdout, "%10d ", pool[i].genome[j]);
            fprintf(proteinoutfile, "%10d ", pool[i].genome[j]);
#endif
#ifdef BIT
            if (get_bitgene(&pool[i].bitgenome[j]) == 1)
#else
            if (pool[i].genome[j] == 1)
#endif
            {
				/* JK: should be replaced later
                if (strcmp(alphabet, "MV2000") == 0)
                    fprintf(selectedout, "%s+", prots->protein_entries[j].aaseq);
                else
				*/
                fprintf(selectedout, "%s+", prots->protein_entries[j].SAseq);
            }
        }
        fprintf(selectedout, "\n");

        fprintf(stdout, "\t%6.4f\t%6.4f\t", pool[i].fitness, KL_distance);
        fprintf(proteinoutfile, "\t%6.4f\t%6.4f\t", pool[i].fitness, KL_distance);
        for(k = 0; k < alphabet_array_len; k++)
        {
            if(bg_freq->freq[k] == 0.0)
                continue;
            fprintf(stdout, "%4d", polyfasta_count[k]);
            fprintf(proteinoutfile, "%4d", polyfasta_count[k]);
        }
        fprintf(stdout, "\n");
        fprintf(proteinoutfile, "\n");

    }
    fprintf(stdout, "\n");
    fflush(stdout);
    fprintf(proteinoutfile, "\n");
    fflush(proteinoutfile);

    fclose(selectedout);
}

/*-----------------------------------------------------------------------------*/  
/* parse the protein list file */
void parse_proteinlist(char *listfilename, Prots *prots)
{
    FILE *listfile;
    char line[256];
    unsigned int i;
    unsigned int k = 0;
    int allocated = 64;

    /*----------------------------------------------*/
    /* record the protein names in an array*/
    prots->protein_entries = safe_malloc(allocated * sizeof(ProteinEntry));

    listfile = safe_open(listfilename, "r");
    while(fgets(line, 256, listfile))
    { 
    prots->protein_entries[k].protein_name = safe_malloc(256 * sizeof(char));
    for (i = 0; i < 256 && isprint(line[i]) > 0; ++ i)
        prots->protein_entries[k].protein_name[i] = line[i];
    prots->protein_entries[k].protein_name[i] = '\0';

    if (i > 0)
        ++ k;
    else
        fprintf(stderr, "Warning: empty protein name at line %d\n", k);

    if (k == allocated)
    {
        allocated += 64;
        prots->protein_entries = safe_realloc(prots->protein_entries, allocated * sizeof(ProteinEntry));
    }
    }

    fclose(listfile);

    prots->protnum = k;
}

/*-----------------------------------------------------------------------------*/  
/* fill the protein table*/
char *fill_protein_table(Prots *prots, FreqSet *freqSet, char *alphabet, int *total_len)
{
	/*int i;*/
    int k = 0;
    int count_aa[26];
    int count_SA[33];
    char *fastafilename = 0;
    char *SAfastafilename = 0;
    FILE *tableout;
    FILE *fastafile;
    FILE *SAfastafile;
    int total_aalen = 0;
    int total_SAlen = 0;
    char *aasetfasta; /* aa sequence set */
    char *SAsetfasta; /* structure sequence set */
    FILE *setout;
    int nsymb = 0; /* local variable to satisfy word_entropy interface */
    float KL_distance = 0.;

    prots->aa_entropy_sum = 0.;
    prots->SA_entropy_sum = 0.;

    tableout = safe_open("protein_table.dat", "w");
    fprintf(tableout, "Basename  \taa_score.\tsa_score.\t");
	/* JK: removed for simplicity
    for(i = 0; i < 26; i++)
    {
        if(MV2000[i].freq == 0.0)
            continue;
        fprintf(tableout, "%4c", MV2000[i].aa);
    }
    fprintf(tableout, "\t\t");
    for(i = 0; i < 33; i++)
    {
        if(Camproux2004[i].freq == 0.0)
            continue;
        fprintf(tableout, "%4c", Camproux2004[i].aa);
    }
    fprintf(tableout, "\n");*/

    for(k = 0; k < prots->protnum; k ++)
    {
        fastafilename = safe_malloc(strlen(prots->protein_entries[k].protein_name) + 7);
		printf("%d %s\n", k, prots->protein_entries[k].protein_name); /* JK: debug */
        /*sprintf(fastafilename, "%s.AA.fasta", prots->protein_entries[k].protein_name);*/
        sprintf(fastafilename, "%s.fasta", prots->protein_entries[k].protein_name);/* JK: debug */
        fastafile = safe_open(fastafilename, "r");

        SAfastafilename = safe_malloc(strlen(prots->protein_entries[k].protein_name) + 10);
        sprintf(SAfastafilename, "%s.SA.fasta", prots->protein_entries[k].protein_name);
        SAfastafile = safe_open(SAfastafilename, "r");

        read_sequence(fastafile, SAfastafile, prots, k);

        fclose(fastafile);
		free(fastafilename);
        fclose(SAfastafile);
		free(SAfastafilename);

        total_aalen += strlen(prots->protein_entries[k].aaseq);
        total_SAlen += strlen(prots->protein_entries[k].SAseq);
    }

    aasetfasta = safe_malloc(sizeof(char) * (total_aalen + prots->protnum + 1));
    SAsetfasta = safe_malloc(sizeof(char) * (total_SAlen + prots->protnum + 1));

    strcpy(aasetfasta, prots->protein_entries[0].aaseq);
    strcpy(SAsetfasta, prots->protein_entries[0].SAseq);
        
    for(k = 1; k < prots->protnum; k++)
    {
        strcat(aasetfasta, "-");
        strcat(aasetfasta, prots->protein_entries[k].aaseq);
        strcat(SAsetfasta, "-");
        strcat(SAsetfasta, prots->protein_entries[k].SAseq);
    }

    for(k = 0; k < prots->protnum; k++)
    {
#ifndef COMPRESS_SCORE
        prots->protein_entries[k].aaseq_score = score_seq(prots->protein_entries[k].aaseq, prots->protein_entries[k].aaseq, aasetfasta, freqSet, 26, &KL_distance);
        prots->protein_entries[k].aaseq_entropy = word_entropy(prots->protein_entries[k].aaseq, prots->protein_entries[k].aaseq, &nsymb);
        prots->protein_entries[k].SAseq_score = score_seq(prots->protein_entries[k].SAseq, prots->protein_entries[k].SAseq, SAsetfasta, freqSet, 33, &KL_distance);
        prots->protein_entries[k].SAseq_entropy = word_entropy(prots->protein_entries[k].SAseq, prots->protein_entries[k].SAseq, &nsymb);
        prots->aa_entropy_sum += prots->protein_entries[k].aaseq_entropy;
        prots->SA_entropy_sum += prots->protein_entries[k].SAseq_entropy;
#endif
#ifdef COMPRESS_SCORE
        prots->protein_entries[k].aaseq_score = score_compress(prots->protein_entries[k].aaseq, strlen(prots->protein_entries[k].aaseq), strlen(prots->protein_entries[k].aaseq));
        prots->protein_entries[k].SAseq_score = score_compress(prots->protein_entries[k].SAseq, strlen(prots->protein_entries[k].SAseq), strlen(prots->protein_entries[k].SAseq));
#endif

        get_counts(prots->protein_entries[k].aaseq, count_aa, 26);
        get_counts(prots->protein_entries[k].SAseq, count_SA, 33);

        fprintf(tableout, "%10s\t%8.3f\t%8.3f\t", 
            prots->protein_entries[k].protein_name,
            prots->protein_entries[k].aaseq_score, prots->protein_entries[k].SAseq_score);

		/* JK: removed for simplicity
        for(i = 0; i < 26; i++)
        {
            if(MV2000[i].freq == 0.0)
                continue;
            fprintf(tableout, "%4d", count_aa[i]);
        }
        fprintf(tableout, "\t\t");
        for(i = 0; i < 33; i++)
        {
            if(Camproux2004[i].freq == 0.0)
                continue;
            fprintf(tableout, "%4d", count_SA[i]);
        }*/
        fprintf(tableout, "\t%s\t%s\n", 
            prots->protein_entries[k].aaseq, prots->protein_entries[k].SAseq);
        fflush(tableout);
    }
    fclose(tableout);

    setout = safe_open("originalset.seq", "w");

	/* JK: should be replaced later
    if (strcmp(alphabet, "MV2000") == 0)
    {
        fprintf(setout, "%s", aasetfasta);
        fclose(setout);
        *total_len = total_aalen;
		free(SAsetfasta);
        return aasetfasta;
    }
    else 
    { */
        fprintf(setout, "%s", SAsetfasta);
        fclose(setout);
        *total_len = total_SAlen;
		free(aasetfasta);
        return SAsetfasta;
   /* }*/
}

/*-----------------------------------------------------------------------------*/  
/* generate random background distribution */
void generate_background(Pool *pool, int ix, Prots *prots, char *alphabet, Pool *bg_pool)
{
    int numselected_prot = 0;
    int selected_len = 0;
    int random_len;
    int set_lengths[genenum];
    int i,j;
    FILE *lengthlog;

    /* read the length and composition of the best individual */
    for(j = 0; j < genenum; j++)
    {
#ifdef BIT
        if(get_bitgene(&pool[ix].bitgenome[j]) == 1) 
#else
        if(pool[ix].genome[j] == 1) 
#endif
        {
            ++ numselected_prot;
			/* JK : should be replaced later
            if (strcmp(alphabet, "MV2000") == 0)
                selected_len += strlen(prots->protein_entries[j].aaseq);
            else */
                selected_len += strlen(prots->protein_entries[j].SAseq);
        }
		/* JK : should be replaced later
        if (strcmp(alphabet, "MV2000") == 0)
            set_lengths[j] = strlen(prots->protein_entries[j].aaseq);
        else */
            set_lengths[j] = strlen(prots->protein_entries[j].SAseq);
    }

    /* length occurencies are recorded in a file */
    lengthlog = safe_open("background.len","w");

    /* loop over the number of desired random individuals */
    for(i = 0; i < numran_pool; i++)
    {
        random_len = 0;
        bg_pool[i].fitness = 0; /* initialise fitness */

        /* initialize genomes at random */
        for (j = 0; j < genenum; ++j) /* genes */
        {
#ifdef BIT
            set_bitgene(&bg_pool[i].bitgenome[j], ((rand()/(double)INT_MAX) < ((double) numselected_prot/genenum)) ? 1 : 0);
            random_len += set_lengths[j] * get_bitgene(&bg_pool[i].bitgenome[j]);
#else
            bg_pool[i].genome[j] = ((rand()/(double)INT_MAX) < ((double) numselected_prot/genenum)) ? 1 : 0;
            random_len += set_lengths[j] * bg_pool[i].genome[j];
#endif
        }

        /* check that the individual has length similar to the best individual */
        while(abs(random_len - selected_len) > (selected_len * 0.05)) 
        {
            /* if too short: randomly enlongate it */
            if (random_len < selected_len)
            {
#ifdef BIT
                do{
                    j = (int) (rand()/(double)INT_MAX * genenum);
                }while(get_bitgene(&bg_pool[i].bitgenome[j]) == 1);
                set_bitgene(&bg_pool[i].bitgenome[j], 1);
                random_len += set_lengths[j];
#else
                do{
                    j = (int) (rand()/(double)INT_MAX * genenum);
                }while(bg_pool[i].genome[j] == 1);
                bg_pool[i].genome[j] = 1;
                random_len += set_lengths[j];
#endif
            }
            /* if too long: randomly shorten it */
            if (random_len > selected_len)
            {
#ifdef BIT
                do{
                    j = (int) (rand()/(double)INT_MAX * genenum);
                }while(get_bitgene(&bg_pool[i].bitgenome[j]) == 0);
                set_bitgene(&bg_pool[i].bitgenome[j], 0);
                random_len -= set_lengths[j];
#else                
                do{
                    j = (int) (rand()/(double)INT_MAX * genenum);
                }while(bg_pool[i].genome[j] == 0);
                bg_pool[i].genome[j] = 0;
                random_len -= set_lengths[j];
#endif
            }
        }
        fprintf(lengthlog, "%3d: %3d %3d %5.3f\n", i, random_len, selected_len, (float) random_len/selected_len);
        fflush(lengthlog);
    }
    fclose(lengthlog);

}

/*-----------------------------------------------------------------------------*/  
/* parse the command line arguments*/
void parse_args(int argc, char **argv, char *listfilename, char *alphabet, float *ptr_minset_perc, char *freqfilename)
{
        int c;
        opterr = 0; 
		extern int k_word_len;
		extern int freqset;

        while ((c = getopt (argc, argv, "hf:a:p:k:c:")) != -1) 
                switch (c)
                {
                        case 'f': 
                                strcpy(listfilename,optarg); 
                                break; 
                        case 'a': 
                                strcpy(alphabet,optarg); 
                                break; 
                        case 'p': 
                                *ptr_minset_perc = atof(optarg);
                                break; 
                        case 'k': 
                                k_word_len = atoi(optarg);
                                break; 
                        case 'c': 
                                strcpy(freqfilename, optarg);
								++ freqset;
                                break; 
                        case 'h': 
                                fprintf(stderr, "Usage: ./ga -f listfilename -a alphabet -p minset_perc <-k k_word-length> <-c freqfilename>\n");
                                exit(0);
                        case '?': 
                                if (isprint (optopt)) 
                                        fprintf (stderr, "Unknown option `-%c'.\n", optopt); 
                                else 
                                        fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt); 
                        default:
                                fprintf(stderr, "Usage: ./ga -f listfilename -a alphabet -p minset_perc <-k k_word-length> <-c freqfilename>\n");
                                exit(0);

                }
}

/*_____________________________________________________________________________*/
/** load the frequency set data into freqSet */
void init_freq_data(FreqSet *freqSet, char *setName, char *freqFileName) {

	extern int freqset; /*flag indicating uploaded frequency set */
	FILE *freqFile; /* file containing uploaded frequency set */

    if (freqset) {
        fprintf(stderr, "code character frequency file: %s\n", freqFileName);
        freqFile = safe_open(freqFileName, "r");
        read_frequencies(freqFile, freqSet);
        fclose(freqFile);
    } else {
        init_freq_set(freqSet, setName);
    }
}

/*-----------------------------------------------------------------------------*/
/* main function */
int main(int argc, char **argv)
{
    char listfilename[200];
    char freqfilename[200];
    /*char alphabet[200];*/
    float minset_perc = 0.0;
    int nselected = 0;

    /*----------------------------------------------------------------------------*/
    /* JK: replaced by routines 'freq.[ch]', 'freq_data.h' and 'getfreqs.[ch]'
	AminoacidFreq MV2000[26];
    AminoacidFreq Camproux2004[33];*/

    /*----------------------------------------------------------------------------*/

    /*AminoacidFreq *bg_freq; JK: replaced by 'freqs' routines */
	FreqSet *freqSet; /* character frequencies of specified alphabet */
	
    /*int alphabet_array_len;*/
    int *polyfasta_count;
    int total_len;
    char *psetfasta; /* pointer fo sequence set */
    float KL_distance = 1.0;

#ifdef HOMSTRAD
    FILE *homstradlist;
    char *homstradlistname = "homstrad.dist.list";
    char *homstrad = "/home/jkleinj/data/homstrad_1104/";
#endif
#ifdef CATH
    FILE *cathlist;
    char *cathlistname = "cath.pairs.test.list";
    char *cath = "/mb/databases/cath/v2.4";
#endif

    FILE *gaoutfile;
    /*JK FILE *pooloutfile;*/
    FILE *proteinoutfile;
	FILE *messagefile;
	char messagefilename[12] = "message.php";
    /* digit encoding of 'pool' output files: */
    /* 1.: this repat, 2. repga, 3. this fraction, 4. jackknife */
    char gaoutfilename[13] = "0_0.0_0.ga";
    /*JK char pooloutfilename[13] = "0_0.0_0.pool";*/
    char proteinoutfilename[13] = "0_0.0_0.prot";
    char selectedsetfilename[13] = "0_0.0_0.seq";
#ifdef ALIGNMODE
    char pdb0filename[200] = ""; /* RUN input file */
    char pdb1filename[200] = ""; /* RUN input file */
    char s[50];
#endif

    Prots prots; /* list of proteins */

    /*JK char cleanupcmnd[500];*/ /* command string to clean up */

    int f, i, l, ix, iy, x, z, j;

#ifdef ALIGNMODE
    Parameter parameter[parnum]; /* RUN parameters */
    genenum = parnum;
#endif
    Fam fam; /* family list */

    Pool pool[popsize]; /* gene pool */

#ifdef GENERATE_BACKGROUND
    Pool bg_pool[numran_pool]; /* random gene pool */
    FILE *bgoutfile;
#endif

    int average[popsize]; /* average values for equilibrium */

	clock_t start, end; /* measure process time */

    /*----------------------------------------------------------------------------*/
    /* JK: replaced by 'freq' routines
     set_aa_seq_freq(&MV2000[0]);
     set_aa_str_freq(&Camproux2004[0]); */
	freqSet = safe_malloc(sizeof(FreqSet));
	freqSet->setname = safe_malloc(64 * sizeof(char));

    /*----------------------------------------------------------------------------*/
    /* read the fasta file of the protein collection */
    parse_args(argc, argv, listfilename, freqSet->setname, &minset_perc, freqfilename); /* parse the command line argument */
    if (strlen(listfilename) == 0)
    {
		fprintf(stderr, "Usage: ./ga -f listfilename -a alphabet -p minset_perc <-k k_word-length> <-c freqfilename>\n");
		exit(0);
    }
    if (minset_perc <= 0 || minset_perc > 100)
    {
		fprintf(stderr, "Usage: ./ga -f listfilename -a alphabet -p minset_perc <-k k_word-length> <-c freqfilename>\n");
		exit(0);
    }

    /* JK: replaced by 'freq' routines
	if (strcmp(alphabet, "MV2000") == 0)
    {
        bg_freq = MV2000;
        alphabet_array_len = 26;
    } else if (strcmp(alphabet, "Camproux2004") == 0)
    {
        bg_freq = Camproux2004;
        alphabet_array_len = 33;
    } else
    {
        fprintf(stderr, "Alphabet %s not implemented\n", alphabet);
        fprintf(stderr, "   available Alphabet:\n");
        fprintf(stderr, "   MV2000\n");
        fprintf(stderr, "   Camproux2004\n");
        exit(99);
    }
	*/

	init_freq_data(freqSet, &(freqSet->setname[0]), &(freqfilename[0]));

    /*polyfasta_count = safe_malloc(sizeof(int) * alphabet_array_len);*/
    polyfasta_count = safe_malloc(sizeof(int) * freqSet->codeRange);

    parse_proteinlist(listfilename, &prots);

    psetfasta = fill_protein_table(&prots, freqSet, &(freqSet->setname[0]), &total_len);

    genenum = prots.protnum;
    nselected = floor(genenum * minset_perc/100);

    for(i = 0; i < popsize; ++i)
	{
#ifdef BIT
		pool[i].bitgenome = safe_malloc(sizeof(Bitgene) * genenum); /* bitgene pool */
#else
        pool[i].genome = safe_malloc(sizeof(int) * genenum); /* gene pool */
#endif
	}

#ifdef GENERATE_BACKGROUND
    for(i = 0; i < numran_pool; ++i)
    {
#ifdef BIT
        bg_pool[i].bitgenome = safe_malloc(sizeof(int) * genenum); /* gene pool */
        for (j = 0; j < genenum; ++j) /* genes */
            set_bitgene(&bg_pool[i].bitgenome[j], 0); /* initialize to zero */
#else
        bg_pool[i].genome = safe_malloc(sizeof(int) * genenum); /* gene pool */
        for (j = 0; j < genenum; ++j) /* genes */
            bg_pool[i].genome[j] = 0; /* initialize to zero */
#endif /* (BIT)*/
    }
#endif /* (GENERATE_BACKGROUND) */

    /*----------------------------------------------------------------------------*/
    /* repeat entire GA */
    for (z = 0;  z < repga; ++z)
    {
        /* set the starting point for pseudorandom integers */
        srand((unsigned int)(time(NULL)%10000));

        /* repeat for each jackknife fraction */
        for (x = 0;  x < jackknife; ++x)
        {
            /*----------------------------------------------------------------------------*/
            init_random_pool(&pool[0]); /* initialise random pool and fitness */
            /*init_seeded_pool(&pool[0]);*/ /* initialise seeded pool and fitness */
            /*print_pool(&pool[0], outfilename);*/ /* print out pool */
            /*print_pool(&pool[0], outfilename);*/ /* print out pool */

            /* at the start constrain the newly generated pool... */
            for (iy = 0; iy < popsize; ++iy)
                constrain_genome(&pool[0], iy, nselected);

            /*----------------------------------------------------------------------------*/
            /* output file */
            gaoutfilename[0] = z + '0';
            gaoutfilename[2] = repga + '0';
            gaoutfilename[4] = x + '0';
            gaoutfilename[6] = jackknife + '0';
            /*JK pooloutfilename[0] = z + '0';
            pooloutfilename[2] = repga + '0';
            pooloutfilename[4] = x + '0';
            pooloutfilename[6] = jackknife + '0';*/
            proteinoutfilename[0] = z + '0';
            proteinoutfilename[2] = repga + '0';
            proteinoutfilename[4] = x + '0';
            proteinoutfilename[6] = jackknife + '0';
            selectedsetfilename[0] = z + '0';
            selectedsetfilename[2] = repga + '0';
            selectedsetfilename[4] = x + '0';
            selectedsetfilename[6] = jackknife + '0';
            /*JK sprintf(cleanupcmnd, "rm %s %s %s %s", gaoutfilename, pooloutfilename, proteinoutfilename, selectedsetfilename);
            system(cleanupcmnd);*/
            gaoutfile = safe_open(gaoutfilename, "a");
            /*JK pooloutfile = safe_open(pooloutfilename, "ab");*/
            proteinoutfile = safe_open(proteinoutfilename, "a");
            /*JK fwrite(&genenum, sizeof(int), 1, pooloutfile);
            fwrite(&popsize, sizeof(int), 1, pooloutfile);
            fwrite(&loop, sizeof(int), 1, pooloutfile);
            fflush(pooloutfile);*/

			start = clock();
            /*----------------------------------------------------------------------------*/
            for (l = 0, ix = 0; l < loop; ++l, ix = fitmate)
            {
				if (((z+1) * (x+1) * (l+1)) % 100 == 0)
				{
					end = clock();
					messagefile = safe_open(messagefilename, "w");
					fprintf(messagefile, "process %d / %d\n, estimated runtime %11.1f s",
						((z+1) * (x+1) * (l+1)), (repga * jackknife * loop), 
						(float)(repga * jackknife * loop) * ((float)end - (float)start) / (CLOCKS_PER_SEC * (float)((z+1) * (x+1) * (l+1))));
					fclose(messagefile);
				}

                for (pool[ix].fitness = 0; ix < popsize; ++ix)
                {
#ifdef ALIGNMODE
                    fprintf(stderr, "repeat %d/%d, db-frac %d/%d, generation %d/%d, genome %d/%d\n",
                        z, repga, x, jackknife, l, loop, ix, popsize);
                    fflush(stderr);
                    fprintf(gaoutfile, "repeat %d/%d, db-frac %d/%d, generation %d/%d, genome %d/%d\n",
                        z, repga, x, jackknife, l, loop, ix, popsize);
                    fflush(gaoutfile);
#endif

                    /*----------------------------------------------------------------------------*/
                    /* run benchmark */

                    /* transform genes to RUN parameters */
#ifdef ALIGNMODE
                    set_parameters(&pool[0], ix, &parameter[0], gaoutfile);
#endif
#ifndef ALIGNMODE
                    fam.nfam = 1;
#endif

                    for (f = x; f < fam.nfam; f += jackknife)
                    {
                        /* file names */
#ifdef HOMSTRAD
                        sprintf(pdb0filename, "%s%s/%s.0.pdb",
                            homstrad, fam.family[f], fam.family[f]);
                        sprintf(pdb1filename, "%s%s/%s.1.pdb",
                            homstrad, fam.family[f], fam.family[f]);
#endif
#ifdef CATH
                        sprintf(pdb0filename, "%s/%s",
                            cath, fam.family[f], fam.family[f]);
                        ++ f;
                        sprintf(pdb1filename, "%s/%s",
                            cath, fam.family[f]);
#endif

                        /* function 'mem_fitness' assigns already known fitness value */
                        /* otherwise calculate fitness  */
                        if (mem_fitness(&pool[0], ix) == 0)
                        {
#ifdef ALIGNMODE
                            /*pool[ix].fitness += aliss(pdb0filename, pdb1filename, &parameter[0]);*/
                            /*pool[ix].fitness += 0;*/
#endif
#ifndef ALIGNMODE
                            calculate_fitness(&pool[0], ix, &prots, freqSet, freqSet->codeRange, &(freqSet->setname[0]), polyfasta_count, total_len, psetfasta, &KL_distance);
#endif
                        }
                    }
                    fprintf(stderr, "fitness=%f\n", pool[ix].fitness);
                    fflush(stdout);
#ifdef ALIGNMODE
                    fprintf(gaoutfile, "fitness=%f\n", pool[ix].fitness);
                    fflush(gaoutfile);
#endif

                }

                sort_fitness(&pool[0]); /* sort pool by fitness */
#ifdef ALIGNMODE
                /*JK print_pool_ascii(&pool[0], pooloutfile);*/ /* print out pool */
#endif
#ifndef ALIGNMODE
                /*JK print_pool_bin(&pool[0], pooloutfile);*/ /* print out pool */
#endif
                if (l == loop-1)
                {
                    print_selected_proteins_lastrun(pool, &prots, proteinoutfile, freqSet, freqSet->codeRange, &(freqSet->setname[0]), polyfasta_count, total_len, selectedsetfilename, psetfasta); /* print selected protein list */
#ifdef GENERATE_BACKGROUND
                    generate_background(pool, 0, &prots, &(freqSet->setname[0]), bg_pool);

                    bgoutfile = safe_open("background_dist.prot", "w");
                    print_selected_proteins_lastrun(bg_pool, &prots, bgoutfile, freqSet, freqSet->codeRange, &(freqSet->setname[0]), polyfasta_count, total_len, "background.seq", psetfasta);
                    fclose(bgoutfile);
#endif
                }
                
                for (iy = fitmate; iy < popsize; ++iy)
                {
                    if (breedmode == 0)
                    {
                        breed_crossover(&pool[0], iy); /* crossover fittest genomes */
#ifndef ALIGNMODE
                        constrain_genome(&pool[0], iy, nselected);
#endif
                    }
                    if (breedmode == 1)
                    {
                        equilibrium(&pool[0], &average[0]); /* create equilibrium values */
                        breed_equilibrium(&pool[0], &average[0], iy);
                    }
                }
            }
            fclose(gaoutfile);
            /*fclose(pooloutfile);*/
			fclose(proteinoutfile);
        }
    }

    /* free memory */
    for(i = 0; i < popsize; ++i)
	{
#ifdef BIT
        free(pool[i].bitgenome);
#else
        free(pool[i].genome);
#endif
	}
#ifdef ALIGNMODE
    for (i = 0; i < fam.nfam; ++ i)
        free(fam.family[i]);
#endif
    for(i = 0; i < prots.protnum; i++)
    {
        free(prots.protein_entries[i].protein_name);
        free(prots.protein_entries[i].aa_description);
        free(prots.protein_entries[i].aaseq);
        free(prots.protein_entries[i].SA_description);
        free(prots.protein_entries[i].SAseq);
	}
    free(prots.protein_entries);
    free(psetfasta);
#ifdef GENERATE_BACKGROUND
    for(i = 0; i < numran_pool; ++i)
	{
#ifdef BIT
        free(&bg_pool[i].bitgenome);
#else
        free(bg_pool[i].genome);
#endif /* (BIT) */
	}
#endif /* (GENERATE_BACKGROUND) */

	free(polyfasta_count);

	free(freqSet->freq);
	free(freqSet);

    return 0;
}
                                                                                   
/*-----------------------------------------------------------------------------*/  
