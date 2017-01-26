/*===============================================================================
freqs.h : code character frequency definitions
(C) 2007 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef FREQS_H
#define FREQS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*____________________________________________________________________________*/
/* structures */

/* Frequency set */
typedef struct {
    char *setname;
    int nFreq;
	int codeRange;
    char *codeOrder;
    float *freq;
} FreqSet;

/*____________________________________________________________________________*/
/* prototypes */
int init_freq_set(FreqSet *selected_set, char *setname);

#endif

