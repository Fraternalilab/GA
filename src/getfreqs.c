/*==============================================================================
getfreqs.c : Read code character frequencies
(C) 2007 and Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "getfreqs.h"
#include "safe.h"

#define MAX(a,b)	a > b ? a : b
#define MIN(a,b)	a < b ? a : b

/* reading the code characters ('char', first column) and 
	their frequencies ('float', second column);
	the two columns are separated by white space;
	the characters define the entire code alphabet;
	characters are expected to be within [a-zA-Z];
	frequencies are expected to sum up to 1, otherwise a warning is issued */
void read_frequencies(FILE *freqFile, FreqSet *freqSet)
{
	unsigned int nLine = 0; /* number of encountered lines */
	unsigned int allocated = 64;
	unsigned int rangeMax = -1; /* ASCII extrema of code characters */
	unsigned int rangeMin = 1000;
	float freqSum = 0; /* sum of frequency values */
	freqSet->codeOrder = safe_malloc(allocated * sizeof(char));
	freqSet->freq = safe_malloc(allocated * sizeof(float));

	while(! feof(freqFile)) {
		/* strict format checking */
		assert (fscanf(freqFile, "%c%f\n",
		&(freqSet->codeOrder[nLine]),
		&(freqSet->freq[nLine])) == 2);

		/* determine ASCII range of code characters */
		assert((freqSet->codeOrder[nLine] - 'A') >= 0);
		rangeMax = MAX(rangeMax, (freqSet->codeOrder[nLine] - 'A'));
		rangeMin = MIN(rangeMin, (freqSet->codeOrder[nLine] - 'A'));

		freqSum += freqSet->freq[nLine];
		++ nLine;

		if (nLine == allocated) {
			allocated += 64;
			freqSet->codeOrder = safe_realloc(freqSet->codeOrder, (allocated + 1) * sizeof(char));
			freqSet->freq = safe_realloc(freqSet->freq, allocated * sizeof(float));
		}
	}

	freqSet->codeOrder[nLine] = '\0';
    freqSet->nFreq = nLine;
	freqSet->codeRange = rangeMax - rangeMin + 1;
	if (freqSum != 1.)
	{
		fprintf(stderr, "Error: Frequencies should add up to 1, but actual value is %f\n", freqSum);
		exit(1);
	}
}

