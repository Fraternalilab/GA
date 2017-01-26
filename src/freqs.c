/*=============================================================================
freqs.c : code character frequency definitions
(C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "freqs.h"
#include "freqs_data.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* initialise fragment sets with constant values for the available alphabets */
int init_freq_set(FreqSet *selected_set, char *setname){

    int i;

    ConstantFreqSet *constant_set;

    for (constant_set = constant_freq_sets; constant_set < constant_freq_sets + sizeof constant_freq_sets / sizeof *constant_freq_sets; constant_set ++) 
        if (strcmp(constant_set->setname, setname) == 0){
            selected_set->setname = constant_set->setname;
            selected_set->nFreq = constant_set->nFreq;
            selected_set->codeRange = constant_set->codeRange;
            selected_set->codeOrder = constant_set->codeOrder;

			selected_set->freq = safe_malloc(selected_set->nFreq * sizeof(float));
            for (i = 0; i < selected_set->nFreq; ++ i)
				selected_set->freq[i] = constant_set->freq[i];
            return 0;
        }
    fprintf(stderr, "\nExiting: Frequency set '%s' not implemented!\nAvailable sets:\n", setname);

    for (constant_set = constant_freq_sets; constant_set < constant_freq_sets + sizeof constant_freq_sets / sizeof *constant_freq_sets; constant_set ++)
        fprintf(stderr, "\t%s\n", constant_set->setname);

    exit(1);
}

