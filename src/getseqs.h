/*===============================================================================
alignment.c : Read alignment of FASTA sequences
(C) 2004 Jens KLeinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

int read_sequence(FILE *aafile, FILE *SAfile, Prots *prots, int k);

#endif
