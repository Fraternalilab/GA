/*==============================================================================
alignment.c : Read alignment of FASTA sequences
(C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "ga.h"
#include "safe.h"

/*-----------------------------------------------------------------------------*/
static int read_sequence_name(FILE *aafile, FILE *SAfile, ProteinEntry *protein_entries)
{
    /* Reads a line starting with '>' and followed by the sequence name */
    /* Leading white space is skipped.  Reading proceeds until the end  */
    /* of the line is reached.  The name is read into sequence->name.   */

    int ch, length = 0, allocated = 64;

    while ((ch = getc(aafile)) != EOF && isspace(ch))
	;

    if (ch != '>')
		return 0;

    protein_entries->aa_description = safe_malloc(allocated);

    do
	{
		protein_entries->aa_description[length ++] = ch;

		if (length == allocated)
			protein_entries->aa_description = safe_realloc(protein_entries->aa_description, allocated += 64);
    } while ((ch = getc(aafile)) != EOF && ch != '\n' && isprint(ch));

    protein_entries->aa_description[length] = '\0';

	/*-----------------------------------------------------------------------------*/
	ch = 0; length = 0; allocated = 64;
    while ((ch = getc(SAfile)) != EOF && isspace(ch))
	;

    if (ch != '>')
		return 0;

    protein_entries->SA_description = safe_malloc(allocated);

    do
	{
		protein_entries->SA_description[length ++] = ch;

		if (length == allocated)
			protein_entries->SA_description = safe_realloc(protein_entries->SA_description, allocated += 64);
    } while ((ch = getc(SAfile)) != EOF && ch != '\n' && isprint(ch));

    protein_entries->SA_description[length] = '\0';

	return 1;
}

/*-----------------------------------------------------------------------------*/
static void read_sequence_residues(FILE *aafile, FILE *SAfile, ProteinEntry *protein_entries)
{
    /* Reads the residues in a sequence, up to (but not including) the      */
    /* next sequence header (starting with '>'), or up to end of file.      */
    /* Residues may span multiple lines.  White space and gaps are skipped. */
    /* Nonalpha characters are rejected, resulting in an error message.     */
    /* Alpha characters are NOT converted to upper case.  The string is read*/
    /* into sequence->residues and zero-terminated; the length is stored    */
    /* into sequence->length.                                               */

    int ch, length = 0, allocated = 64;

    protein_entries->aaseq = safe_malloc(allocated);

    while ((ch = getc(aafile)) != EOF && ch != '>')
		if (isalpha(ch))
		{
			protein_entries->aaseq[length ++] = toupper(ch);

			if (length == allocated)
				protein_entries->aaseq = safe_realloc(protein_entries->aaseq, allocated += 64);
		}
		else if (ch != '.' && !isspace(ch))
		{
			fprintf(stderr, "illegal character '%c' in protein sequence\n", ch);
			exit(1);
		}

    if (ch == '>')
		ungetc('>', aafile);

    if (length == 0)
	{
		fprintf(stderr, "zero-sized sequence\n");
		exit(1);
    }

    protein_entries->aaseq[length] = '\0';

	/*-----------------------------------------------------------------------------*/
	ch = 0; length = 0; allocated = 64;
    protein_entries->SAseq = safe_malloc(allocated);

    while ((ch = getc(SAfile)) != EOF && ch != '>')
		if (isalpha(ch))
		{
			/* SA alphabet contains case-sensitive coding*/
			/*protein_entries->aaseq[length ++] = toupper(ch);*/ 
			protein_entries->SAseq[length ++] = ch;

			if (length == allocated)
				protein_entries->SAseq = safe_realloc(protein_entries->SAseq, allocated += 64);
		}
		else if (ch != '.' && !isspace(ch))
		{
			fprintf(stderr, "illegal character '%c' in protein sequence\n", ch);
			exit(1);
		}

    if (ch == '>')
		ungetc('>', SAfile);

    if (length == 0)
	{
		fprintf(stderr, "zero-sized sequence\n");
		exit(1);
    }

    protein_entries->SAseq[length] = '\0';
}

/*-----------------------------------------------------------------------------*/
int read_sequence(FILE *aafile, FILE *SAfile, Prots *prots, int k)
{
    if (read_sequence_name(aafile, SAfile, &prots->protein_entries[k]))
	{
		read_sequence_residues(aafile, SAfile, &prots->protein_entries[k]);
		return 1;
    }
	else
		return 0;
}

