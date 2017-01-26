/*==============================================================================
read_poolbin.c 
Copyright (C) 2006 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/*-----------------------------------------------------------------------------*/  
/* parse the command line arguments*/
void parse_args(int argc, char **argv, char *filename)
{
		int c;
		opterr = 0; 

		while ((c = getopt (argc, argv, "hf:")) != -1) 
				switch (c)
				{
						case 'f': 
								strcpy(filename,optarg); 
								break; 
						case 'h': 
								fprintf(stdout, "Usage: ./read_poolbin -f listfilename\n");
								exit(99);
						case '?': 
								if (isprint (optopt)) 
										fprintf (stdout, "Unknown option `-%c'.\n", optopt); 
								else 
										fprintf (stdout, "Unknown option character `\\x%x'.\n", optopt); 
						default:
								fprintf(stdout, "Usage: ./read_poolbin -f listfilename\n");
								exit(99);

				}
}

/*-----------------------------------------------------------------------------*/
/* main function */

int main(int argc, char **argv)
{
	FILE *file;
	char filename[200];
	int in, i, j, k;
	char ci;
	float fin;
	int genenum;
	int popsize;
	int loop;

	parse_args(argc, argv, filename);

	if (strlen(filename) == 0)
	{
			fprintf(stdout, "Usage: ./read_poolbin -f listfilename\n");
			exit(1);
	}

	file=fopen(filename,"rb");

	fread(&genenum, sizeof(int), 1, file);
	fread(&popsize, sizeof(int), 1, file);
	fread(&loop, sizeof(int), 1, file);

	fprintf(stdout, "genenum: %d\n", genenum);
	fprintf(stdout, "popsize: %d\n", popsize);
	fprintf(stdout, "loop:    %d\n", loop);

	for (i = 0; i < loop; ++i)
	{
		for (j = 0; j < popsize; ++j)
		{
			fread(&in, sizeof(int), 1, file);
			fprintf(stdout, "%4d: ", in);
			for (k = 0; k < genenum; ++k)
			{
				fread(&ci, sizeof(char), 1, file);
				fprintf(stdout, "%c", ci);
			}
			fread(&fin, sizeof(float), 1, file);
			fprintf(stdout, "\t%f", fin);
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	}

	fclose(file);
	return 0;
}
