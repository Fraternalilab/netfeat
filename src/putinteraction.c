/*==============================================================================
putinteraction.c : interaction output as list or matrix
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "putinteraction.h"

/*____________________________________________________________________________*/
/** print interaction matrix in ASCII format */
void print_interaction_matrix_ascii(char *outFileName, Ints *ints, Arg *arg)
{
	FILE *outFile;
	unsigned int i, j;

	if (! arg->silent) fprintf(stdout, "\t\tPrinting interaction matrix\n");

	outFile = safe_open(outFileName, "w");

	/* print number of nodes to first line */
	fprintf(outFile, "%d\n", ints->N);
	/* print matrix values*/
	/* columns */
	for (i = 0; i < ints->N; ++ i) {
		if (i > 0)
			fprintf(outFile, "\n");
		/* progress info (for large files) */

		/*if (i % 1000 == 0)
			fprintf(stderr, "%d / %d\n", i, ints->N);*/

		/* rows */
		for (j = 0; j < ints->N; ++ j)
			fprintf(outFile, "%1d ", ints->c[i][j]);
	}
	fprintf(outFile, "\n");
	fclose(outFile);
}

/*____________________________________________________________________________*/
/** print interaction matrix in binary format */
void print_interaction_matrix_binary(char *outFileName, Ints *ints)
{
	FILE *outFile;
	int matSize = ints->N * ints->N * sizeof(int);

	if ((outFile = fopen(outFileName, "wb")) == 0 ||
		 fwrite(ints->c, matSize, 1, outFile) != 1)
			Warning("Could not write interaction matrix file\n");

	fclose(outFile);
}


/*____________________________________________________________________________*/
/** print interaction list in ASCII format */
void print_interaction_list_ascii(char *outFileName, Ints *ints, Arg *arg)
{
	FILE *outFile;
	unsigned int i, j;

	if (! arg->silent) fprintf(stdout, "\t\tPrinting interaction list\n");

	outFile = safe_open(outFileName, "w");

	/* print matrix values*/
	/* columns */
	for (i = 0; i < (ints->N - 1); ++ i)
		/* rows */
		for (j = i + 1; j < ints->N; ++ j)
			/* interactions */
			if (ints->c[i][j])
				fprintf(outFile, "%8d\tpp\t%8d\n", i, j);
	
	fclose(outFile);
}

