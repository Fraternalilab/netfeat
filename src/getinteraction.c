/*==============================================================================
getinteraction.c : interaction input as list or matrix
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "getinteraction.h"

/*____________________________________________________________________________*/
/** read node list */
void read_node_list(char *protFileName, Ints *ints, Arg *arg)
{
	FILE *protFile;
    int nLine = 0;
	int allocated = 64;

	if (! arg->silent) fprintf(stdout, "\tReading node list %s\n",
							protFileName);

	/* allocate node name array */
	ints->node = safe_malloc(allocated * sizeof(char [64]));
	protFile = safe_open(protFileName, "r");

    while(fscanf(protFile, "%s\n", &(ints->node[nLine][0])) == 1) { 
		/*fprintf(stderr, "%s:%d: %d %s\n", __FILE__, __LINE__, nLine, &(ints->node[nLine][0]));*/
        ++ nLine;

		if (nLine == allocated) {
			allocated += 64;
			ints->node = safe_realloc(ints->node, allocated * sizeof(char [64]));
		}
    }   
	fclose(protFile);
	ints->N = nLine;
	assert(ints->N > 1);

	if (! arg->silent)
		fprintf(stdout, "\t\tnumber of nodes N = %d\n", ints->N);
}

/*____________________________________________________________________________*/
/** read interaction list */
void read_interaction_list(char *intsFileName, Ints *ints, Arg *arg)
{
	FILE *intsFile;
    int nLine = 0;
	int allocated = 64;
	char dummi[8], dummy[8], dummj[8];

	if (! arg->silent) fprintf(stdout, "\tReading interaction list %s\n",
							intsFileName);

	ints->nInteractionListed = 0;
	ints->protein1 = safe_malloc(allocated * sizeof(char [64]));
	ints->protein2 = safe_malloc(allocated * sizeof(char [64]));
	intsFile = safe_open(intsFileName, "r");

    while((fscanf(intsFile, "%s %s %s\n",
		&(ints->protein1[nLine][0]), &(dummy[0]), &(ints->protein2[nLine][0]))) == 3) {
		/*fprintf(stderr, "%s:%d: %d\t%s\t%s\n",
			__FILE__, __LINE__,
			nLine, &(ints->protein1[nLine][0]), &(ints->protein2[nLine][0]));*/
		++ nLine;

		if (nLine == allocated) {
			allocated += 64;
			ints->protein1 = safe_realloc(ints->protein1, allocated * sizeof(char [64]));
			ints->protein2 = safe_realloc(ints->protein2, allocated * sizeof(char [64]));
		}
    }   

	ints->nInteractionListed = nLine;
	assert(ints->nInteractionListed > 0);

	fclose(intsFile);

	if (! arg->silent)
		fprintf(stdout, "\t\tnumber of interactions = %d\n", nLine);
}

/*____________________________________________________________________________*/
/** assign interaction matrix */
void assign_interaction_matrix(int n, Ints *ints, Arg *arg)
{
	unsigned int i, j;
	int p0, p1;

	/* match node numbers and fill interaction matrix */
	/* for all interactions */
	if (! arg->silent) fprintf(stdout, "Computing interaction matrix %d\n", n);

	for (i = 0, ints->nInteraction = 0; i < ints->nInteractionListed; ++ i) {
		/* exclude self-interactions */
		if (! arg->selfInteraction)
			if (strcmp(&(ints->protein1[i][0]), &(ints->protein2[i][0])) == 0)
				continue;

		p0 = -1;
		p1 = -1;
		/* look through node list */
		for (j = 0; j < ints->N; ++ j) {
			/* match node ID numbers of both interacting nodes */
			if (strcmp(&(ints->protein1[i][0]), &(ints->node[j][0])) == 0)
				p0 = j;
			if (strcmp(&(ints->protein2[i][0]), &(ints->node[j][0])) == 0)
				p1 = j;
			if (p0 >= 0 && p1 >= 0) {
				/* assign interaction to matrix */
					/* count interaction as binary [0,1] */
					ints->c[p0][p1] = 1;
					ints->nInteraction += ints->c[p0][p1]; /* count interactions */
					ints->c[p1][p0] = 1;
					ints->nInteraction += ints->c[p1][p0]; /* count interactions */
				break;
			}
		}
		if (p0 < 0 || p1 < 0)
			fprintf(stderr, "Warning: at least one interaction partner of pair %s %s not in node list\n",
				&(ints->protein1[i][0]), &(ints->protein2[i][0]));
	}
	if (! arg->silent)
		fprintf(stdout, "\t\tnumber of interactions = %d\n", ints->nInteraction);
}

/*____________________________________________________________________________*/
/** input is interaction matrix */
int read_interaction_matrix(char *matFileName, Ints *ints, Arg *arg, int nInts)
{
	unsigned int i, j;
	FILE *matInFile = 0;
	int scan = 0;

	/* read input matrix */
	matInFile = safe_open(matFileName, "r");

	/* first line of input matrix is dimension (of square matrix) */
	scan = fscanf(matInFile, "%d", &(ints->N)); /* get number of nodes */
	assert(ints->N > 1); /* safety line */
	/* read only number of interactions if 'nInts' flag is '1' */
	if (nInts) {
		fclose(matInFile);
		return(0);
	}

	if (! arg->silent) fprintf(stdout, "\tReading interaction matrix %s of size N = %d\n",
							matFileName, ints->N);

	/* allocate node name array */
	ints->node = safe_malloc(ints->N * sizeof(int));

	/* read matrix values */
	for (i = 0, ints->nInteraction = 0; i < ints->N; ++ i) {
		/*strcpy(&(ints->node[i][0]), (char)(i + 48));*/ /* set node name array */
		for (j = 0; j < ints->N; ++ j) {
			assert(fscanf(matInFile, "%d", &(ints->c[i][j])) == 1); 
			ints->nInteraction += ints->c[i][j]; /* count interactions */
		}
	}

	fclose(matInFile);

	/* assert symmetry of input matrix */
	if (symmetry_mat2D_int(ints->c, ints->N, ints->N) != 0)
		Error("Input matrix not symmetric!");

	/* assert binary of input matrix */
	if (binary_mat2D_int(ints->c, ints->N, ints->N) != 0)
		Error("Input matrix not binary!");

	if (! arg->silent)
		fprintf(stdout, "\t\tnumber of interactions = %d\n", ints->nInteraction);

	return(0);
}


