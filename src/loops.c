/*==============================================================================
loops.c : network loops 
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "loops.h"
#include "math.h"

int nLoops = 0;

#define MIN_DEPTH 	3	
#define MAX_DEPTH	7
#define NON_SET -1
#define BLOCKED 1

/*____________________________________________________________________________*/
/* implementation of Tiernan algorithm for elementary circuit counting */
void loops_Tiernan(Ints *ints, Arg *arg)
{
	/* Tiernan vars */
	int* P;
	int** H;
	int k,N;
	int i ,ii, j;
	float tp,old_tp; /* time progress */
	int *CC = &(ints->loop_counter[0]);  
	int depth = arg->loopDepth;

	/*____________________________________________________________________________*/
	/* EC1 initialize */
	printf("Progress(%%) :             ");
	N = ints->N;
	P = malloc(N * sizeof(int));
	H = malloc(N * sizeof(int *));
	for (i = 0; i < N; ++ i) {
		P[i] = 0;
		CC[i] = 0;
		H[i] =  malloc(N*sizeof(int));
		for (ii = 0; ii < N; ++ ii)
			H[i][ii] = NON_SET;	
	}
	k = 0;
	P[0] = 0;

	/*____________________________________________________________________________*/
	/* EC2 path extension */
EC2:	
	/* search for next node */
	for (j = P[0] + 1; j < N ; ++ j) {
	
		/* SME::Limit Path Length */
		if (k == depth)
			break;
		
		/* test if j is connected to the last node P[k]	*/
		if (ints->c[P[k]][j] == 0)
			continue;
 		
		/* test that j is not in P(ath) */
		for (i = 1; i <= k; ++ i) 
			if (P[i] == j)
				break;
		if  (P[i] == j)
				continue;
		
		/* test that j is not in block list */
		if (H[P[k]][j] == BLOCKED)
			continue;

		/* 3 if node found extend path */
		k = k + 1;
		P[k] = j;
		goto EC2;
	}
	/* path cannot be extended */

	/*____________________________________________________________________________*/
	/* EC3:	Circuit Confirmation */
	if (ints->c[P[k]][P[0]] == 1) {
		/* increase Circuit Counter for depth k */
		CC[k] += 1;

		/* time progress indicator */
		tp = (P[0] * 1.0 / N) + (P[1] * 1.0 / (N * N)) + (P[2] * 1.0 / (N * N * N) ) ;
		printf("\b\b\b\b\b\b\b\b\b\b\b\b%12f",sqrt(tp) * 100);
	}
	
	/* if all initial node P[0] circuits are counted go to next inital node */
	if (k == 0)
		goto EC5;

	/* clean block list H for P[k] */
	for (i = 0; i < N; ++ i)
		H[P[k]][i] = NON_SET;

	/* Block P[k] */
	H[P[k-1]][P[k]] = BLOCKED;
	P[k] = 0;

	/* backtrack and continue counting */
	k = k - 1;
	goto EC2;

	/*____________________________________________________________________________*/
	/* EC5: Advance Initial Vertex */
EC5:
	
	//If all initial nodes are processed terminate execution
	if (P[0] == N - 1)
		goto EC6;
	
	//cleanup and count from next initial node
	P[0] = P[0] + 1;
	k = 0;
	for (i = 0; i < N; ++ i) {
		for (ii = 0; ii < N; ++ ii)
			H[i][ii] = NON_SET;
	}
	goto EC2;

	/*____________________________________________________________________________*/
	/* EC6: Terminate */
EC6:
		printf("\b\b\b\b\b\b\b\b\b\b\b\b%12f",100.0); /* finished 100% progress */
		/* print_path(CC,depth); */
		free(P);
		for (i = 0; i < N; ++ i)
			free(H[i]);
		free(H);
}

/*____________________________________________________________________________*/
/* count graph loops */
void loops(Ints *ints, Arg *arg)
{
	int i;
	printf("\tLoops : ");
	
	ints->loop_counter = malloc(ints->N * sizeof(int));
	loops_Tiernan(ints, arg);

	printf("\n");
}

/*____________________________________________________________________________*/
/* display a path in the graph */
void print_path(int* P, int k) {
	int i;

	for (i = 0; i <= k; ++ i)
			printf("%d ",P[i]);
		printf("\n");
}

/*____________________________________________________________________________*/
/* print loops to output file */
void print_loops(Ints *ints, Arg *arg)
{
	FILE *outFile;
	unsigned int i, j;
	 
	outFile = safe_open("loops.dat", "w");

	fprintf(outFile,"loop search depth: %d\n", arg->loopDepth);
	fflush(stdout);
	/* Tiernan is a directional graph algorithm, */
	/*    loop number be must divided by 2 for undirected graphs */ 
	for (i = 3; i <= arg->loopDepth; ++ i) {
		fprintf(outFile, "%8d\t%d\n", i, (int)floor(.5 * ints->loop_counter[i])); 
	}
	fclose(outFile);
}

