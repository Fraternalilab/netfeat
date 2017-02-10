/*==============================================================================
putprop.c : print network properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "putprop.h"

/*____________________________________________________________________________*/
/** print degree per node */
void print_k_list(int n, Ints *ints)
{
	char outFileName[64];
	sprintf(outFileName, "k_list.%c.dat", (n+48));
	print_array_int(outFileName, ints->k_list, ints->N);
}

/*____________________________________________________________________________*/
/** print degree distribution */
void print_k(int n, Ints *ints)
{
	char outFileName[64];
	sprintf(outFileName, "k.%c.dat", (n+48));
	print_array_int(outFileName, ints->k, (ints->k_max + 1));
}

/*____________________________________________________________________________*/
/** print degree probability distribution */
void print_p(int n, Ints *ints)
{
	char outFileName[64];
	sprintf(outFileName, "p.%c.dat", (n+48));
	print_array_double(outFileName, ints->p, (ints->k_max + 1));
}

/*____________________________________________________________________________*/
/** print connection probability matrix */
void print_P_conn(int n, Ints *ints)
{
	char outFileName[64];
	sprintf(outFileName, "P_conn.%c.dat", (n+48));
	print_mat2D_double_lowlim(outFileName, ints->W, 1, (ints->k_max + 1), 1, (ints->k_max + 1));
}

/*____________________________________________________________________________*/
/** print degree correlation matrix */
void print_Pi_corr(int n, int level, Ints *ints)
{
	char outFileName[64];
	sprintf(outFileName, "Pi_corr.%c.%c.mat", (n+48), (level+48));
	if (level == 0)
		print_mat2D_doublee_lowlim(outFileName, ints->Pi, 1, (ints->k_max + 1), 1, (ints->k_max + 1));
	/*
	else
		print_mat2D_double_lowlim(outFileName, ints->Pi_smooth[level-1], 1, (ints->k_max + 1), 1, (ints->k_max + 1));
	*/
}

/*____________________________________________________________________________*/
/** print other useful information */
void print_etc(int n, Ints *ints)
{
	FILE *outFile;
	char outFileName[64];
	sprintf(outFileName, "etc.%c.dat", (n+48));

	outFile = safe_open(outFileName, "w");

	fprintf(outFile, \
		"nodes %d\ninteractions %d\nk_max %d\n<k> %lf\n<kk> %lf\nassortativity %lf\n",
		ints->N, ints->nInteraction, ints->k_max, ints->k_av, ints->k_var, ints->R);

	fclose(outFile);
}

/*____________________________________________________________________________*/
/** print entropy values */
void print_entropy(int n, Ints *ints)
{
	FILE *outFile;
	char outFileName[64];
	sprintf(outFileName, "entropytab.%c.dat", (n+48));

	outFile = safe_open(outFileName, "w");

	fprintf(outFile, \
		"k_av %lf\n\nN_entropy %lf\np_entropy %lf\nPi_entropy %lf\nS_entropy %lf\nC_complexity %lf\n"
		"\nN_entropy_pl %lf\np_entropy_pl %lf\nPi_entropy_pl %lf\nS_entropy_pl %lf\nC_complexity_pl %lf\n",
		ints->k_av, ints->S_N, ints->C_p, ints->C_Pi, ints->S, ints->C,
		ints->S_N_pl, ints->C_p_pl, ints->C_Pi_pl, ints->S_pl, ints->C_pl);

	fclose(outFile);
}

