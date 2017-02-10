/*==============================================================================
ints.h : interaction data structure
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef INTS_H
#define INTS_H

/*____________________________________________________________________________*/
/* interaction structure */
typedef struct  
{
	/* input data */
	int N; /* number of nodes */
	char (*node)[64]; /* array of node ID numbers */
	int nInteractionListed; /* number of listed node interactions */
	int nInteraction;
	char (*protein1)[64]; /* array of ineraction proteins */
	char (*protein2)[64]; /* first (1) and second (2) interactor */
	int **c; /* interaction matrix */
	int *k_list; /* list of 'degree assigned to node i' */

	/* degree stats */
	int k_max; /* maximum degree in interaction matrix */
	double k_av; /* degree average */
	double k_var; /* degree variance */
	double k_sq_av; /* degree square average */
	double logk_av; /* <log k> */
	int *k; /* degree frequency distribution k(c) : 
				frequency of of all degrees k_i over all vertices i */

	/* degree-degree statistics */
	double *w; /*  sum over all node pairs with degrees (k,k1), independent of connectivity */
	double **W; /* probability for randomly drawn nodes with degrees (k,k1) */
	double *p; /* degree probability p(k) : p_i(c) = 1/N k_i(c) */
	double **Pi; /* degree correlation matrix */

	/* entropy */
	double S_N; /* N-dependent 'zero' entropy S_N */
	double C_p; /* complexity C_P of degree distribution */
	double C_Pi; /* complexity C_Pi of degree correlation */
	double S; /* Shannon entropy of network : S[p,Pi] = S_N - C_p - C_Pi */
	double C; /* complexity of network: C[p,Pi] = C_p + C_Pi */
	double S_N_pl; /* the same properties as above, but exressed 'per link' instead of 'per node' */
	double C_p_pl;
	double C_Pi_pl;
	double S_pl;
	double C_pl;
	double R; /* assortativity */
    int loop_search_depth;
    int* loop_counter; /* loop_count[k] = number of loops (elementary circuits) of length k where k <= loop_search_depth. */ 

	/* for distance computation, i.e. expanded to the larger matrix */
	double *p_dist; /* p for distance computation */
	double *p_smooth; /* p_dist in smoothened form */
	double **Pi_dist; /* Pi for distance computation */
	double **Pi_smooth; /* Pi_dist in smoothened form */
	double *w_dist; /* degree distribution */
	double *w_smooth; /* w in smoothened form */
	double **W_dist; /* degree correlation matrix */
	double **W_smooth; /* W in smoothened form */
} Ints;

/* distance structure */
typedef struct  
{
	/* Kullback-Leibler distance terms */
	double D_A_p; /* mutual information of function p for network A */
	double D_A_Pi; /* mutual information of function Pi for network A */
	double D_B_p; /* mutual information of function p for network B */
	double D_B_Pi; /* mutual information of function Pi for network B */
	double D_AB; /* Kullback-Leibler distance between two networks A and B */

	/* network-pair specific constants */
	int nSmoothLevel; /* number of smoothing levels */
	int k_max_pair; /* maximum k of two networks */
	int k_min_pair;
	int N_max_pair; /* maximum N of two networks */
	int N_min_pair; /* minimum N of two networks */
	double fillval; /* value to fill in matrix zeros; can be overwritten by 'arg->fillval' */
} Dist;

#endif

