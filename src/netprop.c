/*==============================================================================
netprop.c : macroscopic network properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netprop.h"

/*____________________________________________________________________________*/
/** compute log factorial function:
	log k! */
__inline__ static double log_factorial(int k)
{
	unsigned int i;
	double sum = 0.;
	double log();

	if (k > 1)
		for(i = 2; i <= k; ++ i)
			sum += log(i);

	return sum;
}

/*____________________________________________________________________________*/
/** compute log Poisson function:
	log Pi(k) = k*log(<k>) - <k> - SUM_{l=1}^{k}(log(l)) */
__inline__ static double log_poisson(int k, double k_av)
{
	unsigned int i;
	double sum = 0.;
	double log();

	return (k * log(k_av) - k_av - log_factorial(k));
}

/*____________________________________________________________________________*/
/** print degree (int) matrix in ASCII format */
void print_degree_matrix_ascii_int(char *outFileName, int **mat, int k_max)
{
	FILE *outFile;
	unsigned int i, j;

	fprintf(stdout, "Printing degree matrix '%s'\n", outFileName);
	outFile = safe_open(outFileName, "w");

	/* print matrix dimension in first line */
	fprintf(outFile, "%d\n", (k_max + 1));

	/* print matrix values*/
	/* columns */
	for (i = 0; i <= k_max; ++ i) {
		if (i > 0)
			fprintf(outFile, "\n");

		/* progress info (for large files) */
		/* if (i % 1000 == 0)
			fprintf(stderr, "%d / %d\n", i, k_max); */

		/* rows */
		for (j = 0; j <= k_max; ++ j)
			fprintf(outFile, "%d ", mat[i][j]);
	}
	fprintf(outFile, "\n");
	
	fclose(outFile);
}

/*____________________________________________________________________________*/
/** print degree (double) matrix in ASCII format */
void print_degree_matrix_ascii_double(char *outFileName, double **mat, int k_max)
{
	FILE *outFile;
	unsigned int i, j;

	fprintf(stdout, "Printing degree matrix '%s'\n", outFileName);
	outFile = safe_open(outFileName, "w");

	/* print matrix dimension in first line */
	fprintf(outFile, "%d\n", (k_max + 1));

	/* print matrix values*/
	/* columns */
	for (i = 0; i <= k_max; ++ i) {
		if (i > 0)
			fprintf(outFile, "\n");

		/* progress info (for large files) */
		/* if (i % 1000 == 0)
			fprintf(stderr, "%d / %d\n", i, k_max); */

		/* rows */
		for (j = 0; j <= k_max; ++ j)
			fprintf(outFile, "%6.5f ", mat[i][j]);
	}
	fprintf(outFile, "\n");
	
	fclose(outFile);
}

/*____________________________________________________________________________*/
/** equality test with rounding error tolerance for double values */
__inline__ static int approximately_equal(double a, double b)
{
    return ((a == b) || (fabsf(a - b) <= .00001F * fabsf(a) + .00001F));
}

/*____________________________________________________________________________*/
/** degree_distribution (manuscript section 2.1, equation 1) :
	k_av, k_var, k_max,
	compute degree frequencies k_i, the overall frequency distribution k[k_i] 
	and probability distribution p[k_i] of degrees;
	initialise for range [0, k_max] */
void degree_statistics(Ints *ints, Arg *arg)
{
	unsigned int i, j, k, k_max;
	int k_i; /* degree k of node i: number of links to this node */
	int checksum_i = 0; /* for consistency check */
	double checksum_f = 0.;

	/** initialise variables */
	ints->k_av = 0.;
	ints->k_var = 0.;
	ints->logk_av = 0.;
	ints->k_max = 0;

	/*____________________________________________________________________________*/
	/** determine k_max */
	for (i = 0; i < ints->N; ++ i) {
		k = 0;
		for (j = 0; j < ints->N; ++ j) {
			k += ints->c[i][j]; /* degree 'k_i': sum up all interactions of node 'i' */
		}
		k_max = k_max > k ? k_max : k; /* k_max */
	}
	ints->k_max = k_max;

	/*____________________________________________________________________________*/
	/** initialise arrays */
	/* degree distribution */
	ints->k = safe_malloc((ints->k_max + 1) * sizeof(int));
	init_array_int(ints->k, (ints->k_max + 1), 0);
	/* empirical degree distribution */
	ints->p = safe_malloc((ints->k_max + 1) * sizeof(double));
	init_array_double(ints->p, (ints->k_max + 1), 0.);
	/* initialise degree list array */
	ints->k_list = safe_malloc(ints->N * sizeof(int));
	init_array_int(ints->k_list, ints->N, 0);

	/*____________________________________________________________________________*/
	/** k_av, k_var, k_max, degree frequency distribution */
	for (i = 0; i < ints->N; ++ i) {
		for (j = 0, k_i = 0; j < ints->N; ++ j) {
			k_i += ints->c[i][j]; /* degree 'k_i': sum up all interactions of node 'i' */
		}
		++ ints->k[k_i]; /* k_i(c) : the frequency of degree 'k_i' */
		ints->k_list[i] = k_i; /* degree list: assign degree 'k_i' to node 'i' */
		ints->k_av += ((k_i - ints->k_av) / (i + 1)); /* average degrees */
		ints->k_var += (((k_i * k_i) - ints->k_var) / (i + 1)); /* degree variation */
	}

#ifdef DEBUG
	for (i = 0; i < ints->N; ++ i)
		fprintf(stderr, "%s:%d: %d %d\n", __FILE__, __LINE__, i, ints->k_list[i]);
#endif

	/*____________________________________________________________________________*/
	/** degree probability distribution : p(i) = 1/N k(i) */
	for (i = 0; i <= ints->k_max; ++ i) {
		if (ints->k[i] > 0)
			ints->p[i] = (double)ints->k[i] / (double)ints->N;

		/** variables for consistency checks (see below) */
		checksum_i += (ints->k[i] * i);
		checksum_f += ints->p[i];
#ifdef DEBUG
		if (ints->p[i] > 1e-6) {
		fprintf(stderr, "%s:%d: %d %d %f\n", 
			__FILE__, __LINE__, i, ints->k[i], ints->p[i]);
		}
#endif
	}

	/*____________________________________________________________________________*/
	/** consistency checks */
	/* sum over all [degree_probability * degree] equals number of interactions */	
#ifdef DEBUG
	fprintf(stderr, "%s:%d: checksum_i %d\tints->nInteraction %d\n",
				__FILE__, __LINE__, checksum_i, ints->nInteraction);
#endif
	assert(checksum_i == ints->nInteraction);
	/* the probability distribution is normalised */
	assert(approximately_equal(checksum_f, 1.));

	/* print result */
	if (! arg->silent) fprintf(stdout, "\t\taverage degree <k> = %lf\n\t\tdegree variance <kk> = %lf\n\t\tmaximal degree k_max = %d\n",
								ints->k_av, ints->k_var, ints->k_max);
	if (ints->k_max <= 0) Error("\t\tmaximal degree <= 0: something wrong here?\n");
}

/*____________________________________________________________________________*/
/** wiring (manuscript section 2.1, equation 3) :
	W: probability for two randomly drawn nodes with degrees (k,k') 
	to be connected; symmetric matrix */
void wiring(Ints *ints, Arg *arg)
{
	unsigned int i, j;
	/* relative degree probability */
	ints->w = safe_malloc((ints->N + 1) * sizeof(double));
	init_array_double(ints->w, (ints->N + 1), 0.);
	/* probability for randomly drawn nodes with degrees (k,k1) */
	ints->W = alloc_mat2D_double(ints->W, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->W, (ints->k_max + 1), (ints->k_max + 1), 0.);

	/*____________________________________________________________________________*/
	/** w */ 
    for (i = 0; i < ints->N; ++ i) {
		ints->w[i] = ints->p[i] * (double)i / ints->k_av;
    }

    /*____________________________________________________________________________*/
    /** W : loops run over nodes, but W is function of degrees */
    for (i = 0; i <= ints->N; ++ i) {
        for (j = 0; j <= ints->N; ++ j) {
			if (ints->c[i][j] == 1) {
				ints->W[ints->k_list[i]][ints->k_list[j]] += (1. / (ints->k_av * (double)ints->N));
			}
        }
    }

	if (symmetry_mat2D_double(ints->W, (ints->k_max + 1), (ints->k_max + 1)) != 0)
		Error("Matrix W not symmetric");
}

/*____________________________________________________________________________*/
/** degree-degree correlation function Pi[k_i][k_j] 
	(manuscript section 2.1, equation 6) */
void degree_degree_distribution_eq6(Ints *ints, Arg *arg)
{
	unsigned int i, j;

	ints->Pi = alloc_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1), 0.);

	/* loop starts at '1': degree '0' has no correlation,
		because a node with degree '0' has no link */
	for (i = 1; i <= ints->k_max; ++ i) {
		for (j = 1; j <= ints->k_max; ++ j) {
			if (ints->W[i][j] > 0.) {
				ints->Pi[i][j] = ints->W[i][j] * ints->k_av * (ints->N - 1) / (i * j);
			}
		}
	}
	if (symmetry_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1)) != 0)
		Error("Matrix Pi not symmetric");
}

/*____________________________________________________________________________*/
/** manuscript section 3.4, equation 35: Shannon entropy per node
	and section 4.1, equation 39: complexity */
void entropy_pn(int n, Ints *ints, Arg *arg)
{
	FILE *entropyOutFile = 0;
	char entropyOutFileName[64];
	unsigned int i, j;
	double log();

	/* initialise values */
	ints->S_N = 0.; /* N-dependent 'zero' entropy S_N */
	ints->C_p = 0.; /* complexity C_P of degree distribution */
	ints->C_Pi = 0.; /* complexity C_Pi of degree correlation */
	ints->S = 0.; /* total Shannon entropy of network : S = S_N - C_p - C_Pi */

	/*____________________________________________________________________________*/
	/* individual computation of the three terms of equation 35 */
	/** term 1 : N-dependent 'zero' entropy S_0 */
	ints->S_N = .5 * ints->k_av  * (log(ints->N / ints->k_av) + 1);

	/*____________________________________________________________________________*/
	/** term 2 : complexity C_p of degree distribution */ 
	for (i = 0; i <= ints->k_max; ++ i)
		if (ints->p[i] > 0.)
			ints->C_p += (ints->p[i] * (log(ints->p[i]) - log_poisson(i, ints->k_av)));

	/*____________________________________________________________________________*/
	/** term 3 : complexity C_Pi of degree correlation (wiring complexity) */
	for (i = 1; i <= ints->k_max; ++ i)
		for (j = 1; j <= ints->k_max; ++ j)
			if ((ints->p[i] > 0.) && (ints->p[j] > 0.) && (ints->Pi[i][j] > 0.))
				ints->C_Pi += (ints->p[i]  * ints->p[j] * i * j * \
							ints->Pi[i][j] * log(ints->Pi[i][j]));
	ints->C_Pi /= (2 * ints->k_av);

	/*____________________________________________________________________________*/
	/** total entropy */
	ints->S = ints->S_N - ints->C_p - ints->C_Pi;

	/* complexity */
	ints->C = ints->C_p + ints->C_Pi;

	/*fprintf(stderr, "%s:%d: term1 = %e, term2 = %e, term3 = %e, total = %e\n",
		__FILE__, __LINE__, ints->S_N, ints->C_p, ints->C_Pi, ints->S);*/

	/*____________________________________________________________________________*/
	/** output to file */
	sprintf(entropyOutFileName, "entropy.%c.dat", (n+48));
	entropyOutFile = safe_open(entropyOutFileName, "w");
	fprintf(entropyOutFile, "S_0 %f\nC_p %f\nC_Pi %f\nS %f\n",
		ints->S_N, ints->C_p, ints->C_Pi, ints->S);
	fclose(entropyOutFile);
	
	/*____________________________________________________________________________*/
	/** output to stdout */
	if (! arg->silent)
		fprintf(stdout, "\n\t\t(per node)\n"
						"\t\tsize entropy S_0 = %f\n"
						"\t\tdegree complexity C_p = %f\n"
						"\t\tPi complexity C_Pi = %f\n"
						"\t\ttotal entropy S = %f\n",
			ints->S_N, ints->C_p, ints->C_Pi, ints->S);
}

/*____________________________________________________________________________*/
/** Shannon entropy per link (not in manuscript):
	the same as the Shannon entropy above,
	but entropy expressed 'per link (pl)' instead of 'per node',
	which removes any artifact coming from fluctuations in link numbers;
	the change is a division by 1/(.5*k_av) */
void entropy_pl(int n, Ints *ints, Arg *arg)
{
	FILE *entropyplOutFile = 0;
	char entropyplOutFileName[64];
	unsigned int i, j;
	double log();
	double scalingFactor = .5 * ints->k_av;

	/* initialise values */
	ints->S_N_pl = 0.;
	ints->C_p_pl = 0.;
	ints->C_Pi_pl = 0.;
	ints->S_pl = 0.;
	ints->C_pl = 0.;

	/*____________________________________________________________________________*/
	ints->S_N_pl = ints->S_N / scalingFactor;
	ints->C_p_pl = ints->C_p / scalingFactor;
	ints->C_Pi_pl = ints->C_Pi / scalingFactor;

	/*____________________________________________________________________________*/
	ints->S_pl = ints->S / scalingFactor;
	ints->C_pl = ints->C / scalingFactor;

	/*____________________________________________________________________________*/
	/** output to file */
	sprintf(entropyplOutFileName, "entropypl.%c.dat", (n+48));
	entropyplOutFile = safe_open(entropyplOutFileName, "w");
	fprintf(entropyplOutFile, "S_0_pl %f\nC_p_pl %f\nC_Pi_pl %f\nS_pl %f\n",
		ints->S_N_pl, ints->C_p_pl, ints->C_Pi_pl, ints->S_pl);
	fclose(entropyplOutFile);
	
	/*____________________________________________________________________________*/
	/** output to stdout */
	if (! arg->silent)
		fprintf(stdout, "\n\t\t(per link)\n"
						"\t\tsize entropy S_0_pl = %f\n"
						"\t\tdegree complexity C_p_pl = %f\n"
						"\t\tPi complexity C_Pi_pl = %f\n"
						"\t\ttotal entropy S_pl = %f\n",
			ints->S_N_pl, ints->C_p_pl, ints->C_Pi_pl, ints->S_pl);
}

/*____________________________________________________________________________*/
/** assortativity : following Ton's implementation
	Holme and Zhao, "Exploring the assortativity-clustering space of a 
	network's degree sequence", Physical Review E 75, 046111 (2007) */
void assortativity(Ints *ints, Arg *arg)
{
	unsigned int i, j;
	double correlation = 0.;
	double average = 0.;
	double variance = 0.;
	double norm = 0.;
	int k_i = 0.;
	int k_j = 0.;
	int k_i_sq = 0.;
	int k_j_sq = 0.;
	double pow();
	ints->R = 0.;

    for (i = 0; i < ints->N; ++ i) {
		for (j = 0; j < ints->N; ++ j) {
			if ((ints->k_list[i] > 0) && (ints->k_list[j] > 0)) {
				k_i = (double)ints->k_list[i];
				k_j = (double)ints->k_list[j];
				k_i_sq = pow(k_i, 2);
				k_j_sq = pow(k_j, 2);

				if (ints->Pi[k_i][k_j] > 0.) {
					correlation += (k_i_sq *k_j_sq * ints->Pi[k_i][k_j]);
					average = average + 0.5 * (k_i + k_j) * k_i * k_j * ints->Pi[k_i][k_j];
					variance = variance + 0.5 * (k_i_sq + k_j_sq) * k_i * k_j * ints->Pi[k_i][k_j];
					norm = norm + (k_i * k_j * ints->Pi[k_i][k_j]);
				}
			}
		}
	}

	correlation /= norm;
	average /= norm;
	variance /= norm;

	ints->R = (correlation - (average * average)) / (variance - (average * average));

	if (! arg->silent)
		fprintf(stdout, "\n\t\tassortativity R = %lf\n", ints->R);
}

/*____________________________________________________________________________*/
/** compute (single) network properties */
void network_properties(int n, Ints *ints, Arg *arg)
{

}

