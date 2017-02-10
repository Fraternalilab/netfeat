/*==============================================================================
netprop.c : macroscoPic network properties
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
/** assign maximal degree from matrix to 'k_max',
	which is needed for all loops over the range of degrees [0, k_max] */
void max_degree(Ints *ints, Arg *arg)
{
	unsigned int i, j, k;
	ints->k_max = 0;

	for (i = 0; i < ints->N; ++ i) {
		for (j = 0, k = 0; j < ints->N; ++ j) {
			k += ints->c[i][j];
		}
		ints->k_max = ints->k_max > k ? ints->k_max : k;
	}
	if (! arg->silent) fprintf(stdout, "\tmaximal degree k_max = %d\n", ints->k_max);
	if (ints->k_max <= 0) Error("\tmaximal degree <= 0: something wrong here?\n");
}

/*____________________________________________________________________________*/
/** average degree from p(k) array */
void average_degree(Ints *ints, Arg *arg)
{
	unsigned int i;

	for (i = 0, ints->k_av = 0.; i <= ints->k_max; ++ i)
		ints->k_av += (ints->p[i] * i);

	if (! arg->silent)
		fprintf(stdout, "\taverage degree <k> = %f\n", ints->k_av);
	if (ints->k_av <= 0) Error("\taverage degree <= 0: something wrong here?\n");
}

/*____________________________________________________________________________*/
/** equality test with rounding error tolerance for double values */
__inline__ static int approximately_equal(double a, double b)
{
    return ((a == b) || (fabsf(a - b) <= .00001F * fabsf(a) + .00001F));
}

/*____________________________________________________________________________*/
/** manuscript section 2.1, equation 1 :
	compute degree frequencies k_i, the overall frequency distribution k[k_i] 
	and probability distribution p[k_i] of degrees;
	initialise for range [0, k_max] */
void equation_1(Ints *ints, Arg *arg)
{
	unsigned int i, j;
	int k_i; /* degree k of node i: number of links to this node */
	int checksum_i = 0; /* for consistency check */
	double checksum_f = 0.;

	/** initialise variables */
	ints->k_av = 0.;
	ints->k_var = 0.;
	ints->logk_av = 0.;
	/** initialise arrays */
	/* degree distribution */
	ints->k = safe_malloc((ints->k_max + 1) * sizeof(int));
	init_array_int(ints->k, (ints->k_max + 1), 0);
	/* emPirical degree distribution */
	ints->p = safe_malloc((ints->k_max + 1) * sizeof(double));
	init_array_double(ints->p, (ints->k_max + 1), 0.);
	/* initialise degree list array */
	ints->k_list = safe_malloc(ints->N * sizeof(int));
	init_array_int(ints->k_list, ints->N, 0);

	/*____________________________________________________________________________*/
	/** compute degree distribution */
	for (i = 0; i < ints->N; ++ i) {
		for (j = 0, k_i = 0; j < ints->N; ++ j) {
			k_i += ints->c[i][j]; /* degree 'k_i': sum up all interactions of node 'i' */
		}
		++ ints->k[k_i]; /* k_i(c) : the frequency of degree 'k_i' */
		ints->k_list[i] = k_i; /* degree list: assign degree 'k_i' to node 'i' */
		ints->k_av += ((k_i - ints->k_av) / (i + 1)); /* average degrees */
		ints->k_var += (((k_i * k_i) - ints->k_var) / (i + 1)); /* degree variation */
	}

	/*for (i = 0; i < ints->N; ++ i)
		fprintf(stderr, "%s:%d: %d %d\n", __FILE__, __LINE__, i, ints->k_list[i]);*/

	/*____________________________________________________________________________*/
	/** compute emPirical degree probability distribution : p(i) = 1/N k(i) */
	for (i = 0; i <= ints->k_max; ++ i) {
		if (ints->k[i] > 0)
			ints->p[i] = (double)ints->k[i] / (double)ints->N;

		/** variables for consistency checks (see below) */
		checksum_i += (ints->k[i] * i);
		checksum_f += ints->p[i];
		/*fprintf(stderr, "%s:%d: %d %d %f\n", 
			__FILE__, __LINE__, i, ints->k[i], ints->p[i]);*/
	}

	/*____________________________________________________________________________*/
	/** consistency checks */
	/* sum over all [degree_probability * degree] equals number of interactions */	
	assert(checksum_i == ints->nInteraction);
	/* the probability distribution is normalised */
	assert(approximately_equal(checksum_f, 1.));

	/* print result */
	if (! arg->silent) fprintf(stdout, "\taverage degree <k> = %lf\n\tdegree variance <kk> = %lf\n",
		ints->k_av, ints->k_var);
}

/*____________________________________________________________________________*/
/** manuscript section 2.1, equation 3 :
	P_conn_k_k1: probability for two randomly drawn nodes with degrees (k,k') 
	to be connected; symmetric matrix */
void equation_3(Ints *ints, Arg *arg)
{
	unsigned int i, j;
	/** initialise arrays */
	/* sum over all connected node pairs with degrees (k,k1) */
	ints->conn_k_k1 = alloc_mat2D_int(ints->conn_k_k1, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_int(ints->conn_k_k1, (ints->k_max + 1), (ints->k_max + 1), 0);
	/* sum over all node pairs with degrees (k,k'), independent of their connectivity */
	ints->all_k_k1 = alloc_mat2D_int(ints->all_k_k1, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_int(ints->all_k_k1, (ints->k_max + 1), (ints->k_max + 1), 0);
	/* probability for randomly drawn nodes with degrees (k,k1) */
	ints->P_conn_k_k1 = alloc_mat2D_double(ints->P_conn_k_k1, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->P_conn_k_k1, (ints->k_max + 1), (ints->k_max + 1), 0.);

	/*____________________________________________________________________________*/
	/**  probability for two randomly drawn nodes with degrees (k,k') to be 
		connected, un-normalised */
    for (i = 0; i < ints->N; ++ i) {
        for (j = 0; j < i; ++ j) { /* i != j */
            /* sum over connected node paris with degrees (k,k') */
			if (ints->c[i][j] == 1) {
				++ ints->conn_k_k1[ints->k_list[i]][ints->k_list[j]];
				++ ints->conn_k_k1[ints->k_list[j]][ints->k_list[i]];
			}
            /* sum over all node pairs with degrees (k,k') */
            ++ ints->all_k_k1[ints->k_list[i]][ints->k_list[j]];
            ++ ints->all_k_k1[ints->k_list[j]][ints->k_list[i]];

            /*fprintf(stderr, "%s:%d: %c(%d)\t%c(%d)\t%d\t%f\t%d\n",
				 __FILE__, __LINE__,
				i+65, ints->k_list[i], j+65, ints->k_list[j], ints->intMat[i][j],
				ints->P_conn_k_k1[ints->k_list[i]][ints->k_list[j]],
				ints->all_k_k1[ints->k_list[i]][ints->k_list[j]]);*/
        }
    }

    /*____________________________________________________________________________*/
    /** ratio of conn_k_k1 and all_k_k1 */
    for (i = 0; i <= ints->k_max; ++ i) {
        for (j = 0; j <= ints->k_max; ++ j) {
			if (ints->all_k_k1[i][j] > 0) {
				ints->P_conn_k_k1[i][j] = \
					(double)ints->conn_k_k1[i][j] / (double)ints->all_k_k1[i][j];
			} else {
				ints->P_conn_k_k1[i][j] = -1.;
            }

			/*fprintf(stderr, "%s:%d: %d\t%d\t%d\t%d\t%f\n",
				__FILE__, __LINE__,
				i, j, ints->conn_k_k1[i][j], ints->all_k_k1[i][j],
				ints->P_conn_k_k1[i][j]);*/
        }
    }
	if (symmetry_mat2D_double(ints->P_conn_k_k1, (ints->k_max + 1), (ints->k_max + 1)) != 0)
		Error("Matrix P_conn_k_k1 not symmetric");
}

/*____________________________________________________________________________*/
/** manuscript section 2.1, equation 6 : 
	Pi: degree correlation function; symmetric matrix */
void equation_6(Ints *ints, Arg *arg)
{
	unsigned int i, j;

	ints->Pi = alloc_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1), 0.);

	/* loop starts at '1': degree '0' has no correlation,
		because a node with degree '0' has no link */
	for (i = 1; i <= ints->k_max; ++ i) {
		for (j = 1; j <= ints->k_max; ++ j) {
			if (ints->P_conn_k_k1[i][j] > 0.) {
				ints->Pi[i][j] = \
					ints->P_conn_k_k1[i][j] * ints->k_av * (ints->N - 1) / (i * j);
			} else {
				if (ints->p[i] == 0 || ints->p[j] == 0)
					ints->Pi[i][j] = 1.;
			}
			/*fprintf(stderr, "%s:%d: %d\t%d\t%f\t%f\n",
				__FILE__, __LINE__, i, j, ints->P_conn_k_k1[i][j], ints->Pi[i][j]);*/
		}
	}
	if (symmetry_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1)) != 0)
		Error("Matrix Pi not symmetric");
}

/*____________________________________________________________________________*/
/** manuscript section 2.2, equation 8 : 
	Pi: degree correlation function; symmetric matrix */
void equation_8(Ints *ints, Arg *arg)
{
	unsigned int i, j;
	double value, factor;

	ints->Pi = alloc_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1), 0.);

	factor = ints->k_av * (1. - 1. / (double)ints->N)/ (double)ints->N;

	/* loop starts at '1': degree '0' has no correlation,
		because a node with degree '0' has no link */
	for (i = 1; i <= ints->k_max; ++ i) {
		for (j = 1; j <= i; ++ j) {
			if (ints->conn_k_k1[i][j] > 0) {
				if (i != j) {
					value = ints->conn_k_k1[i][j] * factor / (double)(i * j * ints->p[i] * ints->p[j]);
				} else {
					value = ints->conn_k_k1[i][j] * factor / (double)(i * j * ints->p[i] * (ints->p[j] - 1. / (double)ints->N));
				}
				ints->Pi[i][j] = ints->Pi[j][i] =	value;
			} else {
				if (ints->p[i] == 0 || ints->p[j] == 0)
					ints->Pi[i][j] = ints->Pi[j][i] = 1.;
			}
			/*fprintf(stderr, "%s:%d: %d\t%d\t%f\t%f\n",
				__FILE__, __LINE__, i, j, ints->P_conn_k_k1[i][j], ints->Pi[i][j]);*/
		}
	}
	if (symmetry_mat2D_double(ints->Pi, (ints->k_max + 1), (ints->k_max + 1)) != 0)
		Error("Matrix Pi not symmetric");
}

/*____________________________________________________________________________*/
/** manuscript section 3.4, equation 35: Shannon entropy 
	and section 4.1, equation 39: complexity */
void equation_35(int n, Ints *ints, Arg *arg)
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
		fprintf(stdout, "\n\t(per node)\n"
						"\tsize entropy S_0 = %f\n"
						"\tdegree complexity C_p = %f\n"
						"\tPi complexity C_Pi = %f\n"
						"\ttotal entropy S = %f\n",
			ints->S_N, ints->C_p, ints->C_Pi, ints->S);
}

/*____________________________________________________________________________*/
/** per-link Shannon entropy (not in manuscript):
	the same as the Shannon entropy above,
	but entropy expressed 'per link (pl)' instead of 'per node',
	which removes any artifact coming from fluctuations in link numbers;
	the change is a division by 1/(.5*k_av) */
void equation_35_pl(int n, Ints *ints, Arg *arg)
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
		fprintf(stdout, "\n\t(per link)\n"
						"\tsize entropy S_0_pl = %f\n"
						"\tdegree complexity C_p_pl = %f\n"
						"\tPi complexity C_Pi_pl = %f\n"
						"\ttotal entropy S_pl = %f\n",
			ints->S_N_pl, ints->C_p_pl, ints->C_Pi_pl, ints->S_pl);
}

/*____________________________________________________________________________*/
/** mutual information of function p */
__inline__ static double mip(Ints *ints0, Ints *ints1, int k_max, Arg *arg)
{
	unsigned int i;
	double D_p = 0.;
	double log();

	for (i = 1; i <= k_max; ++ i) {
		if ((ints0->p[i] > 0.) && (ints1->p[i] > 0.)) {
			D_p += (double)(ints0->p[i] * log(ints0->p[i] / ints1->p[i]));

			/* fprintf(stderr, "%s:%d: %d\t%le\t%le\t%le\n",
			__FILE__, __LINE__,
			i, ints0->p[i], log(ints0->p[i]), ints1->p[i]);*/
		}
	}

	return D_p;
}

/*____________________________________________________________________________*/
/** smoothed mutual information of function p */
__inline__ static double mip_smooth(Ints *ints0, Ints *ints1, int k_max, Arg *arg)
{
	unsigned int i;
	double D_p = 0.;
	double log();

	for (i = 1; i <= k_max; ++ i) {
		if ((ints0->p_smooth[i] > 0.) && (ints1->p_smooth[i] > 0.)) {
			D_p += (double)(ints0->p_smooth[i] * log(ints0->p_smooth[i] / ints1->p_smooth[i]));

			/* fprintf(stderr, "%s:%d: %d\t%le\t%le\t%le\n",
			__FILE__, __LINE__,
			i, ints0->p_smooth[i], log(ints0->p_smooth[i]), ints1->p_smooth[i]); */
		}
	}

	return D_p;
}

/*____________________________________________________________________________*/
/** mutual information of function Pi */
__inline__ static double miPi(Ints *ints0, Ints *ints1, int k_max, Arg *arg)
{
	unsigned int i, j;
	double D_Pi = 0.;
	double log();

	for (i = 1; i <= k_max; ++ i) {
		for (j = 1; j <= k_max; ++ j) {
			if ((ints0->p[i] > 0.) && (ints0->p[j] > 0.) && \
				(ints0->Pi[i][j] > 0.) && (ints1->Pi[i][j] > 0.)) {
			D_Pi += \
				(double)((ints0->p[i] * ints0->p[j] * i * j) / (4 * ints0->k_av)) * \
				ints0->Pi[i][j] * \
				log(ints0->Pi[i][j] / ints1->Pi[i][j]);

			/* fprintf(stderr, "%s:%d: %d\t%le\t%le\t%le\t%le\n",
				__FILE__, __LINE__,
				i, ints0->p[i], ints0->p[j], ints0->Pi[i][j],
				ints1->Pi[i][j], D_Pi); */
			}
		}
	}

	return D_Pi;
}

/*____________________________________________________________________________*/
/** smoothed mutual information of function Pi */
__inline__ static double miPi_smooth(Ints *ints0, Ints *ints1, int k_max, Arg *arg)
{
	unsigned int i, j;
	double D_Pi = 0.;
	double log();

	for (i = 1; i <= k_max; ++ i) {
		for (j = 1; j <= k_max; ++ j) {
			if ((ints0->p_smooth[i] > 0.) && (ints0->p_smooth[j] > 0.) && \
				(ints0->Pi_smooth[i][j] > 0.) && (ints1->Pi_smooth[i][j] > 0.)) {
			D_Pi += \
				(double)((ints0->p_smooth[i] * ints0->p_smooth[j] * i * j) / (4 * ints0->k_av)) * \
				ints0->Pi_smooth[i][j] * \
				log(ints0->Pi_smooth[i][j] / ints1->Pi_smooth[i][j]);

			/* fprintf(stderr, "%s:%d: %d\t%le\t%le\t%le\t%le\n",
				__FILE__, __LINE__,
				i, ints0->p_smooth[i], ints0->p_smooth[j], ints0->Pi_smooth[i][j],
				ints1->Pi_smooth[i][j], D_Pi); */
			}
		}
	}

	return D_Pi;
}

/*____________________________________________________________________________*/
/** manuscript section 4.2, equation 42 : 
	Kullback-Leibler distance between two networks */
void equation_42(Ints *intsA, Ints *intsB, Dist *dist, Arg *arg, int level)
{
	int i;
	FILE *KLdistanceOutFile = 0;
	char KLdistanceOutFileName[64];
	unsigned int k_max;

	/* initialise values */
	dist->D_A_p = 0.; /* mutual information of function p for network A */
	dist->D_A_Pi = 0.; /* mutual information of function Pi for network A */
	dist->D_B_p = 0.; /* mutual information of function p for network B */
	dist->D_B_Pi = 0.; /* mutual information of function Pi for network B */
	dist->D_AB = 0.; /* Kullback-Leibler distance between two networks A and B */

	/* minimal degree of the two networks to avoid zero probabilities */
	k_max = intsA->k_max < intsB->k_max ? intsA->k_max : intsB->k_max;
	if (! arg->silent)
		fprintf(stdout, "\tk_max for KL-distance set to min(%d, %d) = %d\n",
			intsA->k_max, intsB->k_max, k_max);

	/*____________________________________________________________________________*/
	/** individual computation of the four terms of equation 42 */	
	/** term 1 : mutual information of function p for network A */
	dist->D_A_p = .5 * mip(intsA, intsB, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 2 : mutual information of function Pi for network A */
	dist->D_A_Pi = miPi(intsA, intsB, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 3 : mutual information of function p for network B */
	dist->D_B_p = .5 * mip(intsB, intsA, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 4 : mutual information of function Pi for network B */
	dist->D_B_Pi = miPi(intsB, intsA, k_max, arg);

	/*____________________________________________________________________________*/
	/** total Kullback-Leibler distance between networks A and B */
	dist->D_AB = dist->D_A_p + dist->D_A_Pi + dist->D_B_p + dist->D_B_Pi;
	/*fprintf(stderr, "%s:%d: term1 = %8.2e, term2 = %8.2e, term3 = %8.2e, term4 = %8.2e, total = %8.2e\n",
		__FILE__, __LINE__,
		dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);*/


	/*____________________________________________________________________________*/
	/** output of KL distance for raw Pi matrix */
	/* output to file */
	sprintf(KLdistanceOutFileName, "KL_distance.%c.dat", (level+48));
	KLdistanceOutFile = safe_open(KLdistanceOutFileName, "w");
	fprintf(KLdistanceOutFile, "k_max %d\nD_A_p %8.2e\nD_A_Pi %8.2e\nD_B_p %8.2e\nD_B_Pi %8.2e\nD %8.2e\n",
		k_max, dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);
	fclose(KLdistanceOutFile);

	/* output to stdout */
	if (! arg->silent)
		fprintf(stdout, "\tp-dependent KL distance D_A_p = %8.2e\n"
						"\tp&Pi-dependent KL distance D_A_Pi = %8.2e\n"
						"\tp-dependent KL distance D_B_p = %8.2e\n"
						"\tp&Pi-dependent KL distance D_B_Pi = %8.2e\n"
						"\tKL distance D = %8.2e\n",
			dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);
}


/*____________________________________________________________________________*/
/** complete Kullback-Leibler distance between two networks: rho */
/*  refers to: degree statistics from B, but Pi kernel from A    */
void find_rho(Ints *intsA, Ints *intsB, double *rho, int k_max)
{
	int i,j;
	int n = 0;
	double *rho_new = 0;
	double epsilon = 1.;
	double diff = 0.;
	double oneminus = 0.;
	double val = 0.;

	rho_new = safe_malloc((k_max + 1)* sizeof(double));
    
    /*puts("   preparation ...");*/
    /*for (i = ; i <= k_max; ++ i)
        rho[i] = 1.0;*/
   init_array_double(rho, (k_max + 1), 1.);
 
    /*puts("   mapPing ...");*/
    while (1) {
        ++ n; /* loop counter */
        /* printf("   eps=%1.12lf, n=%2d\n",epsilon,n); */
        diff = 0.0;
		oneminus = 1.0 - epsilon;

        for (i = 1; i <= k_max; ++ i) {
            val = 0.; 
            for (j = 1; j <= k_max; ++ j) {
				val += (intsA->Pi[i][j] * intsB->p[j] * (double)j / (rho[j] * intsB->k_av));
			}
            rho_new[i] = (epsilon * val) + (oneminus * rho[i]);
            val = rho[i] - rho_new[i]; 
            diff += (val * val); 
        }

		/* termination criterion */
        if (diff < 1e-13)
			break;

        for (i = 1; i <= k_max; ++ i)
			rho[i] = rho_new[i];

		/* re-parametrise */
		if (n > 100) {
			n = 0;
			epsilon *= 0.9;
		}
	}
	free(rho_new);
}

/*____________________________________________________________________________*/
/** complete Kullback-Leibler distance between two networks : interference term */
double interference_term(Ints *intsA, Ints *intsB)
{
	int i;
	double D_wiring3 = 0.;
	double D_wiring4 = 0.;
	double *rho;
	unsigned int k_max;

	/* minimal degree of the two networks to avoid zero probabilities */
	k_max = intsA->k_max < intsB->k_max ? intsA->k_max : intsB->k_max;

	/*puts("testing rho calculation ...");*/

	rho = safe_malloc(((intsA->k_max) + 1)* sizeof(double));
    D_wiring3 = 0.0;
	find_rho(intsA, intsA, rho, intsA->k_max);
    for(i = 1; i <= k_max; ++ i)
		D_wiring3 += (intsA->p[i] * (double)i * log(rho[i]) / intsA->k_av);
	free(rho);

	rho = safe_malloc(((intsB->k_max) + 1)* sizeof(double));
	find_rho(intsB, intsB, rho, intsB->k_max);
    for(i = 1; i <= k_max; ++ i)
		D_wiring3 += (intsB->p[i] * (double)i * log(rho[i]) / intsB->k_av);
	free(rho);

    fprintf(stderr, "%s:%d: error level: %lf\n",
		__FILE__, __LINE__, D_wiring3);

	/* calculating first correction */
    /*puts("calculating first rho array ...");*/
	rho = safe_malloc(((intsB->k_max) + 1)* sizeof(double));
    D_wiring3 = 0.0;
	find_rho(intsA, intsB, rho, intsB->k_max);
    for(i = 1; i <= k_max; ++ i)
		D_wiring3 += (intsB->p[i] * (double)i * log(rho[i]));
	free(rho);

	/* calculating second correction */
    /*puts("calculating second rho array ...");*/
	rho = safe_malloc(((intsA->k_max) + 1)* sizeof(double));
    D_wiring4 = 0.0;
	find_rho(intsB, intsA, rho, intsA->k_max);
    for(i = 1; i <= k_max; ++ i)
		D_wiring4 += (intsA->p[i] * (double)i * log(rho[i]));
	free(rho);

    fprintf(stderr, "%s:%d: interference term: %lf\n",
		__FILE__, __LINE__,
		(0.5 * D_wiring3) + (0.5 * D_wiring4)); 

	return ((0.5 * D_wiring3) + (0.5 * D_wiring4));	
}

/*____________________________________________________________________________*/
/** complete Kullback-Leibler distance between two networks */
void KBcompl(Ints *intsA, Ints *intsB, Dist *dist, Arg *arg, int level)
{
	int i;
	FILE *KLdistanceOutFile = 0;
	char KLdistanceOutFileName[64];
	unsigned int k_max;

	/* initialise values */
	dist->D_A_p = 0.; /* mutual information of function p for network A */
	dist->D_A_Pi = 0.; /* mutual information of function Pi for network A */
	dist->D_B_p = 0.; /* mutual information of function p for network B */
	dist->D_B_Pi = 0.; /* mutual information of function Pi for network B */
	dist->D_AB = 0.; /* Kullback-Leibler distance between two networks A and B */
	double interferterm = 0.;  /* interference term for Kullback-Leibler distance */

	/* minimal degree of the two networks to avoid zero probabilities */
	k_max = intsA->k_max < intsB->k_max ? intsA->k_max : intsB->k_max;
	if (! arg->silent)
		fprintf(stdout, "\tk_max for KL-distance set to min(%d, %d) = %d\n",
			intsA->k_max, intsB->k_max, k_max);

	/*____________________________________________________________________________*/
	/** individual computation of the four terms of equation 42 */	
	/** term 1 : mutual information of function p for network A */
	dist->D_A_p = .5 * mip(intsA, intsB, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 2 : mutual information of function Pi for network A */
	dist->D_A_Pi = miPi(intsA, intsB, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 3 : mutual information of function p for network B */
	dist->D_B_p = .5 * mip(intsB, intsA, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 4 : mutual information of function Pi for network B */
	dist->D_B_Pi = miPi(intsB, intsA, k_max, arg);

	/*____________________________________________________________________________*/
	/** completion term for Kullback-Leibler distance */
	interferterm = interference_term(intsA, intsB);

	/*____________________________________________________________________________*/
	/** total Kullback-Leibler distance between networks A and B */
	dist->D_AB = dist->D_A_p + dist->D_A_Pi + dist->D_B_p + dist->D_B_Pi + interferterm;
	/*fprintf(stderr, "%s:%d: term1 = %8.2e, term2 = %8.2e, term3 = %8.2e, term4 = %8.2e, total = %8.2e\n",
		__FILE__, __LINE__,
		dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);*/


	/*____________________________________________________________________________*/
	/** output of KL distance for raw Pi matrix */
	/* output to file */
	sprintf(KLdistanceOutFileName, "KL_distance.%c.dat", (level+48));
	KLdistanceOutFile = safe_open(KLdistanceOutFileName, "w");
	fprintf(KLdistanceOutFile, "k_max %d\nD_A_p %8.2e\nD_A_Pi %8.2e\nD_B_p %8.2e\nD_B_Pi %8.2e\nD %8.2e\n",
		k_max, dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);
	fclose(KLdistanceOutFile);

	/* output to stdout */
	if (! arg->silent)
		fprintf(stdout, "\tp-dependent KL distance D_A_p = %8.2e\n"
						"\tp&Pi-dependent KL distance D_A_Pi = %8.2e\n"
						"\tp-dependent KL distance D_B_p = %8.2e\n"
						"\tp&Pi-dependent KL distance D_B_Pi = %8.2e\n"
						"\tKL distance D = %8.2e\n",
			dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);
}

/*____________________________________________________________________________*/
/* (unused here) Ton's degree complexity function */
void degree_complexity(double *Pk, int *degrees, double kav, int N)  
{
	int i;  
	double Complexity, log(), logfactorial();

    Complexity = 0.0;

	for (i = 0; i < N; i ++) {
			Complexity = Complexity + log(Pk[degrees[i]]) - (double)degrees[i] * \
				log(kav) + kav + log_factorial(degrees[i]);
			/*fprintf(stderr, "%s:%d: i=%d, log(Pk)=%lf, deg=%d, log(kav)=%lf, \
					log_fac(deg)=%lf, Complexity=%lf\n",
					__FILE__, __LINE__,
					i+1, log(Pk[degrees[i]]), degrees[i], log(kav), 
					log_factorial(degrees[i]), Complexity);*/
			}
    Complexity = Complexity / (double)N;
}

/*____________________________________________________________________________*/
/* (unused here) Ton's Pi complexity function */
void Pi_complexity(double **Prob, int *degrees, double kav, int N)
{
	double log(), Complexity, value;
	int i, j;
    Complexity = 0.0; 
    for (i = 0; i < N; i ++)
		for (j = 0; j < N; j ++)
			if((degrees[i] > 0) && (degrees[j] > 0)) {
				value = Prob[degrees[i]][degrees[j]];
				if(value > 0.0) {
					Complexity = Complexity + (double)(degrees[i] * degrees[j]) * value * log(value);
					/*fprintf(stdout, "%s:%d: i=%d, deg_i=%d, deg_j=%d, value=%lf, Complexity=%lf\n",
						__FILE__, __LINE__,
						i+1, degrees[i], degrees[j],  value, Complexity);*/
				}
			}
    Complexity = Complexity / (2.0 * kav * (double)(N * N));
    /*fprintf(stderr, "%s:%d: Pi-complexity: %lf\n",
		__FILE__, __LINE__, Complexity);*/
}

/*____________________________________________________________________________*/
/* mobility */
void mobility(Ints *ints)
{
	unsigned int i, j;
	int sum = 0;
	int sum_square = 0;
	int **source_mat = 0;
	int **target_mat = 0;
	int size = 0;

 /* initialise mobility terms */
	for (i = 0; i < 7; ++ i)
		ints->n[i] = 0.;

	/* allocate matrices */
	source_mat = alloc_mat2D_int(source_mat, ints->N, ints->N);
	target_mat = alloc_mat2D_int(target_mat, ints->N, ints->N);

	/* sum over degree distribution */
	for (i = 0; i <= ints->k_max; ++ i) {
		sum += ints->k[i];
		sum_square += pow(ints->k[i], 2);
	}

	/* term 1 */
	ints->n[1] = .25 * (double)pow(sum, 2);
	
	/* term 2 */
	ints->n[2] = .25 * (double)sum;

	/* term 3 */
	ints->n[3] = .5 * (double)sum_square;

	/* term 4 */
	for (i = 0; i < ints->N; ++ i)
		for (j = 0; j < ints->N; ++ j)
			if (ints->c[i][j] > 0)
				ints->n[4] += ints->k_list[i] * ints->k_list[j];
	ints->n[4] = .5 * (double)ints->n[4];

	size = sizeof(ints->c);
	assert(sizeof(source_mat) == size);
	assert(sizeof(target_mat) == size);
	copy_mat2D_int(target_mat, ints->c, ints->N, ints->N);

	for (i = 0; i < 4; ++ i) {
		/*multiply_mat2D_int(target_mat, source_mat, ints->N, ints->N, source_mat, ints->N, ints->N);*/
		copy_mat2D_int(source_mat, target_mat, ints->N, ints->N);

		if (i == 3)
			ints->n[5] = .25 * (double)trace_mat2D_int(target_mat, ints->N);

		if (i == 2)
			ints->n[6] = .5 * (double)trace_mat2D_int(target_mat, ints->N);
	}

	ints->n[0] = ints->n[1] + ints->n[2] - ints->n[3] - ints->n[4] + ints->n[5] + ints->n[6];

	fprintf(stdout, \
		"\tterm1 %lf\n\tterm2 %lf\n\tterm3 %lf\n\tterm4 %lf\n"
		"\tterm5 %lf\n\tterm6 %lf\n\tmobility %lf\n",
		ints->n[1], ints->n[2], ints->n[3], ints->n[4], \
		ints->n[5], ints->n[6], ints->n[0]);

	free_mat2D_int(source_mat, ints->N);
	free_mat2D_int(target_mat, ints->N);
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
		fprintf(stdout, "\n\tassortativity R = %lf\n", ints->R);
}

/*____________________________________________________________________________*/
/* fill a small constant values to remove zero values */
void fill_zero(Ints *ints, Dist *dist, Arg *arg)
{
    unsigned int i, j;
    double norm = 0.;

    for (i = 1; i <= dist->k_max_pair; ++ i) { 
        /* fill in p array */
        ints->p_dist[i] += arg->fillval;
        norm += ints->p_dist[i];

        for (j = 1; j <= dist->k_max_pair; ++ j) { 
            /* fill in Pi matrix */
            ints->Pi_dist[i][j] += arg->fillval;
        }    
    }    

    for (i = 1; i <= dist->k_max_pair; ++ i) { 
        /* normalise p_dist array */
        ints->p_dist[i] /= norm;
    }    
}

/*____________________________________________________________________________*/
/** compile network properties */
void network_properties(int n, Ints *ints, Arg *arg)
{
	/*____________________________________________________________________________*/
	/* initialise values */
	ints->k_max = 0;
	ints->k_av = 0.;
	ints->P_conn = 0.;

	/*____________________________________________________________________________*/
	/** preparatory routine: compute maximal degree */
	max_degree(ints, arg);

	/*____________________________________________________________________________*/
	/** manuscript section 2.1, equation 1 :
		compute degree frequencies k_i, the overall frequency distribution k[k_i] 
		and probability distribution p[k_i] of degrees */
	equation_1(ints, arg);

	/*____________________________________________________________________________*/
	/** manuscript section 2.1, equation 3 :
		P_conn_k_k1 probability for two randomly drawn nodes with 
		degrees (k,k') to be connected */
	equation_3(ints, arg);

	/*____________________________________________________________________________*/
	/* equations 6 and 8 yield the same results */
	/** manuscript section 2.1, equation 6 : 
		degree correlation function Pi_k_k' */
	equation_6(ints, arg);
	/** manuscript section 2.2, equation 8 : 
		degree correlation function Pi_k_k' */
	/*equation_8(ints, arg);*/

	/*____________________________________________________________________________*/
	/** manuscript section 3.4, equation 35 (and 39) : 
		Shannon entropy */
	equation_35(n, ints, arg); /* entroPic properties per node */
	equation_35_pl(n, ints, arg); /* entroPic properties per link */
	degree_complexity(ints->p, ints->k_list, ints->k_av, ints->N);
	/*Pi_complexity(ints->Pi, ints->k_list, ints->k_av, ints->N);*/
	assortativity(ints, arg);
}

/*____________________________________________________________________________*/
/** compile network pair distance */
void network_pair_distance(Ints *ints0, Ints *ints1, Dist *dist, Arg *arg)
{
	int i = 0;

	/*____________________________________________________________________________*/
	/** manuscript section 4.2, equation 42 : 
		Kullback-Leibler distance between two networks */

	/*
	if (arg->smooth) { 
		smoothen_p(0, ints0, arg);
		smoothen_Pi(0, ints0, arg);
		smoothen_p(1, ints1, arg);
		smoothen_Pi(1, ints1, arg);
	}
	*/

	/*equation_42(ints0, ints1, dist, arg, 0);*/ /* old version */
	KBcompl(ints0, ints1, dist, arg, 0);
}

