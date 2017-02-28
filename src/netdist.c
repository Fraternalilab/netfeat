/*==============================================================================
netdist :  network distances
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netdist.h"

/*____________________________________________________________________________*/
/** mutual information of function p */
__inline__ static double mip(Ints *ints0, Ints *ints1, int k_max_pair, Arg *arg)
{
	unsigned int i;
	double D_p = 0.;
	double log();

	for (i = 1; i <= k_max_pair; ++ i) {
		if ((ints0->p_dist[i] > 0.) && (ints1->p_dist[i] > 0.)) {
			D_p += (double)(ints0->p_dist[i] * log(ints0->p_dist[i] / ints1->p_dist[i]));

			/* fprintf(stderr, "%s:%d: %d\t%le\t%le\t%le\n",
			__FILE__, __LINE__,
			i, ints0->p[i], log(ints0->p[i]), ints1->p[i]);*/
		}
	}

	return D_p;
}

/*____________________________________________________________________________*/
/** smoothened mutual information of function p */
__inline__ static double mip_smooth(Ints *ints0, Ints *ints1, int k_max_pair, Arg *arg)
{
	unsigned int i;
	double D_p = 0.;
	double log();

	for (i = 1; i <= k_max_pair; ++ i) {
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
__inline__ static double miPi(Ints *ints0, Ints *ints1, int k_max_pair, Arg *arg)
{
	unsigned int i, j;
	double D_Pi = 0.;
	double log();

	for (i = 1; i <= k_max_pair; ++ i) {
		for (j = 1; j <= k_max_pair; ++ j) {
			if ((ints0->p_dist[i] > 0.) && (ints0->p_dist[j] > 0.) && \
				(ints0->Pi_dist[i][j] > 0.) && (ints1->Pi_dist[i][j] > 0.)) {
				D_Pi += \
				(double)((ints0->p_dist[i] * ints0->p_dist[j] * i * j) / (4 * ints0->k_av)) * \
				ints0->Pi_dist[i][j] * \
				log(ints0->Pi_dist[i][j] / ints1->Pi_dist[i][j]);
			}
		}
	}

	return D_Pi;
}

/*____________________________________________________________________________*/
/** smoothened mutual information of function Pi */
__inline__ static double miPi_smooth(Ints *ints0, Ints *ints1, int k_max_pair, Arg *arg)
{
	unsigned int i, j;
	double D_Pi = 0.;
	double log();

	for (i = 1; i <= k_max_pair; ++ i) {
		for (j = 1; j <= k_max_pair; ++ j) {
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
/** Kullback-Leibler distance between two networks
	(manuscript section 4.2, equation 42) */
void KB(Ints *ints0, Ints *ints1, Dist *dist, Arg *arg)
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
	k_max = ints0->k_max < ints1->k_max ? ints0->k_max : ints1->k_max;
	if (! arg->silent)
		fprintf(stdout, "\tk_max for KL-distance set to min(%d, %d) = %d\n",
			ints0->k_max, ints1->k_max, k_max);

	/*____________________________________________________________________________*/
	/** individual computation of the four terms of equation 42 */	
	/** term 1 : mutual information of function p for network A */
	dist->D_A_p = .5 * mip(ints0, ints1, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 2 : mutual information of function Pi for network A */
	dist->D_A_Pi = miPi(ints0, ints1, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 3 : mutual information of function p for network B */
	dist->D_B_p = .5 * mip(ints1, ints0, k_max, arg);

	/*____________________________________________________________________________*/
	/** term 4 : mutual information of function Pi for network B */
	dist->D_B_Pi = miPi(ints1, ints0, k_max, arg);

	/*____________________________________________________________________________*/
	/** total Kullback-Leibler distance between networks A and B */
	dist->D_AB = dist->D_A_p + dist->D_A_Pi + dist->D_B_p + dist->D_B_Pi;
	/*fprintf(stderr, "%s:%d: term1 = %8.2e, term2 = %8.2e, term3 = %8.2e, term4 = %8.2e, total = %8.2e\n",
		__FILE__, __LINE__,
		dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);*/


	/*____________________________________________________________________________*/
	/** output of KL distance for raw Pi matrix */
	/* output to file */
	sprintf(KLdistanceOutFileName, "KL_distance.dat");
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
void find_rho(Ints *ints0, Ints *ints1, double *rho, int k_max)
{
	int i,j;
	int n = 0;
	int n_all = 0;
	double *rho_new = 0;
	double epsilon = 1.;
	double diff = 0.;
	double oneminus = 0.;
	double val = 0.;

	init_array_double(rho, (k_max + 1), 1.);
	rho_new = safe_malloc((k_max + 1)* sizeof(double));
 
    do {
        ++ n; /* local loop counter */
        ++ n_all; /* global loop counter */
		diff = 0.;
		oneminus = 1. - epsilon;

        for (i = 1; i <= k_max; ++ i) {
            val = 0.; 
            for (j = 1; j <= k_max; ++ j) {
				val += (ints0->Pi_dist[i][j] * ints1->p_dist[j] * (double)j / (rho[j] * ints1->k_av));
				/*if (isnan(val)) {
					fprintf(stderr, "%s:%d: Pi=%e, p=%e, j=%d, rho=%e, k_av=%d\n",
						__FILE__, __LINE__, ints0->Pi_dist[i][j], ints1->p_dist[j], j, rho[j], ints1->k_av);
				}*/
			}
            rho_new[i] = (epsilon * val) + (oneminus * rho[i]);
            val = rho[i] - rho_new[i]; 

            diff += (val * val); 
			/*fprintf(stderr, "%s:%d: diff=%e, eps=%e, rho=%e, rho_new=%e, val=%e, n=%2d\n",
				__FILE__, __LINE__, diff, epsilon, rho[i], rho_new[i], val, n);*/
        }

		if (diff < 1e-13 || isnan(diff)) {
			if (isnan(val))
				Warning("'val' is NaN\n");	
			break;
		}

		copy_arrays_double(rho, rho_new, (k_max + 1));

		/* re-parametrise */
		if (n > 100) {
			n = 0;
			epsilon *= 0.9;
		}
	} while (n_all < 10000);

	free(rho_new);
}

/*____________________________________________________________________________*/
/** complete Kullback-Leibler distance between two networks : interference term */
double interference_term(Ints *ints0, Ints *ints1, Dist *dist)
{
	int i;
	double D_wiring3 = 0.;
	double D_wiring4 = 0.;
	double *rho;

	rho = safe_malloc(((ints0->k_max) + 1)* sizeof(double));

	/*____________________________________________________________________________*/
	/* consistency check */
	/* ints0 */
    D_wiring3 = 0.;
	find_rho(ints0, ints0, rho, dist->k_max_pair);
    for(i = 1; i <= dist->k_min_pair; ++ i)
		D_wiring3 += (ints0->p_dist[i] * (double)i * log(rho[i]) / ints0->k_av);

	/* ints1 */
	find_rho(ints1, ints1, rho, dist->k_max_pair);
    for(i = 1; i <= dist->k_min_pair; ++ i)
		D_wiring3 += (ints1->p_dist[i] * (double)i * log(rho[i]) / ints1->k_av);

    fprintf(stderr, "%s:%d: error level: %lf\n",
		__FILE__, __LINE__, D_wiring3);

	/*____________________________________________________________________________*/
	/* calculating first correction */
    D_wiring3 = 0.;
	find_rho(ints0, ints1, rho, dist->k_max_pair);
    for(i = 1; i <= dist->k_min_pair; ++ i)
		D_wiring3 += (ints1->p_dist[i] * (double)i * log(rho[i]));

	/* calculating second correction */
    D_wiring4 = 0.;
	find_rho(ints1, ints0, rho, dist->k_max_pair);
    for(i = 1; i <= dist->k_min_pair; ++ i)
		D_wiring4 += (ints0->p_dist[i] * (double)i * log(rho[i]));

    fprintf(stderr, "%s:%d: interference term: %lf\n",
		__FILE__, __LINE__,
		(0.5 * D_wiring3) + (0.5 * D_wiring4)); 

	free(rho);

	return ((0.5 * D_wiring3) + (0.5 * D_wiring4));	
}

/*____________________________________________________________________________*/
/** complete Kullback-Leibler distance between two networks */
void KBcompl(Ints *ints0, Ints *ints1, Dist *dist, Arg *arg)
{
	int i;
	FILE *KLdistanceOutFile = 0;
	char KLdistanceOutFileName[64];

	/* initialise values */
	dist->D_A_p = 0.; /* mutual information of function p for network A */
	dist->D_A_Pi = 0.; /* mutual information of function Pi for network A */
	dist->D_B_p = 0.; /* mutual information of function p for network B */
	dist->D_B_Pi = 0.; /* mutual information of function Pi for network B */
	dist->D_AB = 0.; /* Kullback-Leibler distance between two networks A and B */
	double interferterm = 0.;  /* interference term for Kullback-Leibler distance */

	/*____________________________________________________________________________*/
	/** individual computation of the four terms of equation 42 */	
	/** term 1 : mutual information of function p for network A */
	dist->D_A_p = .5 * mip(ints0, ints1, dist->k_max_pair, arg);

	/*____________________________________________________________________________*/
	/** term 2 : mutual information of function Pi for network A */
	dist->D_A_Pi = miPi(ints0, ints1, dist->k_max_pair, arg);

	/*____________________________________________________________________________*/
	/** term 3 : mutual information of function p for network B */
	dist->D_B_p = .5 * mip(ints1, ints0, dist->k_max_pair, arg);

	/*____________________________________________________________________________*/
	/** term 4 : mutual information of function Pi for network B */
	dist->D_B_Pi = miPi(ints1, ints0, dist->k_max_pair, arg);

	/*____________________________________________________________________________*/
	/** completion term for Kullback-Leibler distance */
	interferterm = interference_term(ints0, ints1, dist);
	/*interferterm = 0.;*/

	/*____________________________________________________________________________*/
	/** total Kullback-Leibler distance between networks A and B */
	dist->D_AB = dist->D_A_p + dist->D_A_Pi + dist->D_B_p + dist->D_B_Pi + interferterm;
	/*fprintf(stderr, "%s:%d: term1 = %8.2e, term2 = %8.2e, term3 = %8.2e, term4 = %8.2e, total = %8.2e\n",
		__FILE__, __LINE__,
		dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);*/

	/*____________________________________________________________________________*/
	/** output of KL distance for raw Pi matrix */
	/* output to file */
	sprintf(KLdistanceOutFileName, "KL_distance.dat");
	KLdistanceOutFile = safe_open(KLdistanceOutFileName, "w");
	fprintf(KLdistanceOutFile, "k_max %d\nD_A_p %8.2e\nD_A_Pi %8.2e\nD_B_p %8.2e\nD_B_Pi %8.2e\nD %8.2e\n",
		dist->k_max_pair, dist->D_A_p, dist->D_A_Pi, dist->D_B_p, dist->D_B_Pi, dist->D_AB);
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
/** compute network pair distance */
void network_pair_distance(Ints *ints, Dist *dist, Arg *arg)
{
	unsigned int n;
	Ints *Pints = 0; /* pointer to various interation data */
	Ints *Pints0 = &(ints[0]);
	Ints *Pints1 = &(ints[1]);

	if (! arg->silent) fprintf(stdout, "\nComparing networks\n");

	/*____________________________________________________________________________*/
	/** define min/max and fill values */
	/* max/min of the individual k_max */
	dist->k_max_pair = Pints0->k_max > Pints1->k_max ? Pints0->k_max : Pints1->k_max;
	dist->k_min_pair = Pints0->k_max < Pints1->k_max ? Pints0->k_max : Pints1->k_max;

	/* max/min of the individual N */
	dist->N_max_pair = Pints0->N > Pints1->N ? Pints0->N : Pints1->N;
	dist->N_min_pair = Pints0->N < Pints1->N ? Pints0->N : Pints1->N;

	/* 'fillval' is a constant to fill in '0' matrix values */
	if (arg->fillval)
		/* 'fillval' from command line */
		dist->fillval = arg->fillval;
	else
		/* default 'fillval' adapts to network properties */
		dist->fillval = 1 / (dist->k_max_pair * dist->k_max_pair * sqrt(dist->N_max_pair));

	/*____________________________________________________________________________*/
	/* prepare matrices for network comparison */
	for (n = 0; n < arg->nNet; ++ n) {
		Pints = &(ints[n]);

		/* smoothen if specified */
		if (arg->smooth) {
			smoothen_p(Pints, arg);
			smoothen_Pi(Pints, arg);
			smoothen_W(Pints, arg);
			recompute_w(Pints, dist, arg);
		}

		/*____________________________________________________________________________*/
		/** initialise 'p_dist' and 'Pi_dist' for distance computation; 
			these matrices adopt the dimension of the larger network,
			are filled with a small number 'fillval' to avoid extreme distances
			(finite network effect) and optionally smoothened for the same reason */
		/* 'p_dist' */
		Pints->p_dist = safe_malloc((dist->k_max_pair + 1) * sizeof(double));
		init_array_double(Pints->p_dist, (dist->k_max_pair + 1), dist->fillval);

		/* 'Pi_dist' */
		Pints->Pi_dist = alloc_mat2D_double(Pints->Pi_dist, (dist->k_max_pair + 1), (dist->k_max_pair + 1));
		init_mat2D_double(Pints->Pi_dist, (dist->k_max_pair + 1), (dist->k_max_pair + 1), dist->fillval);

		/** fill matrices for distance computation with values */
		if (arg->smooth) {
			add_array_double(Pints->p_dist, Pints->p_dist, Pints->p_smooth, (Pints->k_max + 1));
			add_mat2D_double(Pints->Pi_dist, Pints->Pi_dist, Pints->Pi_smooth, (Pints->k_max + 1), (Pints->k_max + 1));
		} else {
			add_array_double(Pints->p_dist, Pints->p_dist, Pints->p, (Pints->k_max + 1));
			add_mat2D_double(Pints->Pi_dist, Pints->Pi_dist, Pints->Pi, (Pints->k_max + 1), (Pints->k_max + 1));
		}

		/*____________________________________________________________________________*/
		/** initialise 'w_dist' and 'W_dist' for distance computation; 
			these matrices adopt the dimension of the larger network,
		/* 'w_dist' */
		/*
		Pints->w_dist = safe_malloc((dist->k_max_pair + 1) * sizeof(double));
		init_array_double(Pints->w_dist, (dist->k_max_pair + 1), 0.);
		*/
		/* 'W_dist' */
		/*
		Pints->W_dist = alloc_mat2D_double(Pints->W_dist, (dist->k_max_pair + 1), (dist->k_max_pair + 1));
		init_mat2D_double(Pints->W_dist, (dist->k_max_pair + 1), (dist->k_max_pair + 1), 0.);
		*/
		/** fill matrices for distance computation with values */
		/*
		if (arg->smooth) {
			add_array_double(Pints->w_dist, Pints->w_dist, Pints->w_smooth, (Pints->k_max + 1));
			add_mat2D_double(Pints->W_dist, Pints->W_dist, Pints->Pi_smooth, (Pints->k_max + 1), (Pints->k_max + 1));
		} else {
			add_array_double(Pints->w_dist, Pints->w_dist, Pints->w, (Pints->k_max + 1));
			add_mat2D_double(Pints->W_dist, Pints->W_dist, Pints->W, (Pints->k_max + 1), (Pints->k_max + 1));
		}
		/* normalise W */
		/*
		norm_corr_mat(Pints->W, Pints->W, dist->k_max_pair);
		*/
	}

	/*____________________________________________________________________________*/
	/** Kullback-Leibler distance between two networks
		(manuscript section 4.2, equation 42) */
	KB(Pints0, Pints1, dist, arg); /* p/Pi version */
	/*
	KBcompl(Pints0, Pints1, dist, arg);*/ /* w/W version */
}

