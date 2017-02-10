/*==============================================================================
netcons.c : consistency checks of network properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netcons.h"

/*____________________________________________________________________________*/
/** compute log poisson function:
	log Pi(k) = k*log(<k>) - <k> - SUM_{l=1}^{k}(log(l)) */
__inline__ static double log_poisson(int k, double k_av)
{
	unsigned int i;
	double sum;

	for(i = 2, sum = 0.; i <= k; ++ i)
		sum += log(i);

	return (k * log(k_av) - k_av - sum);
}

/*____________________________________________________________________________*/
/** equality test with rounding error tolerance for double values */
__inline__ static int approximately_equal(double a, double b)
{
	return ((a == b) || (fabsf(a - b) <= .00001F * fabsf(a) + .00001F));
}

/*____________________________________________________________________________*/
/** in regular random graph: Pi[i][j] = 1/(i*j) */
int cons_randGraph_Pi(Ints *ints)
{
	unsigned int i, j;

	for (i = 0; i <= ints->k_max; ++ i) {
		for (j = 0; j <= ints->k_max; ++ j) {
			if (i > 0 && j > 0) {
				if (approximately_equal(ints->Pi[i][j], (1./(i * j))))
					continue;
				else
					fprintf(stderr, "%s:%d: %d %d, %f == %f\n", 
						__FILE__, __LINE__,  i, j,
						ints->Pi[i][j], (1./(i * j)));
			}
		}
	}

	return 0;
}

/*____________________________________________________________________________*/
/** in regular random graph: p(k) = Pi(k) */
int cons_randGraph_p(Ints *ints)
{
	unsigned int i;

	for (i = 0; i <= ints->k_max; ++ i) {
		if (approximately_equal(ints->p[i], log_poisson(i, ints->k_av)))
			continue;
		else
			fprintf(stderr, "%s:%d: %d, %f == %f\n",
				__FILE__, __LINE__,  i, 
				ints->p[i], log_poisson(i, ints->k_av));
	}

	return 0;
}

/*____________________________________________________________________________*/
/** in regular random graph: C(p, Pi) = -1/<k> <log k> */
int cons_randGraph_C(Ints *ints)
{
	if (approximately_equal(ints->C, (-1/ints->k_av * ints->logk_av)))
		;
	else
		fprintf(stderr, "%s:%d: %f == %f\n",
			__FILE__, __LINE__, 
			ints->C, (-1/ints->k_av * ints->logk_av));

	return 0;
}


