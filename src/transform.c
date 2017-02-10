/*==============================================================================
transform.c : transform network matrices
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "transform.h"

/*____________________________________________________________________________*/
/* smoothen Pi */
void smoothen_Pi(Ints *ints, Arg *arg)
{
	int level, k, l, i, j, intw;
	double **SmoothPi = 0;
	double **SmoothPiDat = 0;
	double width, value, norm, exp(), weight;
	int kcut = ints->k_max;
	FILE *outFile = 0;
	char outFileName[64];
	FILE *outFile1 = 0;
	char outFileName1[64];

	if (! arg->silent) fprintf(stdout, "\tSmoothening Pi(k,k')\n");
	ints->Pi_smooth = alloc_mat2D_double(ints->Pi_smooth, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->Pi_smooth, (ints->k_max + 1), (ints->k_max + 1), 0.);

	/* matrix of smoothed data */
	SmoothPi = alloc_mat2D_double(SmoothPi, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(SmoothPi, (ints->k_max + 1), (ints->k_max + 1), 0.);

	for (level = 1; level <= arg->smooth; ++ level) {
		intw = (int)(3 * level + .5);
		width = (double)pow(level, 2);

		for (k = 1; k <= ints->k_max; ++ k) {
			for (l = 1; l <= k; ++ l) {
				value = norm = 0.;
				for (i = k - intw; i <= k + intw; ++ i) {
					for (j = l - intw; j <= l + intw; ++ j){
						if ((i >= 0) && (j >= 0) && (i <= ints->k_max) && (j <= ints->k_max)) {
							weight = (double)(pow((k-i), 2)) + (double)(pow((l-j), 2));
							weight = exp(-0.5 * weight / width);
							value += (weight * (double)ints->Pi[i][j]);
							norm += weight;
						}
					}
					SmoothPi[k][l] = SmoothPi[l][k] = value / norm;
				}
			}
		}
	}

	copy_mat2D_double(ints->Pi_smooth, SmoothPi, (ints->k_max + 1), (ints->k_max + 1));

	free_mat2D_double(SmoothPi, (ints->k_max + 1));
}

/*____________________________________________________________________________*/
/* smoothen p(k) */
void smoothen_p(Ints *ints, Arg *arg)
{
	int level, k, i, intw;
	double *Smooth_p = 0;
	double *Smooth_pDat = 0;
	double width, value, norm, exp(), weight;
	FILE *outFile = 0;
	char outFileName[64];

	if (! arg->silent) fprintf(stdout, "\tSmoothening p(k)\n");
	ints->p_smooth = safe_malloc((ints->k_max + 1) * sizeof(double));
	init_array_double(ints->p_smooth, (ints->k_max + 1), 0.);

	/* matrix of smoothed data */
	Smooth_p = safe_malloc((ints->k_max + 1) * sizeof(double));
	init_array_double(Smooth_p, (ints->k_max + 1), 0.);

	for (level = 1; level <= arg->smooth; ++ level) {
		intw = (int)(3 * level + .5);
		width = (double)pow(level, 2);

		for (k = 0; k <= ints->k_max; ++ k) {
			value = norm = 0.;
			for (i = k - intw; i <= k + intw; ++ i) {
				if ((i >= 0) && (i <= ints->k_max)) {
					weight = (double) pow((k-i), 2);
					weight = exp(-0.5 * weight / width);
					value = value + weight * (double)ints->p[i];
					norm = norm + weight;
				}
			}
			Smooth_p[k] = value / norm;
		}

		/*
		sprintf(outFileName, "p.%c.%c.dat", (n+48), (level+48));
		outFile = safe_open(outFileName, "w");

		print_array_double(outFileName, Smooth_p, (ints->k_max + 1));

		fclose(outFile);
		*/
	}

	copy_array_double(ints->p_smooth, Smooth_p, (ints->k_max + 1)); 

	free(Smooth_p);
	free(Smooth_pDat);
}

/*____________________________________________________________________________*/
/* smoothen W(k,k') */
void smoothen_W(Ints *ints, Arg *arg)
{
	if (! arg->silent) fprintf(stdout, "\tSmoothening W(k,k')\n");
	ints->W_smooth = alloc_mat2D_double(ints->W_smooth, (ints->k_max + 1), (ints->k_max + 1));
	init_mat2D_double(ints->W_smooth, (ints->k_max + 1), (ints->k_max + 1), 0.);

}

/*____________________________________________________________________________*/
/* recompute w(k) after smoothening of W */
void recompute_w(Ints *ints, Dist *dist, Arg *arg)
{
	unsigned int i, j;
	double sum1 = 0.;
	double sum2 = 0.;

	if (! arg->silent) fprintf(stdout, "\tRecomputing w(k)\n");
	ints->w_smooth = safe_malloc((ints->k_max + 1) * sizeof(double));
	init_array_double(ints->w_smooth, (ints->k_max + 1), 0.);

    for (i = 1; i <= dist->k_max_pair; ++ i) {         
		for (j = 1; j <= dist->k_max_pair; ++ j) {
			ints->w_smooth[i] += ints->W_smooth[i][j];
		}
		ints->p[i] = ints->w_smooth[i] * ints->k_av / (double)i;
		sum1 += ints->p[i];
    }

	/* p[0] ??? */
	if ((ints->p[0] = 1. - sum1) < dist->fillval)
		ints->p[0] = dist->fillval;
}

/*____________________________________________________________________________*/
/* normalise correlation matrix */
void norm_corr_mat(double **matNorm, double **matIn, int x)
{
	unsigned int i, j;
	double norm = 0.;

    for (i = 1; i <= x; ++ i)
		for (j = 1; j <= x; ++ j)
			norm += matIn[i][j];

    for (i = 1; i <= x; ++ i)
		for (j = 1; j <= x; ++ j)
			matNorm[i][j] = matIn[i][j] / norm;
}

