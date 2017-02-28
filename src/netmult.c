/*==============================================================================
netmult.c : network multiplication
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netmult.h"

/*____________________________________________________________________________*/
/** self-multiplication of the binary interaction matrix */
int network_mults(int n, Ints *ints, Arg *arg)
{
	int i;
	int **mx_c, **my_c; /* multiplied interaction matrices */
	char outFileName[64];

	if (! arg->silent) fprintf(stdout, "\n\tLoop detection in network %d\n", n);

	/* 'mx' indicates multiplied intermediate matrix */
	mx_c = alloc_mat2D_int(mx_c, ints->N, ints->N);
	/* 'my' indicates multiplied final matrix */
	my_c = alloc_mat2D_int(my_c, ints->N, ints->N);

	copy_mat2D_int(mx_c, ints->c, ints->N, ints->N);

	/* create matrices for loops up to exponent i */
	for (i = 0; i < 3; ++ i) {
		sprintf(&(outFileName[0]), "m_c.%d.%d.dat", (i + 2), n);

		multiply_mat2D_int(my_c, mx_c, ints->N, ints->N, ints->c, ints->N, ints->N);
		copy_mat2D_int(mx_c, my_c, ints->N, ints->N);
		regularise_mat2D_int(my_c, ints->N, ints->N);
		print_mat2D_int(&(outFileName[0]), my_c, ints->N, ints->N);
	}

	free_mat2D_int(mx_c, ints->N);
	free_mat2D_int(my_c, ints->N);

	return 0;
}

