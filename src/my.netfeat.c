/*==============================================================================
netfeat : compute network features
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netfeat.h"

/*____________________________________________________________________________*/
/** netfeat main routine */
int main(int argc, char *argv[])
{
	unsigned int n;

	Ints ints[2]; /* interaction data structures for 1 or 2 networks */
	Ints *Pints = 0; /* pointer to various interation data */
	Ints *Pints0 = &(ints[0]);
	Ints *Pints1 = &(ints[1]);
	Dist dist; /* network distance data structure */
	Arg arg; /* command line arguments */
	Net *pnet; /* pointer to network input arguments */

	arg.nNetMax = 2; /* maximal number of networks */

	/*____________________________________________________________________________*/
	/** read command line arguments */
	parse_args(argc, &(argv[0]), &arg);

	/* for each network */
	for (n = 0; n < arg.nNet; ++ n) {
		if (! arg.silent) fprintf(stdout, "\nProcessing network %d\n", n);

		/* point to the currently processed interactions and network */
		Pints = &(ints[n]);
		pnet  = &(arg.net[n]);

		/*____________________________________________________________________________*/
		/* two input data formats are accepted:
			1. interaction matrix 
			2. interaction list and node list
		   see the README for the format specification
		*/
		/* matrix is input */
		if (arg.matInFlag) {

			/* read number of interactions only */
			read_interaction_matrix(pnet->matIn, Pints, &arg, 1);

			/* allocate interaction matrix */
			Pints->c = alloc_mat2D_int(Pints->c, Pints->N, Pints->N);
			/* initialise matrix with zero */
			init_mat2D_int(ints->c, ints->N, ints->N, 0);

			/* read interactions */
			read_interaction_matrix(pnet->matIn, Pints, &arg, 0);

			if (arg.mtol) {
				if (! arg.silent) fprintf(stdout, "\nTransformation matrix->list %s\n",
										pnet->listOut);
				print_interaction_list_ascii(pnet->listOut, &ints);
			}
		/* node list and interaction list are input */
		} else {
			read_node_list(pnet->protList, Pints, &arg);
			read_interaction_list(pnet->intsList, Pints, &arg);

			/* allocate interaction matrix */
			Pints->c = alloc_mat2D_int(Pints->c, Pints->N, Pints->N);
			/* initialise matrix with zero */
			init_mat2D_int(ints->c, ints->N, ints->N, 0);

			assign_interaction_matrix(n, Pints, &arg);

			if (arg.ltom) {
				if (! arg.silent) fprintf(stdout, "\nTransformation list->matrix %s\n",
									pnet->matOut);
				print_interaction_matrix_ascii(pnet->matOut, &ints);
			}
		}

		/*____________________________________________________________________________*/
		if (! arg.noNetProp) {
			if (! arg.silent) fprintf(stdout, "\nComputing network properties\n");
			network_properties(n, Pints, &arg);
			/*network_mults(n, Pints, &arg);*/
		}

		/*____________________________________________________________________________*/
		/*
		if (arg.consistency) {
			if (! arg.silent) fprintf(stdout, "\nChecking network consistency\n");
			assert(cons_randGraph_Pi(Pints) == 0);
			assert(cons_randGraph_p(Pints) == 0);
			assert(cons_randGraph_C(Pints) == 0);
		}
		*/
		/*____________________________________________________________________________*/
		/** print output */
		if (! arg.silent) fprintf(stdout, "\nPrinting output\n");
		print_k_list(n, Pints);
		print_k(n, Pints);
		print_p(n, Pints);
		print_P_conn(n, Pints);
		print_Pi_corr(n, 0, Pints);
		print_etc(n, Pints);
		print_entropy(n, Pints);

		/*____________________________________________________________________________*/

	}

	if (arg.compare) {
		if (! arg.silent) fprintf(stdout, "\nComparing networks\n");

		/* max of the individual k_max */
		dist.k_max_pair = Pints0->k_max > Pints1->k_max ? Pints0->k_max : Pints1->k_max;

		for (n = 0; n < arg.nNet; ++ n) {
		/* prepare matrices for network comparison */
			Pints = &(ints[n]);

			/* initialise 'p_dist' and 'Pi_dist' for distance computation */
			/* 'p_dist' */
			Pints->p_dist = safe_malloc((dist.k_max_pair + 1) * sizeof(double));
			init_array_double(Pints->p_dist, (dist.k_max_pair + 1), 0.);
			copy_array_double(Pints->p_dist, Pints->p, (Pints->k_max + 1));

			/* 'Pi_dist' */
			Pints->Pi_dist = alloc_mat2D_double(Pints->Pi_dist, (dist.k_max_pair + 1), (dist.k_max_pair + 1));
			init_mat2D_double(Pints->Pi_dist, (dist.k_max_pair + 1), (dist.k_max_pair + 1), 0.);
			copy_mat2D_double(Pints->Pi_dist, Pints->Pi, (Pints->k_max + 1), (Pints->k_max + 1));

			/* fill in zero values with a small number */
			fill_zero(Pints, dist, arg);
		}

		/* distance */
		network_pair_distance(Pints0, Pints1, &dist, &arg);
	}

	/*____________________________________________________________________________*/
	/** free memory */
	for (n = 0; n < arg.nNetMax; ++ n) {
		free(arg.net[n].matOut);
		free(arg.net[n].listOut);
	}

	for (n = 0; n < arg.nNet; ++ n) {
		Pints = &(ints[n]);

		free_mat2D_int(Pints->c, Pints->N);
		free(Pints->node);
		if (! arg.matInFlag) {
			free(Pints->protein1);
			free(Pints->protein2);
		}
		if (! arg.noNetProp) {
			free(Pints->k);
			free(Pints->p);
			free(Pints->k_list);
			free_mat2D_int(Pints->conn_k_k1, (Pints->k_max + 1));
			free_mat2D_int(Pints->all_k_k1, (Pints->k_max + 1));
			free_mat2D_double(Pints->P_conn_k_k1, (Pints->k_max + 1));
			free_mat2D_double(Pints->Pi, (Pints->k_max + 1));

			if (arg.compare) {
				free(Pints->p_dist);
				free_mat2D_double(Pints->Pi_dist, (dist.k_max_pair + 1));
				if (arg.smooth) {
					free(Pints->p_smooth);
					free_mat2D_double(Pints->Pi_smooth, (Pints->k_max + 1));
				}
			}
		}
	}

	free(arg.net);
	/*____________________________________________________________________________*/
	if (! arg.silent) fprintf(stdout, "\nClean Termination\n\n");
	return 0;
}
