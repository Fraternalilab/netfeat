/*==============================================================================
netfeat : compute network features
(C) 2008-2011 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#include "netfeat.h"

/*____________________________________________________________________________*/
/** netfeat main routine */
int main(int argc, char *argv[])
{
	unsigned int n;

	Ints ints[2]; /* interaction data structures for 1 or 2 networks */
	Ints *pints = 0; /* pointer to various interaction data */
	Ints *pints0 = &(ints[0]);
	Ints *pints1 = &(ints[1]);
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
		pints = &(ints[n]);
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
			read_interaction_matrix(pnet->matIn, pints, &arg, 1);

			/* allocate interaction matrix */
			pints->c = alloc_mat2D_int(pints->c, pints->N, pints->N);
			/* initialise matrix with zero */
			init_mat2D_int(ints->c, ints->N, ints->N, 0);

			/* read interactions */
			read_interaction_matrix(pnet->matIn, pints, &arg, 0);

			if (arg.mtol) {
				if (! arg.silent) fprintf(stdout, "\n\tTransformation of input matrix to output list %s\n",
										pnet->listOut);
				print_interaction_list_ascii(pnet->listOut, pints, &arg);
			}
		/* node list and interaction list are input */
		} else {
			read_node_list(pnet->protList, pints, &arg);
			read_interaction_list(pnet->intsList, pints, &arg);

			/* allocate interaction matrix */
			pints->c = alloc_mat2D_int(pints->c, pints->N, pints->N);
			/* initialise matrix with zero */
			init_mat2D_int(ints->c, ints->N, ints->N, 0);

			assign_interaction_matrix(n, pints, &arg);

			if (arg.ltom) {
				if (! arg.silent) fprintf(stdout, "\n\tTransformation input list to output matrix %s\n",
									pnet->matOut);
				print_interaction_matrix_ascii(pnet->matOut, pints, &arg);
			}
		}

		/*____________________________________________________________________________*/
		if (! arg.noNetProp) {
			if (! arg.silent) fprintf(stdout, "\n\tComputing single network properties\n");

			/*____________________________________________________________________________*/
			/** degree statistics
				k_av : average degree
				k_var : variance of degrees
				k_max : maximal degree
				k[k_i} : degree frequency distribution
				p[k_i] : degree probability distribution */
			degree_statistics(pints, &arg);

			/*____________________________________________________________________________*/
			/** W : probability for two randomly drawn nodes with degrees (k,k') to be connected */
			wiring(pints, &arg);

			/*____________________________________________________________________________*/
			/** Pi[k_i][k_j] : degree-degree correlation function */
			degree_degree_distribution_eq6(pints, &arg);

			/*____________________________________________________________________________*/
			/** Shannon entropy */
			entropy_pn(n, pints, &arg); /* entropic properties per node */
			entropy_pl(n, pints, &arg); /* entropic properties per link */

			/*____________________________________________________________________________*/
			/** assortativity */
			assortativity(pints, &arg);
			network_properties(n, pints, &arg);
			/*network_mults(n, pints, &arg);*/
		}

		/*____________________________________________________________________________*/
		/*
		if (arg.consistency) {
			if (! arg.silent) fprintf(stdout, "\nChecking network consistency\n");
			assert(cons_randGraph_Pi(pints) == 0);
			assert(cons_randGraph_p(pints) == 0);
			assert(cons_randGraph_C(pints) == 0);
		}
		*/
		/*____________________________________________________________________________*/
		/** print output */
		if (! arg.silent) fprintf(stdout, "\nPrinting output\n");
		print_k_list(n, pints);
		print_k(n, pints);
		print_p(n, pints);
		print_P_conn(n, pints);
		print_Pi_corr(n, 0, pints);
		print_etc(n, pints);
		print_entropy(n, pints);
		print_loops(pints, &arg);
	}

	/*____________________________________________________________________________*/
	if (arg.compare) {
		if (! arg.silent) fprintf(stdout, "\nComparing networks\n");

		/*____________________________________________________________________________*/
		/** Kullback-Leibler distance */
		if (! arg.silent) fprintf(stdout, "\nComputing network pair distance\n");
		network_pair_distance(&ints[0], &dist, &arg);
	}

	/*____________________________________________________________________________*/
	/** free memory */
	for (n = 0; n < arg.nNetMax; ++ n) {
		free(arg.net[n].matOut);
		free(arg.net[n].listOut);
	}

	for (n = 0; n < arg.nNet; ++ n) {
		pints = &(ints[n]);

		free_mat2D_int(pints->c, pints->N);
		free(pints->node);
		if (! arg.matInFlag) {
			free(pints->protein1);
			free(pints->protein2);
		}
		if (! arg.noNetProp) {
			free(pints->k);
			free(pints->p);
			free(pints->k_list);
			free(pints->w);
			free_mat2D_double(pints->W, (pints->k_max + 1));
			free_mat2D_double(pints->Pi, (pints->k_max + 1));

			if (arg.compare) {
				free(pints->p_dist);
				free_mat2D_double(pints->Pi_dist, (dist.k_max_pair + 1));
				if (arg.smooth) {
					free(pints->w_smooth);
					free(pints->p_smooth);
					free_mat2D_double(pints->Pi_smooth, (pints->k_max + 1));
					free_mat2D_double(pints->W_smooth, (pints->k_max + 1));
				}
			}
		}
	}

	free(arg.net);
	/*____________________________________________________________________________*/
	if (! arg.silent) fprintf(stdout, "\nClean Termination\n\n");
	return 0;
}
