/*==============================================================================
arg.h : parse command line arguments
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "safe.h"

/*____________________________________________________________________________*/
/* structures */
typedef struct
{
	char *matOut; /* name of interaction matrix outputs */
	char *listOut; /* name of interaction list outputs */
	char *intsList; /* name of node interaction lists */
	char *protList; /* name of node lists */
	char *matIn; /* name of interaction matrices (instead of the two lists above) */
	char *piIn; /* name of PI(k,k') matrices for KL-distance computation only */
	char *pIn; /* name of p(k) arrays for KL-distance computation only */
} Net;

typedef struct  
{
	Net *net; /* file names for input/output of network data */
	int matInFlag; /* flag: indicates to use matrix instead of lists */
	int loopDepth; /* loop search depth */
	int selfInteraction; /* flag: count self-interaction of nodes */
	int noNetProp; /* flag: suppress network property computation */
	int compare; /* flag: compare two networks */
	int consistency; /* flag: consistency check */
	int smooth; /* number of smoothing levels */
	double fillval; /* constant to fill in '0' values */
	int mtol; /* transform matrix to list */
	int ltom; /* transform list to matrix */
	int silent; /* flag: suppress stdout */
	int nNetMax; /* maximal number of networks (to intialise) */
	int nNet; /* actual number of networks */
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif
