/*==============================================================================
netprop.h : macroscopic network properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef NETPROP_H
#define NETPROP_H

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "arg.h"
#include "config.h"
#include "array.h"
#include "error.h"
#include "ints.h"
#include "loops.h"
#include "matrix.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
void max_degree(Ints *ints, Arg *arg);
void average_degree(Ints *ints, Arg *arg);
void network_properties(int n, Ints *ints, Arg *arg);
void mobility(Ints *ints);
void degree_statistics(Ints *ints, Arg *arg);
void wiring(Ints *ints, Arg *arg);
void degree_degree_distribution_eq6(Ints *ints, Arg *arg);
void entropy_pn(int n, Ints *ints, Arg *arg);
void entropy_pl(int n, Ints *ints, Arg *arg);
void assortativity(Ints *ints, Arg *arg);

#endif
