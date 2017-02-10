/*==============================================================================
netdist.h : network distances
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
==============================================================================*/

#ifndef NETDIST_H
#define NETDIST_H

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "arg.h"
#include "array.h"
#include "error.h"
#include "ints.h"
#include "matrix.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
void network_pair_distance(Ints *ints, Dist *dist, Arg *arg);
void KBcompl(Ints *intsA, Ints *intsB, Dist *dist, Arg *arg);

#endif
