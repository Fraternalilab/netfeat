/*==============================================================================
netmult.h : network multiplication
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef NETMULT_H
#define NETMULT_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "error.h"
#include "ints.h"
#include "matrix.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
int network_mults(int n, Ints *ints, Arg *arg);

#endif
