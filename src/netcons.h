/*==============================================================================
netcons.h : consistency checks of network properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef NETCONS_H
#define NETCONS_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "error.h"
#include "ints.h"

/*____________________________________________________________________________*/
/* prototypes */
int cons_randGraph_Pi(Ints *ints);
int cons_randGraph_p(Ints *ints);
int cons_randGraph_C(Ints *ints);

#endif
