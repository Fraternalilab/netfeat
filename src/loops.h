/*==============================================================================
loops.h : network loops 
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/
#ifndef LOOPS_H
#define LOOPS_H

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "arg.h"
#include "array.h"
#include "config.h"
#include "error.h"
#include "ints.h"
#include "matrix.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
void loops_Tiernan(Ints *ints, Arg *arg);
void loops(Ints *ints, Arg *arg);
void print_loops(Ints *ints, Arg *arg);


#endif

