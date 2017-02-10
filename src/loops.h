#ifndef NETPROP_H
#define NETPROP_H

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
void loops_Tiernan(Ints *ints, Arg *arg);
void print_loops(Ints *ints);


#endif

